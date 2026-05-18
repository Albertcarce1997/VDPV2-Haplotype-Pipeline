#!/usr/bin/env python3
"""VDPV2 nucleotide and amino-acid haplotype pipeline.

This script runs a reproducible end-to-end workflow:

1. Optional chopper Q20 QC.
2. FASTQ to FASTA conversion.
3. vsearch dereplication.
4. vsearch reference search with pairwise FASTA output.
5. NT haplotype curation, filtering, assignment, final filtering, and renumbering.
6. AA no-stop translation, mmseqs2 reference alignment, haplotype assignment,
   representative selection, final filtering, and renumbering.

"""

from __future__ import annotations

import argparse
import gzip
import logging
import re
import shlex
import shutil
import subprocess
import sys
import threading
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator


FASTQ_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")
DEFAULT_REFERENCE_HEADER = "PQ816338.1_Poliovirus_2_isolate_ENV_529_3517_PV2_SPA_091624"
MMSEQS_FORMAT_OUTPUT = (
	"query,target,qaln,taln,qstart,qend,tstart,tend,qlen,tlen,"
	"mismatch,gapopen,fident,nident,pident,qcov,tcov,evalue,bits,cigar,qseq,tseq"
)

CODON_TABLE = {
	"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
	"TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
	"TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
	"TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
	"CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
	"CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
	"CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
	"CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
	"ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
	"ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
	"AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
	"AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
	"GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
	"GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
	"GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
	"GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


@dataclass(frozen=True)
class FastaRecord:
	header: str
	sequence: str


@dataclass(frozen=True)
class StageOutputs:
	sample_name: str
	fasta: Path
	derep_fasta: Path
	blast_results: Path
	pairwise_alignments: Path
	nt_raw_haplotypes: Path
	nt_final_fasta: Path
	nt_report: Path
	aa_final_fasta: Path
	aa_report: Path


def configure_logging(verbose: bool) -> None:
	level = logging.DEBUG if verbose else logging.INFO
	logging.basicConfig(level=level, format="%(asctime)s - %(levelname)s - %(message)s")


def quote_command(command: Iterable[object]) -> str:
	return " ".join(shlex.quote(str(part)) for part in command)


def run_command(command: list[object], stdout_path: Path | None = None) -> None:
	logging.info("Running: %s", quote_command(command))
	try:
		if stdout_path is None:
			subprocess.run([str(part) for part in command], check=True)
		else:
			ensure_directory(stdout_path.parent)
			with stdout_path.open("wb") as stdout_handle:
				subprocess.run([str(part) for part in command], check=True, stdout=stdout_handle)
	except subprocess.CalledProcessError as exc:
		raise SystemExit(f"Command failed with exit code {exc.returncode}: {quote_command(command)}") from exc


def require_executable(name: str) -> None:
	if shutil.which(name) is None:
		raise SystemExit(f"Required executable not found on PATH: {name}")


def require_external_tools(run_qc: bool) -> None:
	for executable in ("seqtk", "vsearch", "mmseqs"):
		require_executable(executable)
	if run_qc:
		require_executable("chopper")


def infer_sample_name(path: Path) -> str:
	name = path.name
	lower_name = name.lower()
	for suffix in FASTQ_SUFFIXES:
		if lower_name.endswith(suffix):
			return name[: -len(suffix)]
	return path.stem


def should_reuse(output_path: Path, force: bool) -> bool:
	if output_path.exists() and output_path.stat().st_size > 0 and not force:
		logging.info("Reusing existing output: %s", output_path)
		return True
	return False


def output_is_ready(output_path: Path) -> bool:
	return output_path.exists() and output_path.stat().st_size > 0


def should_reuse_all(output_paths: Iterable[Path], force: bool) -> bool:
	paths = list(output_paths)
	if paths and all(output_is_ready(path) for path in paths) and not force:
		for path in paths:
			logging.info("Reusing existing output: %s", path)
		return True
	return False


def ensure_directory(path: Path) -> None:
	try:
		path.mkdir(parents=True, exist_ok=True)
	except FileExistsError as exc:
		if path.is_dir():
			return
		try:
			path.mkdir(parents=True, exist_ok=True)
		except FileExistsError:
			pass
		if path.is_dir():
			return
		raise SystemExit(
			f"Cannot create output directory '{path}' because that path already exists but is not a directory. "
			"Remove or rename the conflicting path, or choose a different --output-dir."
		) from exc


def resolve_existing_path(path: str | None, base_dir: Path) -> Path | None:
	if path is None:
		return None
	candidate = Path(path)
	if not candidate.is_absolute():
		candidate = base_dir / candidate
	return candidate.resolve()


def find_default_aa_reference(base_dir: Path) -> Path | None:
	root_reference = base_dir / "VDPV2_EU_ref.fasta"
	if root_reference.exists():
		return root_reference.resolve()
	matches = sorted(base_dir.glob("**/VDPV2_EU_ref.fasta"))
	if len(matches) == 1:
		return matches[0].resolve()
	return None


def iter_fasta(path: Path) -> Iterator[FastaRecord]:
	header: str | None = None
	sequence_parts: list[str] = []
	with path.open("r", encoding="utf-8") as handle:
		for raw_line in handle:
			line = raw_line.strip()
			if not line:
				continue
			if line.startswith(">"):
				if header is not None:
					yield FastaRecord(header=header, sequence="".join(sequence_parts))
				header = line[1:]
				sequence_parts = []
			else:
				sequence_parts.append(line)
	if header is not None:
		yield FastaRecord(header=header, sequence="".join(sequence_parts))


def read_fasta(path: Path) -> list[FastaRecord]:
	return list(iter_fasta(path))


def write_fasta(records: Iterable[FastaRecord], path: Path, line_width: int = 80) -> int:
	ensure_directory(path.parent)
	count = 0
	with path.open("w", encoding="utf-8") as handle:
		for record in records:
			count += 1
			handle.write(f">{record.header}\n")
			sequence = record.sequence
			if line_width <= 0:
				handle.write(f"{sequence}\n")
			else:
				for start in range(0, len(sequence), line_width):
					handle.write(f"{sequence[start:start + line_width]}\n")
	return count


def extract_size(header: str, default: int = 1) -> int:
	match = re.search(r"size=(\d+)", header)
	return int(match.group(1)) if match else default


def extract_haplotype(header: str) -> str | None:
	matches = re.findall(r"HT-\d+", header)
	return matches[-1] if matches else None


def haplotype_number(haplotype: str | None) -> int:
	if haplotype is None:
		return 10**12
	match = re.search(r"\d+", haplotype)
	return int(match.group(0)) if match else 10**12


def replace_terminal_haplotype(header: str, new_haplotype: str) -> str:
	if re.search(r"HT-\d+$", header):
		return re.sub(r"HT-\d+$", new_haplotype, header)
	return re.sub(r"HT-\d+", new_haplotype, header, count=1)


def remove_terminal_haplotype(header: str) -> str:
	return re.sub(r"_HT-\d+$", "", header)


def update_header_size(header: str, new_size: int, haplotype: str | None) -> str:
	if ";size=" not in header:
		suffix = f"_{haplotype}" if haplotype else ""
		return f"{header};size={new_size}{suffix}"

	prefix, suffix = header.split(";size=", 1)
	if "_" in suffix:
		_, rest = suffix.split("_", 1)
		return f"{prefix};size={new_size}_{rest}"
	suffix_ht = f"_{haplotype}" if haplotype else ""
	return f"{prefix};size={new_size}{suffix_ht}"


def translate_nt_first_frame(sequence: str) -> str:
	normalized = sequence.upper().replace("U", "T").replace("-", "")
	remainder = len(normalized) % 3
	if remainder:
		normalized += "N" * (3 - remainder)
	amino_acids: list[str] = []
	for start in range(0, len(normalized), 3):
		codon = normalized[start:start + 3]
		amino_acids.append(CODON_TABLE.get(codon, "X"))
	return "".join(amino_acids).rstrip("X")


def contains_first_frame_stop(sequence: str) -> bool:
	return "*" in translate_nt_first_frame(sequence)


def count_records_and_reads(path: Path) -> tuple[int, int]:
	records = 0
	reads = 0
	for record in iter_fasta(path):
		records += 1
		reads += extract_size(record.header)
	return records, reads


def fastq_has_complete_record(path: Path) -> bool:
	opener = gzip.open if path.name.lower().endswith(".gz") else open
	try:
		with opener(path, "rb") as handle:
			lines = [handle.readline() for _ in range(4)]
	except (EOFError, OSError):
		return False
	return all(lines) and lines[0].startswith(b"@") and lines[2].startswith(b"+")


def run_chopper_qc(input_fastq: Path, output_fastq: Path, args: argparse.Namespace) -> Path:
	if output_fastq.exists() and not args.force:
		if output_fastq.stat().st_size > 0 and fastq_has_complete_record(output_fastq):
			logging.info("Reusing existing output: %s", output_fastq)
			return output_fastq
		logging.warning("Existing QC output is incomplete and will be recreated: %s", output_fastq)

	ensure_directory(output_fastq.parent)
	command = [
		"chopper",
		"--quality", str(args.chopper_quality),
		"--headcrop", str(args.chopper_headcrop),
		"--tailcrop", str(args.chopper_tailcrop),
		"--minlength", str(args.chopper_min_length),
		"--threads", str(args.threads),
	]
	if args.chopper_max_length is not None:
		command.extend(["--maxlength", str(args.chopper_max_length)])
	if args.chopper_extra_args:
		command.extend(shlex.split(args.chopper_extra_args))

	logging.info("Running chopper QC: %s", quote_command(command))
	input_opener = gzip.open if input_fastq.name.lower().endswith(".gz") else open
	output_opener = gzip.open if output_fastq.name.lower().endswith(".gz") else open

	with input_opener(input_fastq, "rb") as input_handle, output_opener(output_fastq, "wb") as output_handle:
		process = subprocess.Popen([str(part) for part in command], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
		assert process.stdin is not None
		assert process.stdout is not None

		def write_input() -> None:
			try:
				for chunk in iter(lambda: input_handle.read(1024 * 1024), b""):
					process.stdin.write(chunk)
			except BrokenPipeError:
				pass
			finally:
				try:
					process.stdin.close()
				except BrokenPipeError:
					pass

		writer = threading.Thread(target=write_input)
		writer.start()
		for chunk in iter(lambda: process.stdout.read(1024 * 1024), b""):
			output_handle.write(chunk)
		return_code = process.wait()
		writer.join()

	if return_code != 0:
		try:
			output_fastq.unlink()
		except FileNotFoundError:
			pass
		raise SystemExit(f"chopper failed with exit code {return_code} for {input_fastq}")

	return output_fastq


def convert_fastq_to_fasta(input_fastq: Path, output_fasta: Path, force: bool) -> Path:
	if should_reuse(output_fasta, force):
		return output_fasta
	run_command(["seqtk", "seq", "-a", input_fastq], stdout_path=output_fasta)
	return output_fasta


def dereplicate_fasta(input_fasta: Path, output_fasta: Path, threads: int, force: bool) -> Path:
	if should_reuse(output_fasta, force):
		return output_fasta
	run_command([
		"vsearch",
		"--fastx_uniques", input_fasta,
		"--fastaout", output_fasta,
		"--sizein",
		"--sizeout",
		"--strand", "both",
		"--threads", str(threads),
	])
	return output_fasta


def vsearch_output_paths(output_dir: Path, sample_name: str) -> tuple[Path, Path, Path, Path, Path]:
	alignment_output = output_dir / f"{sample_name}.aln.txt"
	blast_results = output_dir / f"{sample_name}.blast6.tsv"
	matched_fasta = output_dir / f"{sample_name}.match.fasta"
	notmatched_fasta = output_dir / f"{sample_name}.unmatch.fasta"
	pairwise_alignments = output_dir / f"{sample_name}.pairs.fasta"
	return alignment_output, blast_results, matched_fasta, notmatched_fasta, pairwise_alignments


def vsearch_reference_search(
	derep_fasta: Path,
	reference_fasta: Path,
	output_dir: Path,
	sample_name: str,
	threads: int,
	force: bool,
) -> tuple[Path, Path, Path, Path, Path]:
	alignment_output, blast_results, matched_fasta, notmatched_fasta, pairwise_alignments = vsearch_output_paths(output_dir, sample_name)

	expected = [alignment_output, blast_results, matched_fasta, pairwise_alignments]
	if should_reuse_all(expected, force):
		return alignment_output, blast_results, matched_fasta, notmatched_fasta, pairwise_alignments

	run_command([
		"vsearch",
		"--usearch_global", derep_fasta,
		"--db", reference_fasta,
		"--id", "0.90",
		"--strand", "both",
		"--alnout", alignment_output,
		"--blast6out", blast_results,
		"--notmatched", notmatched_fasta,
		"--threads", str(threads),
		"--fastapairs", pairwise_alignments,
		"--top_hits_only",
		"--mincols", "600",
		"--maxhits", "1",
		"--matched", matched_fasta,
		"--sizein",
	])
	return alignment_output, blast_results, matched_fasta, notmatched_fasta, pairwise_alignments


def load_reference_headers(reference_fasta: Path, configured_header: str | None = None) -> set[str]:
	headers = {record.header for record in iter_fasta(reference_fasta)}
	headers.update(header.split()[0] for header in list(headers))
	if configured_header:
		headers.add(configured_header)
		headers.add(configured_header.split()[0])
	return headers


def curate_pairwise_alignments(
	pairwise_fasta: Path,
	output_fasta: Path,
	reference_headers: set[str],
	force: bool,
) -> Path:
	if should_reuse(output_fasta, force):
		return output_fasta

	curated_records: list[FastaRecord] = []
	removed_references = 0
	for record in iter_fasta(pairwise_fasta):
		header_key = record.header.split()[0]
		if record.header in reference_headers or header_key in reference_headers:
			removed_references += 1
			continue
		curated_records.append(FastaRecord(record.header, record.sequence.replace("-", "")))

	written = write_fasta(curated_records, output_fasta)
	logging.info("Curated NT pairwise alignments: wrote %s records, removed %s reference records", written, removed_references)
	return output_fasta


def filter_size_one(input_fasta: Path, output_fasta: Path, force: bool) -> Path:
	if should_reuse(output_fasta, force):
		return output_fasta

	kept: list[FastaRecord] = []
	total = 0
	removed = 0
	for record in iter_fasta(input_fasta):
		total += 1
		if extract_size(record.header) == 1:
			removed += 1
			continue
		kept.append(record)

	write_fasta(kept, output_fasta)
	logging.info("Filtered size=1 NT records: kept %s of %s records", total - removed, total)
	return output_fasta


def assign_haplotypes_by_subsequence(records: list[FastaRecord]) -> dict[str, str]:
	sorted_records = sorted(records, key=lambda record: len(record.sequence), reverse=True)
	assignments: dict[str, str] = {}
	current_haplotype = 1

	for index, record in enumerate(sorted_records):
		if index and index % 1000 == 0:
			logging.info("Assigned %s haplotypes after scanning %s/%s records", current_haplotype - 1, index, len(sorted_records))

		if record.header in assignments:
			continue

		haplotype = f"HT-{current_haplotype}"
		assignments[record.header] = haplotype

		for compare_record in sorted_records[index + 1:]:
			if compare_record.header in assignments:
				continue
			if (
				record.sequence == compare_record.sequence
				or compare_record.sequence in record.sequence
				or record.sequence in compare_record.sequence
			):
				assignments[compare_record.header] = haplotype

		current_haplotype += 1

	return assignments


def assign_nt_haplotypes(input_fasta: Path, output_fasta: Path, force: bool) -> Path:
	if should_reuse(output_fasta, force):
		return output_fasta

	records = read_fasta(input_fasta)
	assignments = assign_haplotypes_by_subsequence(records)
	output_records = [FastaRecord(f"{record.header}_{assignments[record.header]}", record.sequence) for record in records]
	write_fasta(output_records, output_fasta)
	logging.info("Assigned NT haplotypes: %s input records, %s haplotypes", len(records), len(set(assignments.values())))
	return output_fasta


def filter_nt_haplotypes_for_final(input_fasta: Path, output_fasta: Path, min_size: int, force: bool) -> Path:
	if should_reuse(output_fasta, force):
		return output_fasta

	records = read_fasta(input_fasta)
	no_stop_records: list[FastaRecord] = []
	removed_stop_records = 0
	removed_missing_ht = 0

	for record in records:
		if extract_haplotype(record.header) is None:
			removed_missing_ht += 1
			continue
		if contains_first_frame_stop(record.sequence):
			removed_stop_records += 1
			continue
		no_stop_records.append(record)

	haplotype_counts: dict[str, int] = defaultdict(int)
	for record in no_stop_records:
		haplotype = extract_haplotype(record.header)
		if haplotype is not None:
			haplotype_counts[haplotype] += extract_size(record.header)

	keep_haplotypes = {haplotype for haplotype, count in haplotype_counts.items() if count >= min_size}
	curated = [record for record in no_stop_records if extract_haplotype(record.header) in keep_haplotypes]
	write_fasta(curated, output_fasta)

	removed_low_size_records = len(no_stop_records) - len(curated)
	logging.info(
		"Curated NT haplotypes: kept %s records, removed %s stop-containing records, %s low-abundance records, %s records without HT",
		len(curated), removed_stop_records, removed_low_size_records, removed_missing_ht,
	)
	return output_fasta


def build_haplotype_mapping(records: list[FastaRecord]) -> dict[str, dict[str, float | int | str]]:
	haplotype_counts: dict[str, int] = defaultdict(int)
	haplotype_order: list[str] = []
	for record in records:
		haplotype = extract_haplotype(record.header)
		if haplotype is None:
			continue
		if haplotype not in haplotype_counts:
			haplotype_order.append(haplotype)
		haplotype_counts[haplotype] += extract_size(record.header)

	total = sum(haplotype_counts.values())
	sorted_haplotypes = sorted(haplotype_order, key=lambda haplotype: haplotype_counts[haplotype], reverse=True)

	mapping: dict[str, dict[str, float | int | str]] = {}
	for index, old_haplotype in enumerate(sorted_haplotypes, start=1):
		count = haplotype_counts[old_haplotype]
		mapping[old_haplotype] = {
			"new_name": f"HT-{index}",
			"count": count,
			"frequency": (count / total * 100) if total else 0.0,
		}
	return mapping


def renumber_haplotypes(input_fasta: Path, output_fasta: Path, report_file: Path, force: bool) -> Path:
	if should_reuse_all([output_fasta, report_file], force):
		return output_fasta

	records = read_fasta(input_fasta)
	mapping = build_haplotype_mapping(records)
	renumbered_records: list[FastaRecord] = []

	for record in records:
		old_haplotype = extract_haplotype(record.header)
		if old_haplotype is None or old_haplotype not in mapping:
			continue
		new_haplotype = str(mapping[old_haplotype]["new_name"])
		renumbered_records.append(FastaRecord(replace_terminal_haplotype(record.header, new_haplotype), record.sequence))

	renumbered_records.sort(key=lambda record: haplotype_number(extract_haplotype(record.header)))
	write_fasta(renumbered_records, output_fasta)
	write_haplotype_report_from_fasta(output_fasta, report_file)
	logging.info("Renumbered %s haplotypes into %s", len(mapping), output_fasta)
	return output_fasta


def write_haplotype_report_from_fasta(input_fasta: Path, report_file: Path) -> None:
	records = read_fasta(input_fasta)
	haplotype_counts: dict[str, int] = defaultdict(int)
	representative: dict[str, FastaRecord] = {}

	for record in records:
		haplotype = extract_haplotype(record.header)
		if haplotype is None:
			continue
		haplotype_counts[haplotype] += extract_size(record.header)
		if haplotype not in representative or len(record.sequence) > len(representative[haplotype].sequence):
			representative[haplotype] = record

	total = sum(haplotype_counts.values())
	ensure_directory(report_file.parent)
	with report_file.open("w", encoding="utf-8") as handle:
		handle.write("Haplotype Report\n")
		handle.write("================\n\n")
		handle.write(f"Source FASTA\t{input_fasta}\n")
		handle.write(f"Total haplotypes\t{len(haplotype_counts)}\n")
		handle.write(f"Total sequences\t{total}\n\n")
		handle.write("Haplotype\tSequence Count\tFrequency (%)\tRepresentative Length\tRepresentative Header\n")
		for haplotype in sorted(haplotype_counts, key=haplotype_number):
			count = haplotype_counts[haplotype]
			frequency = (count / total * 100) if total else 0.0
			rep = representative[haplotype]
			handle.write(f"{haplotype}\t{count}\t{frequency:.2f}\t{len(rep.sequence)}\t{rep.header}\n")


def strip_nt_haplotype_headers(input_fasta: Path, output_fasta: Path, force: bool) -> Path:
	if should_reuse(output_fasta, force):
		return output_fasta
	records = [FastaRecord(remove_terminal_haplotype(record.header), record.sequence) for record in iter_fasta(input_fasta)]
	write_fasta(records, output_fasta)
	logging.info("Removed NT haplotype suffixes from %s records", len(records))
	return output_fasta


def translate_no_stops(input_fasta: Path, output_fasta: Path, force: bool) -> Path:
	if should_reuse(output_fasta, force):
		return output_fasta

	translated_records: list[FastaRecord] = []
	removed_stops = 0
	for record in iter_fasta(input_fasta):
		amino_acids = translate_nt_first_frame(record.sequence)
		if "*" in amino_acids:
			removed_stops += 1
			continue
		translated_records.append(FastaRecord(record.header, amino_acids))
	write_fasta(translated_records, output_fasta)
	logging.info("Translated AA no-stop records: kept %s, removed %s", len(translated_records), removed_stops)
	return output_fasta


def copy_fasta(input_fasta: Path, output_fasta: Path, force: bool) -> Path:
	if should_reuse(output_fasta, force):
		return output_fasta
	ensure_directory(output_fasta.parent)
	shutil.copyfile(input_fasta, output_fasta)
	return output_fasta


def mmseqs_align_to_reference(
	query_fasta: Path,
	reference_fasta: Path,
	result_tsv: Path,
	aligned_fasta: Path,
	temp_dir: Path,
	threads: int,
	force: bool,
	min_seq_id: float,
	coverage: float,
	sensitivity: float,
) -> tuple[Path, Path]:
	ensure_directory(result_tsv.parent)
	ensure_directory(aligned_fasta.parent)

	if should_reuse_all([result_tsv, aligned_fasta], force):
		return result_tsv, aligned_fasta

	ensure_directory(temp_dir)
	run_command([
		"mmseqs", "easy-search",
		query_fasta,
		reference_fasta,
		result_tsv,
		temp_dir,
		"--format-output", MMSEQS_FORMAT_OUTPUT,
		"-s", str(sensitivity),
		"--alignment-mode", "0",
		"--min-seq-id", str(min_seq_id),
		"--cov-mode", "2",
		"-c", str(coverage),
		"--threads", str(threads),
		"--max-seqs", "1",
	])

	aligned_records: list[FastaRecord] = []
	if result_tsv.exists() and result_tsv.stat().st_size > 0:
		with result_tsv.open("r", encoding="utf-8") as handle:
			for line in handle:
				fields = line.rstrip("\n").split("\t")
				if len(fields) < 3:
					continue
				query_header = fields[0]
				query_alignment = fields[2].replace("-", "")
				aligned_records.append(FastaRecord(query_header, query_alignment))

	write_fasta(aligned_records, aligned_fasta)
	shutil.rmtree(temp_dir, ignore_errors=True)
	logging.info("MMseqs alignment retained %s aligned AA records", len(aligned_records))
	return result_tsv, aligned_fasta


def assign_aa_haplotypes(input_fasta: Path, output_fasta: Path, force: bool) -> Path:
	if should_reuse(output_fasta, force):
		return output_fasta

	records = [record for record in iter_fasta(input_fasta) if len(record.sequence) >= 15]
	sorted_records = sorted(records, key=lambda record: len(record.sequence), reverse=True)
	assignments: dict[str, str] = {}
	output_records: list[FastaRecord] = []
	next_haplotype = 1

	for index, record in enumerate(sorted_records):
		assigned = False
		for previous_record in sorted_records[:index]:
			if previous_record.header not in assignments:
				continue
			if (
				record.sequence == previous_record.sequence
				or record.sequence in previous_record.sequence
				or previous_record.sequence in record.sequence
			):
				assignments[record.header] = assignments[previous_record.header]
				assigned = True
				break
		if not assigned:
			assignments[record.header] = f"HT-{next_haplotype}"
			next_haplotype += 1
		output_records.append(FastaRecord(f"{record.header}_{assignments[record.header]}", record.sequence))

	write_fasta(output_records, output_fasta)
	logging.info("Assigned AA haplotypes: %s records, %s haplotypes", len(output_records), len(set(assignments.values())))
	return output_fasta


def merge_haplotype_reads(input_fasta: Path, output_fasta: Path, force: bool) -> Path:
	if should_reuse(output_fasta, force):
		return output_fasta

	grouped: dict[str, dict[str, object]] = {}
	order: list[str] = []

	for record in iter_fasta(input_fasta):
		haplotype = extract_haplotype(record.header)
		if haplotype is None:
			continue
		if haplotype not in grouped:
			grouped[haplotype] = {"total_size": 0, "representative": record}
			order.append(haplotype)
		grouped[haplotype]["total_size"] = int(grouped[haplotype]["total_size"]) + extract_size(record.header)
		representative = grouped[haplotype]["representative"]
		if isinstance(representative, FastaRecord) and len(record.sequence) > len(representative.sequence):
			grouped[haplotype]["representative"] = record

	output_records: list[FastaRecord] = []
	for haplotype in sorted(order, key=haplotype_number):
		representative = grouped[haplotype]["representative"]
		if not isinstance(representative, FastaRecord):
			continue
		total_size = int(grouped[haplotype]["total_size"])
		output_records.append(FastaRecord(update_header_size(representative.header, total_size, haplotype), representative.sequence))

	write_fasta(output_records, output_fasta)
	logging.info("Collapsed haplotypes into %s representative records", len(output_records))
	return output_fasta


def filter_min_size(input_fasta: Path, output_fasta: Path, min_size: int, force: bool) -> Path:
	if should_reuse(output_fasta, force):
		return output_fasta
	records = read_fasta(input_fasta)
	filtered = [record for record in records if extract_size(record.header) >= min_size]
	write_fasta(filtered, output_fasta)
	logging.info("Filtered records by min size %s: kept %s of %s", min_size, len(filtered), len(records))
	return output_fasta


def write_mmseqs_summary(query_fasta: Path, aligned_fasta: Path, summary_file: Path) -> None:
	initial_clusters, initial_reads = count_records_and_reads(query_fasta)
	aligned_clusters, aligned_reads = count_records_and_reads(aligned_fasta)
	cluster_diff = ((initial_clusters - aligned_clusters) / initial_clusters * 100) if initial_clusters else 100.0
	read_diff = ((initial_reads - aligned_reads) / initial_reads * 100) if initial_reads else 100.0

	ensure_directory(summary_file.parent)
	with summary_file.open("w", encoding="utf-8") as handle:
		handle.write("Metric\tValue\n")
		handle.write(f"Initial clusters\t{initial_clusters}\n")
		handle.write(f"Aligned clusters\t{aligned_clusters}\n")
		handle.write(f"Cluster difference (%)\t{cluster_diff:.2f}\n")
		handle.write(f"Initial sequences\t{initial_reads}\n")
		handle.write(f"Aligned sequences\t{aligned_reads}\n")
		handle.write(f"Sequence difference (%)\t{read_diff:.2f}\n")


def process_nt_haplotypes(
	pairwise_alignments: Path,
	reference_fasta: Path,
	output_root: Path,
	sample_name: str,
	args: argparse.Namespace,
) -> tuple[Path, Path, Path]:
	nt_dir = output_root / "nt"
	ensure_directory(nt_dir)

	curated_pairwise = nt_dir / f"{sample_name}.nt.clean.fasta"
	size_filtered = nt_dir / f"{sample_name}.nt.gt1.fasta"
	raw_haplotypes = nt_dir / f"{sample_name}.nt.rawHT.fasta"
	curated_haplotypes = nt_dir / f"{sample_name}.nt.curHT.fasta"
	representative_haplotypes = nt_dir / f"{sample_name}.nt.rep.fasta"
	final_haplotypes = nt_dir / f"{sample_name}.nt.final.fasta"
	report_file = nt_dir / f"{sample_name}.nt.report.txt"

	if should_reuse_all([raw_haplotypes, final_haplotypes, report_file], args.force):
		logging.info("Completed NT outputs already exist for %s", sample_name)
		return raw_haplotypes, final_haplotypes, report_file

	reference_headers = load_reference_headers(reference_fasta, args.reference_header)
	curate_pairwise_alignments(pairwise_alignments, curated_pairwise, reference_headers, args.force)
	filter_size_one(curated_pairwise, size_filtered, args.force)
	assign_nt_haplotypes(size_filtered, raw_haplotypes, args.force)
	filter_nt_haplotypes_for_final(raw_haplotypes, curated_haplotypes, args.min_haplotype_size, args.force)
	merge_haplotype_reads(curated_haplotypes, representative_haplotypes, args.force)
	renumber_haplotypes(representative_haplotypes, final_haplotypes, report_file, args.force)

	logging.info("Completed NT haplotypes for %s", sample_name)
	return raw_haplotypes, final_haplotypes, report_file


def process_aa_haplotypes(
	nt_raw_haplotypes: Path,
	aa_reference_fasta: Path,
	output_root: Path,
	sample_name: str,
	args: argparse.Namespace,
) -> tuple[Path, Path]:
	aa_dir = output_root / "aa"
	ensure_directory(aa_dir)

	nt_without_ht = aa_dir / f"{sample_name}.aa.nt.fasta"
	translated_sequences = aa_dir / f"{sample_name}.aa.trans.fasta"
	mmseqs_tsv = aa_dir / f"{sample_name}.aa.mmseqs.tsv"
	aligned_fasta = aa_dir / f"{sample_name}.aa.aln.fasta"
	aligned_ht = aa_dir / f"{sample_name}.aa.rawHT.fasta"
	representative_ht = aa_dir / f"{sample_name}.aa.rep.fasta"
	min_size_ht = aa_dir / f"{sample_name}.aa.min{args.min_haplotype_size}.fasta"
	final_haplotypes = aa_dir / f"{sample_name}.aa.final.fasta"
	report_file = aa_dir / f"{sample_name}.aa.report.txt"
	summary_file = aa_dir / f"{sample_name}.aa.mmseqs.txt"
	temp_dir = aa_dir / f"{sample_name}.aa.tmp"

	if should_reuse_all([final_haplotypes, report_file], args.force):
		logging.info("Completed AA outputs already exist for %s", sample_name)
		return final_haplotypes, report_file

	strip_nt_haplotype_headers(nt_raw_haplotypes, nt_without_ht, args.force)
	translate_no_stops(nt_without_ht, translated_sequences, args.force)
	_, aligned_fasta = mmseqs_align_to_reference(
		translated_sequences,
		aa_reference_fasta,
		mmseqs_tsv,
		aligned_fasta,
		temp_dir,
		args.threads,
		args.force,
		args.mmseqs_min_seq_id,
		args.mmseqs_coverage,
		args.mmseqs_sensitivity,
	)
	assign_aa_haplotypes(aligned_fasta, aligned_ht, args.force)
	write_mmseqs_summary(translated_sequences, aligned_fasta, summary_file)
	merge_haplotype_reads(aligned_ht, representative_ht, args.force)
	filter_min_size(representative_ht, min_size_ht, args.min_haplotype_size, args.force)
	renumber_haplotypes(min_size_ht, final_haplotypes, report_file, args.force)

	logging.info("Completed AA haplotypes for %s", sample_name)
	return final_haplotypes, report_file


def process_sample(input_fastq: Path, reference_fasta: Path, aa_reference_fasta: Path, args: argparse.Namespace) -> StageOutputs:
	output_root = args.output_dir.resolve()
	qc_dir = output_root / "qc"
	vsearch_dir = output_root / "vsearch"
	ensure_directory(vsearch_dir)

	sample_name = infer_sample_name(input_fastq)
	working_fastq = input_fastq.resolve()
	if args.run_qc:
		qc_output_name = args.qc_output_name or f"{sample_name}.q20.fastq.gz"
		working_fastq = run_chopper_qc(input_fastq, qc_dir / qc_output_name, args)

	fasta = vsearch_dir / f"{sample_name}.fasta"
	derep_fasta = vsearch_dir / f"{sample_name}.derep.fasta"
	alignment_output, blast_results, matched_fasta, _, pairwise_alignments = vsearch_output_paths(vsearch_dir, sample_name)

	if should_reuse_all([alignment_output, blast_results, matched_fasta, pairwise_alignments], args.force):
		logging.info("Completed vsearch reference-search outputs already exist for %s", sample_name)
	else:
		if should_reuse(derep_fasta, args.force):
			logging.info("Dereplication already exists for %s; resuming at vsearch reference search", sample_name)
		else:
			convert_fastq_to_fasta(working_fastq, fasta, args.force)
			dereplicate_fasta(fasta, derep_fasta, args.threads, args.force)
		_, blast_results, _, _, pairwise_alignments = vsearch_reference_search(
			derep_fasta,
			reference_fasta,
			vsearch_dir,
			sample_name,
			args.threads,
			args.force,
		)

	if not pairwise_alignments.exists() or pairwise_alignments.stat().st_size == 0:
		raise SystemExit(f"Pairwise alignment FASTA is missing or empty for {sample_name}: {pairwise_alignments}")

	nt_raw_haplotypes, nt_final_fasta, nt_report = process_nt_haplotypes(
		pairwise_alignments,
		reference_fasta,
		output_root,
		sample_name,
		args,
	)
	aa_final_fasta, aa_report = process_aa_haplotypes(
		nt_raw_haplotypes,
		aa_reference_fasta,
		output_root,
		sample_name,
		args,
	)

	return StageOutputs(
		sample_name=sample_name,
		fasta=fasta,
		derep_fasta=derep_fasta,
		blast_results=blast_results,
		pairwise_alignments=pairwise_alignments,
		nt_raw_haplotypes=nt_raw_haplotypes,
		nt_final_fasta=nt_final_fasta,
		nt_report=nt_report,
		aa_final_fasta=aa_final_fasta,
		aa_report=aa_report,
	)


def collect_input_fastqs(base_dir: Path, input_paths: list[str]) -> list[Path]:
	if input_paths:
		paths = []
		for input_path in input_paths:
			candidate = Path(input_path)
			if not candidate.is_absolute():
				candidate = base_dir / candidate
			if not candidate.exists():
				raise SystemExit(f"Input FASTQ not found: {candidate}")
			paths.append(candidate.resolve())
		return paths

	discovered = [path.resolve() for path in base_dir.iterdir() if path.is_file() and path.name.lower().endswith(FASTQ_SUFFIXES)]
	if not discovered:
		raise SystemExit("No FASTQ inputs provided and none found in the pipeline directory.")
	return sorted(discovered)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
	script_dir = Path(__file__).resolve().parent
	parser = argparse.ArgumentParser(
		description="Run the VDPV2 NT/AA haplotype pipeline from FASTQ or FASTQ.GZ input.",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
	)
	parser.add_argument("--input", action="append", default=[], help="FASTQ/FASTQ.GZ input. Repeat for multiple samples. If omitted, FASTQ files in the script directory are used.")
	parser.add_argument("--output-dir", type=Path, default=script_dir / "results", help="Directory for all generated outputs.")
	parser.add_argument("--reference", default="Ref_seq.fasta", help="Nucleotide reference FASTA for vsearch.")
	parser.add_argument("--aa-reference", default=None, help="Amino-acid reference FASTA for mmseqs2. Defaults to VDPV2_EU_ref.fasta if found.")
	parser.add_argument("--reference-header", default=DEFAULT_REFERENCE_HEADER, help="Reference header to remove from vsearch pairwise FASTA during NT curation.")
	parser.add_argument("--threads", type=int, default=8, help="Threads for vsearch, chopper, and mmseqs2.")
	parser.add_argument("--min-haplotype-size", type=int, default=20, help="Minimum final haplotype abundance retained for NT and AA outputs.")
	parser.add_argument("--force", action="store_true", help="Overwrite/recreate outputs instead of reusing existing non-empty files.")
	parser.add_argument("--verbose", action="store_true", help="Print debug logging.")

	qc_group = parser.add_argument_group("optional chopper QC")
	qc_group.add_argument("--run-qc", action="store_true", help="Run chopper before FASTQ to FASTA conversion.")
	qc_group.add_argument("--chopper-quality", type=int, default=20, help="Minimum read quality retained by chopper.")
	qc_group.add_argument("--chopper-headcrop", type=int, default=30, help="Bases trimmed from the start of each read by chopper.")
	qc_group.add_argument("--chopper-tailcrop", type=int, default=30, help="Bases trimmed from the end of each read by chopper.")
	qc_group.add_argument("--chopper-min-length", "--chopper-minlength", dest="chopper_min_length", type=int, default=700, help="Minimum read length retained by chopper.")
	qc_group.add_argument("--chopper-max-length", "--chopper-maxlength", dest="chopper_max_length", type=int, default=None, help="Optional maximum read length retained by chopper.")
	qc_group.add_argument("--chopper-extra-args", default="", help="Additional quoted arguments passed to chopper.")
	qc_group.add_argument("--qc-output-name", default=None, help="Optional filename for the QC FASTQ.GZ output when processing one input.")

	mmseqs_group = parser.add_argument_group("mmseqs2 alignment")
	mmseqs_group.add_argument("--mmseqs-min-seq-id", type=float, default=0.95, help="MMseqs2 --min-seq-id value.")
	mmseqs_group.add_argument("--mmseqs-coverage", type=float, default=0.9, help="MMseqs2 -c query coverage value with --cov-mode 2.")
	mmseqs_group.add_argument("--mmseqs-sensitivity", type=float, default=7.5, help="MMseqs2 -s sensitivity value.")

	args = parser.parse_args(argv)
	args.script_dir = script_dir
	if args.qc_output_name and len(args.input) > 1:
		raise SystemExit("--qc-output-name can only be used with a single --input sample.")
	return args


def main(argv: list[str] | None = None) -> int:
	args = parse_args(argv)
	configure_logging(args.verbose)

	base_dir = args.script_dir
	if not args.output_dir.is_absolute():
		args.output_dir = base_dir / args.output_dir

	reference_fasta = resolve_existing_path(args.reference, base_dir)
	if reference_fasta is None or not reference_fasta.exists():
		raise SystemExit(f"Nucleotide reference FASTA not found: {args.reference}")

	aa_reference_fasta = resolve_existing_path(args.aa_reference, base_dir) if args.aa_reference else find_default_aa_reference(base_dir)
	if aa_reference_fasta is None or not aa_reference_fasta.exists():
		raise SystemExit(
			"Amino-acid reference FASTA not found. Put VDPV2_EU_ref.fasta in the pipeline directory "
			"or pass --aa-reference /path/to/VDPV2_EU_ref.fasta."
		)

	require_external_tools(args.run_qc)
	input_fastqs = collect_input_fastqs(base_dir, args.input)
	if args.qc_output_name and len(input_fastqs) > 1:
		raise SystemExit("--qc-output-name can only be used when exactly one input FASTQ is processed.")
	logging.info("Processing %s input FASTQ file(s)", len(input_fastqs))

	outputs: list[StageOutputs] = []
	for input_fastq in input_fastqs:
		logging.info("Starting sample from input: %s", input_fastq)
		outputs.append(process_sample(input_fastq, reference_fasta, aa_reference_fasta, args))

	logging.info("Pipeline complete")
	for output in outputs:
		logging.info("Sample: %s", output.sample_name)
		logging.info("  vsearch pairwise FASTA: %s", output.pairwise_alignments)
		logging.info("  NT final FASTA: %s", output.nt_final_fasta)
		logging.info("  NT report: %s", output.nt_report)
		logging.info("  AA final FASTA: %s", output.aa_final_fasta)
		logging.info("  AA report: %s", output.aa_report)
	return 0


if __name__ == "__main__":
	sys.exit(main())
