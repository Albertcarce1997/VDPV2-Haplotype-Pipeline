# VDPV2 NT/AA Haplotype Pipeline

This repository provides a reproducible pipeline for assigning nucleotide (NT) and amino-acid (AA) haplotypes from VDPV2 FASTQ or FASTQ.GZ sequencing data.

The main script is [VDPV2_haplotype_pipeline.py](VDPV2_haplotype_pipeline.py). It processes one or more samples from raw or already quality-filtered reads through dereplication, reference matching, NT haplotype assignment, AA translation, MMseqs2 reference alignment, and final abundance-based haplotype renumbering.

The sample prefix used in output filenames is the input FASTQ basename with `.fastq.gz`, `.fq.gz`, `.fastq`, or `.fq` removed. The prefix is preserved exactly. For example, `sample_A_q20.fastq.gz` produces files beginning with `sample_A_q20`.

## Installation

Create and activate the Conda environment:

```bash
conda env create -f haplotype_env.yml
conda activate pv_haplotypes
```

You can also use `mamba`:

```bash
mamba env create -f haplotype_env.yml
conda activate pv_haplotypes
```

The environment file [haplotype_env.yml](haplotype_env.yml) includes Python, Biopython, `chopper`, `seqtk`, `vsearch`, and MMseqs2.

After activation, check the command-line tools:

```bash
python --version
which chopper seqtk vsearch mmseqs
```

## Required Inputs

The pipeline requires:

- One or more FASTQ/FASTQ.GZ input files.
- A nucleotide reference FASTA for `vsearch`, provided with `--reference`.
- An amino-acid reference FASTA for MMseqs2, provided with `--aa-reference`.

If the reference files are in the same directory as [VDPV2_haplotype_pipeline.py](VDPV2_haplotype_pipeline.py), the default nucleotide reference is `Ref_seq.fasta`, and the default amino-acid reference is `VDPV2_EU_ref.fasta`.

## Basic Usage

Run from the repository directory:

```bash
python VDPV2_haplotype_pipeline.py \
  --input /path/to/sample.fastq.gz \
  --reference /path/to/Ref_seq.fasta \
  --aa-reference /path/to/VDPV2_EU_ref.fasta \
  --output-dir /path/to/results \
  --threads 8
```

If the FASTQ and reference FASTA files are in the repository directory, relative paths are fine:

```bash
python VDPV2_haplotype_pipeline.py \
  --input sample.fastq.gz \
  --reference Ref_seq.fasta \
  --aa-reference VDPV2_EU_ref.fasta \
  --output-dir results \
  --threads 8
```

If `--input` is omitted, all FASTQ/FASTQ.GZ files in the repository directory are processed.

## Optional Q20 Filtering

If your input FASTQ has already been quality filtered, omit `--run-qc`.

To run Q20 filtering before FASTQ-to-FASTA conversion:

```bash
python VDPV2_haplotype_pipeline.py \
  --input /path/to/raw_sample.fastq.gz \
  --run-qc \
  --chopper-quality 20 \
  --chopper-headcrop 30 \
  --chopper-tailcrop 30 \
  --chopper-min-length 700 \
  --reference /path/to/Ref_seq.fasta \
  --aa-reference /path/to/VDPV2_EU_ref.fasta \
  --output-dir /path/to/results \
  --threads 8
```

These are the default chopper settings used whenever `--run-qc` is enabled:

```text
--quality 20
--headcrop 30
--tailcrop 30
--minlength 700
```

The corresponding pipeline options are `--chopper-quality`, `--chopper-headcrop`, `--chopper-tailcrop`, and `--chopper-min-length`. Additional chopper arguments can be passed as one quoted string if needed:

```bash
--chopper-extra-args "--maxlength 1300"
```

## Output Layout

The default output directory is `results`. A custom directory can be provided with `--output-dir`.

```text
<output_dir>/
  qc/
  vsearch/
  nt/
  aa/
```

### `qc/`

Created only when `--run-qc` is used.

```text
<sample>.q20.fastq.gz
```

### `vsearch/`

FASTQ conversion, dereplication, and nucleotide reference-search outputs.

```text
<sample>.fasta
<sample>.derep.fasta
<sample>.aln.txt
<sample>.blast6.tsv
<sample>.match.fasta
<sample>.unmatch.fasta
<sample>.pairs.fasta
```

The `<sample>.pairs.fasta` file contains the pairwise FASTA records used for NT haplotype processing.

### `nt/`

Nucleotide haplotype outputs.

```text
<sample>.nt.clean.fasta
<sample>.nt.gt1.fasta
<sample>.nt.rawHT.fasta
<sample>.nt.curHT.fasta
<sample>.nt.rep.fasta
<sample>.nt.final.fasta
<sample>.nt.report.txt
```

The primary NT outputs are:

```text
<sample>.nt.final.fasta
<sample>.nt.report.txt
```

### `aa/`

Amino-acid haplotype outputs.

```text
<sample>.aa.nt.fasta
<sample>.aa.trans.fasta
<sample>.aa.mmseqs.tsv
<sample>.aa.mmseqs.txt
<sample>.aa.aln.fasta
<sample>.aa.rawHT.fasta
<sample>.aa.rep.fasta
<sample>.aa.min20.fasta
<sample>.aa.final.fasta
<sample>.aa.report.txt
```

If `--min-haplotype-size` is changed, the `.aa.min20.fasta` filename changes accordingly, for example `.aa.min50.fasta`.

The primary AA outputs are:

```text
<sample>.aa.final.fasta
<sample>.aa.report.txt
```

## Pipeline Details

### Read Preparation

When `--run-qc` is enabled, `chopper` filters raw reads with minimum read quality 20, trims 30 bases from the start and end of each read, and retains reads of at least 700 bases. The resulting FASTQ.GZ is written to `qc/`, and the original input sample prefix is retained for all downstream output filenames.

FASTQ records are converted to FASTA with `seqtk seq -a`.

### Dereplication and Reference Search

FASTA reads are dereplicated with `vsearch --fastx_uniques`, preserving `size=` abundance values in FASTA headers.

Dereplicated sequences are searched against the nucleotide reference FASTA with `vsearch --usearch_global` using these defaults:

```text
--id 0.90
--strand both
--top_hits_only
--mincols 600
--maxhits 1
--sizein
```

### NT Haplotype Assignment

The NT workflow starts from `<sample>.pairs.fasta`.

The pipeline removes the reference sequence from each pairwise alignment record, removes alignment-generated gaps, removes records with `size=1`, and assigns NT haplotypes using exact/subsequence matching. Before final NT renumbering, it removes records that translate with a first-frame stop codon and removes haplotypes below `--min-haplotype-size` total abundance.

For each retained NT haplotype, the pipeline sums the `size=` abundance values and keeps one representative sequence. Final NT haplotypes are sorted by abundance and renamed `HT-1`, `HT-2`, `HT-3`, and so on.

### AA Haplotype Assignment

The AA workflow starts from the raw NT haplotype file, `<sample>.nt.rawHT.fasta`, before the final NT stop-codon and abundance filters.

The pipeline removes NT haplotype labels from headers, translates sequences in the first forward frame, keeps only translations without stop codons, and aligns the translated AA sequences to the AA reference with `mmseqs easy-search`.

Default MMseqs2 settings:

```text
--min-seq-id 0.95
--cov-mode 2
-c 0.9
-s 7.5
--max-seqs 1
```

AA haplotypes are assigned from aligned AA sequences using exact/subsequence matching. For each AA haplotype, the pipeline reports the summed abundance and keeps a representative sequence. Haplotypes below `--min-haplotype-size` are removed, and the retained AA haplotypes are renumbered by abundance.

## Common Options

```bash
--threads 16
```

Set the number of threads used by `chopper`, `vsearch`, and MMseqs2.

```bash
--min-haplotype-size 20
```

Set the minimum final NT and AA haplotype abundance. The default is 20.

```bash
--force
```

Regenerate existing outputs.

By default, existing non-empty output files are reused. This allows interrupted or partial runs to resume from the first missing step. For example, if `<sample>.derep.fasta` exists but the vsearch reference-search outputs are missing, the pipeline resumes from the vsearch reference-search step.

```bash
--reference-header <header_to_remove>
```

Set the nucleotide reference header removed from pairwise FASTA records during NT processing.

```bash
--mmseqs-min-seq-id 0.95 --mmseqs-coverage 0.9 --mmseqs-sensitivity 7.5
```

Adjust MMseqs2 alignment thresholds.
