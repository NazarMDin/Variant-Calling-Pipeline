# Variant Calling Pipeline with DeepVariant and Clair3

This repository contains a Nextflow pipeline for calling short variants from PacBio HiFi reads using **DeepVariant** and **Clair3**, followed by benchmarking against the GIAB HG002 truth set (GRCh38, chr1-22) using **hap.py**.

The pipeline is containerized with **Singularity** (primary) and **Docker** (fallback). Due to system limitations (user namespaces disabled), Singularity failed; therefore, manual Docker commands are provided as a working alternative.

---

## Requirements

- **Nextflow** (>=21.04) – [installation guide](https://www.nextflow.io/docs/latest/getstarted.html)
- **Singularity** (>=3.5) OR **Docker** (>=20.10)
- Healthy computation power recommended 
- Input files (see below)

---

## Input Data

Place the following files in the project root:

- `aligned_reads.bam` + `aligned_reads.bam.bai` – HiFi reads aligned to GRCh38
- `GCF_000001405.40_GRCh38.p14_genomic.fna` – Reference genome
- `GCF_000001405.40_GRCh38.p14_genomic.fna.fai` – Reference index
- `GCF_000001405.40_GRCh38.p14_genomic.dict` – Reference dictionary
- `giab_truth/` – Directory containing:
  - `HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` + `.tbi`
  - `HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed`
- `singularity_images/` – Directory with `.sif` files (DeepVariant, Clair3, hap.py) – *optional if using Docker*

---

## Pipeline Overview

The Nextflow pipeline (`main.nf`) performs the following steps:

1. **DeepVariant** – calls variants using `make_examples`, `call_variants`, `postprocess_variants`
2. **Clair3** – calls variants with HiFi model
3. **Contig Renaming** – converts accession-based contigs (e.g., `NC_000001.11`) to `chr1` format using `bcftools`
4. **Benchmarking** – runs hap.py for each caller against GIAB truth

Results are saved in the `results/` directory.

---

## Usage

### 1. Prepare Input Files

Ensure all required files are present (see above). Then generate a BED file for autosomes:

```bash
awk -v OFS='\t' '{if($1 ~ /^NC_0000[0-9]+/ && $1 !~ /_/) print $1, 0, $2}' GCF_000001405.40_GRCh38.p14_genomic.fna.fai | head -22 > chr1-22.bed
