# Variant Calling Pipeline with DeepVariant and Clair3

This repository provides a complete workflow for calling short variants from PacBio HiFi reads aligned to GRCh38, using **DeepVariant** and **Clair3**, and benchmarking against the **GIAB HG002 truth set** (chr1‑22) with **hap.py**.

The pipeline is implemented in **Nextflow** with container support (Singularity + Docker). Due to system‑specific limitations (user namespaces disabled), a fully manual Docker‑based fallback is also provided.

---

## 📦 Requirements

- **Nextflow** (≥21.04) – [install](https://www.nextflow.io/docs/latest/getstarted.html)
- **Singularity** (≥3.5) **or** **Docker** (≥20.10)
- **minimap2** (for alignment) – `conda install -c bioconda minimap2` or `sudo apt install minimap2`
- **samtools** – `conda install -c bioconda samtools` or `sudo apt install samtools`
- Healthy computational power is recommended
- Input files (see below)

---

## 📥 Input Data

You need the following files (adjust paths as needed):

- **FASTQ** – HiFi reads (e.g., `m21009_241011_231051.hifi_reads.fastq.gz`)
- **Reference genome** – GRCh38 (FASTA) from NCBI: `GCF_000001405.40_GRCh38.p14_genomic.fna`
- **Truth set** – GIAB HG002 v4.2.1 for chr1‑22:
  - VCF: `HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz`
  - Index: `HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi`
  - BED: `HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed`
- **Singularity images** (if using Singularity):
  - `deepvariant_1.6.0.sif`
  - `clair3_latest.sif`
  - `happy.sif` (hap.py)

If you lack the Singularity images, the pipeline can use Docker automatically (via `docker://` URIs). The only exception is `hap.py` – we provide a working Docker image in the manual fallback.

---

## 🚀 Full Workflow

### 1. **Prepare Reference Indexes**

```bash
cd /path/to/project
samtools faidx GCF_000001405.40_GRCh38.p14_genomic.fna
samtools dict GCF_000001405.40_GRCh38.p14_genomic.fna -o GCF_000001405.40_GRCh38.p14_genomic.dict
