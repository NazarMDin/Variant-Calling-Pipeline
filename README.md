# Variant Calling Pipeline with DeepVariant and Clair3

This repository provides a complete, automated workflow for calling short variants from PacBio HiFi reads aligned to GRCh38, using **DeepVariant** (PACBIO model) and **Clair3** (HiFi model), and benchmarking against the **GIAB HG002 truth set** (chr1‑22) with **hap.py** (vcfeval engine).

The pipeline is implemented in **Nextflow** and fully containerized (Singularity + Docker). A manual fallback script is also provided for systems where Nextflow cannot be used.

---

## 📥 Input Data

Place the following files in the project root (adjust paths in `main.nf` and `backup_script.sh` if needed):

| File | Description |
|------|-------------|
| `m21009_241011_231051.hifi_reads.fastq.gz` | Raw PacBio HiFi reads (FASTQ) |
| `aligned_reads.bam` (optional) | Pre‑aligned BAM – if you already have it, set `params.skip_alignment = true` in `main.nf` |
| `GCF_000001405.40_GRCh38.p14_genomic.fna` | GRCh38 reference genome (FASTA) |
| `giab_truth/` | Directory containing:<br>• `HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz` + `.tbi`<br>• `HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed` |

*Note:* The pipeline will automatically generate the required `.fai` and `.dict` files for the reference.

---

## 🚀 Quick Start

### 1. Install Nextflow (if not already)

```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin/

### 1. Install Nextflow (if not already)
