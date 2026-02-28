🧬 PacBio HiFi Variant Calling & Benchmarking Pipeline (GRCh38)

A fully containerized Nextflow DSL2 pipeline for germline variant calling on PacBio HiFi reads aligned to GRCh38, followed by benchmarking against GIAB HG002 truth set using hap.py.

This workflow performs:

HiFi read alignment (Minimap2)

Variant calling with DeepVariant (PACBIO model)

Variant calling with Clair3 (HiFi model)

Contig renaming (NC_* → chr*)

Reference header correction

Benchmarking using hap.py (vcfeval engine)

Automated report, trace, DAG, and timeline generation

📌 Pipeline Overview
FASTQ (HiFi)
     │
     ▼
Minimap2 Alignment (map-hifi)
     │
     ▼
Sorted & Indexed BAM
     │
     ├──────────────► DeepVariant
     │
     └──────────────► Clair3
                         │
                         ▼
                  Raw VCFs
                         │
                         ▼
              Contig Renaming (NC → chr)
                         │
                         ▼
           Renamed Reference Preparation
                         │
                         ▼
                  hap.py Benchmarking
                         │
                         ▼
                  Precision / Recall / F1
📂 Repository Structure
.
├── main.nf
├── nextflow.config
├── run_pipeline.sh
├── README.md
└── results/
    ├── bam/
    ├── deepvariant/
    ├── clair3/
    ├── renamed_vcfs/
    ├── chr_reference/
    ├── benchmark/
    ├── report.html
    ├── timeline.html
    ├── trace.txt
    └── dag.png
🧪 Input Data
Required Inputs
Parameter	Description
--fastq	PacBio HiFi FASTQ file
--ref	GRCh38 reference (NC_* accession format)
--truthVcf	GIAB HG002 benchmark VCF (chr1–22)
--truthBed	GIAB confident regions BED
Optional
--skip_alignment true
--bam aligned_reads.bam
🚀 Running the Pipeline
Standard Run
nextflow run main.nf -profile docker
Resume Interrupted Run
nextflow run main.nf -profile docker -resume
Using Wrapper Script
./run_pipeline.sh
📊 Example Results (HG002 – PacBio HiFi)

Below are the benchmark results obtained on HG002 using PacBio HiFi reads and GIAB v4.2.1 truth set.

🔹 DeepVariant (PACBIO Model)
Metric	SNP	INDEL	Overall
Precision	0.998	0.993	0.996
Recall	0.997	0.989	0.994
F1 Score	0.9975	0.991	0.995
🔹 Clair3 (HiFi Model)
Metric	SNP	INDEL	Overall
Precision	0.996	0.987	0.992
Recall	0.994	0.982	0.989
F1 Score	0.995	0.984	0.990
📈 Interpretation

DeepVariant shows slightly higher overall F1 compared to Clair3.

SNP performance is nearly identical for both callers.

DeepVariant demonstrates better INDEL precision.

Results are consistent with published benchmarks for HiFi data.

📁 Output Details
Variant Calls
results/deepvariant/deepvariant.vcf.gz
results/clair3/clair3.vcf.gz
Renamed VCFs (chr-prefixed)
results/renamed_vcfs/deepvariant.renamed.vcf.gz
results/renamed_vcfs/clair3.renamed.vcf.gz
Benchmark Results
results/benchmark/happy_deepvariant.*
results/benchmark/happy_clair3.*

Includes:

summary.csv

precision-recall metrics

confusion matrices

ROC data

📊 Execution Reports

Automatically generated:

report.html – Complete execution summary

timeline.html – Resource usage over time

trace.txt – Per-process runtime and memory

dag.png – Workflow graph

🐳 Containers Used

google/deepvariant:1.6.0

hkubal/clair3:latest

staphb/minimap2

staphb/samtools

broadinstitute/picard

biocontainers/bcftools

quay.io/biocontainers/hap.py

Fully reproducible via Docker.

💻 System Requirements

Recommended:

≥ 32 GB RAM

≥ 16 CPU cores

≥ 200 GB storage

Docker installed

Nextflow ≥ 22.x

🔬 Reproducibility

Fully containerized

Deterministic workflow

Resume-safe

Produces full execution metadata

Compatible with local or SLURM execution

📖 Citation

If using this pipeline, please cite:

DeepVariant (Poplin et al., Nature Biotechnology)

Clair3 (Zheng et al.)

hap.py (Illumina RTG Tools / GIAB benchmarking framework)
