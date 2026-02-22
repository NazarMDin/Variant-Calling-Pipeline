#!/bin/bash
# Generate accession-to-chr mapping file from reference .fai

FAI="GCF_000001405.40_GRCh38.p14_genomic.fna.fai"
OUT="accession_to_chr.map"

awk '$1 ~ /^NC_0000/ {split($1,a,"."); n=a[1]; sub(/NC_0000/,"",n); n=int(n); if(n>=1 && n<=22) print $1, "chr"n}' $FAI > $OUT
echo "Mapping file created: $OUT"