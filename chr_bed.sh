grep -E '^NC_0000(0[1-9]|1[0-9]|2[0-2])\.' GCF_000001405.40_GRCh38.p14_genomic.fna.fai \
  | cut -f1,2 | awk -v OFS='\t' '{print $1, 0, $2}' > chr1-22.bed