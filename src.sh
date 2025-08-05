#!/bin/bash

# Input VCF (bgzipped and tabix-indexed)
VCF_FILE="GCF_000001405.40.gz"

# Chromosome name mapping file (e.g., NC_000001.11 -> chr1)
RENAME_FILE="RefSeq.seq.accession2UCSC.style.name.txt"

# UCSC-style chromosome names to process

# Create a renamed temporary VCF (on-the-fly pipe)
# echo "🔄 Renaming chromosomes in VCF..."
RENAMED_VCF_TMP="GCF_000001405.40_UCSC.gz"

bcftools annotate \
  --rename-chrs "$RENAME_FILE" \
  --threads 6 \
   "$VCF_FILE" -Oz -o "$RENAMED_VCF_TMP"

 # Index the renamed VCF
 tabix -p vcf "$RENAMED_VCF_TMP"

CHROMOSOMES=(chr{1..22} chrX)
RENAMED_VCF_TMP="GCF_000001405.40_UCSC.gz"
# Output directory
OUT_DIR="dbsnp_by_chr"
mkdir -p "$OUT_DIR"

# # Loop through chromosomes
 for CHR in "${CHROMOSOMES[@]}"; do
   echo "🔍 Processing $CHR..."
   bcftools view \
    -v snps \
    -i 'ID ~ "^rs"' \
    -r "$CHR" \
    --threads 8 \
    "$RENAMED_VCF_TMP" -Ou | \
   bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' | \
   awk -F'\t' '{
    split($4, alts, ",");
    for (i in alts) {
      print $1 "\t" $2 "\t" $3 "\t" alts[i] "\t" $5
    }
   }' > "${OUT_DIR}/${CHR}_dbsnp.tsv"
done

echo "✅ Done! Output written to: $OUT_DIR/"


INPUT_DIR="./dbsnp_by_chr"
OUTPUT_DIR="./cpra_rsID"
mkdir -p "$OUTPUT_DIR"

echo "🔄 Converting dbSNP TSVs to CHR:POS:REF:ALT → rsID format..."

for infile in "$INPUT_DIR"/chr*_dbsnp.tsv; do
  chr_name=$(basename "$infile")
  outfile="$OUTPUT_DIR/$chr_name"

  echo "  → Processing $chr_name"
  awk '{
    gsub(/^chr/, "", $1);   # remove "chr" prefix
    print $1":"$2":"$3":"$4 "\t" $5
  }' "$infile" > "$outfile"
done

echo "✅ Done! Converted files are saved in: $OUTPUT_DIR/"


