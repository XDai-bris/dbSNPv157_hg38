# dbSNP Build 157 (GRCh38) Processing

This repository provides scripts and instructions to download, verify, and process **dbSNP Build 157 VCF data** for *Homo sapiens* (GRCh38 assembly), producing per-chromosome TSV files with variants in `CHR:POS:REF:ALT → rsID` format.

---

## 📥 1. Download dbSNP VCF File

Download the compressed dbSNP VCF file from the NCBI FTP server:

- [`GCF_000001405.40.gz`](https://ftp.ncbi.nih.gov/snp/archive/b157/VCF/GCF_000001405.40.gz)
- [`GCF_000001405.40.gz.tbi`](https://ftp.ncbi.nih.gov/snp/archive/b157/VCF/GCF_000001405.40.gz.tbi)

Place both files in your working directory.

---

## 🧪 2. Verify File Integrity

Run an MD5 checksum to verify file integrity:

```bash
md5sum GCF_000001405.40.gz
```

Compare the result against:

```
6a6f313e92a39c337571174dad12cfe1  GCF_000001405.40.gz
ba10bcbae4f0ad9b01244efdd564d6e2  GCF_000001405.40.gz.tbi
```

---

## 🔄 3. Chromosome Name Mapping

The reference VCF uses RefSeq-style contig names (e.g., `NC_000001.11`). To convert these to UCSC-style (e.g., `chr1`), you’ll need a renaming file.

### Steps:

1. Download the sequence report from:
   [NCBI GCF_000001405.40](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/)

2. Convert to rename format in R:

```r
data <- read.table("./sequence_report.tsv", sep = '\t', header = TRUE)
out <- cbind(data$RefSeq.seq.accession, data$UCSC.style.name)
write.table(out, file = 'RefSeq.seq.accession2UCSC.style.name.txt', quote = FALSE,
            sep = " ", row.names = FALSE, col.names = FALSE)
```

---

## ⚙️ 4. Main Script: Rename, Extract, and Format dbSNP Variants

This section describes the full process of transforming the raw dbSNP VCF file into per-chromosome variant tables in a simplified format (`CHR:POS:REF:ALT → rsID`), suitable for lookups or downstream pipelines.

### ✅ What this script does:
1. **Renames chromosomes** in the dbSNP VCF file from RefSeq-style (e.g., `NC_000001.11`) to UCSC-style (e.g., `chr1`).
2. **Splits the VCF** into individual chromosomes (chr1–chr22, chrX).
3. **Extracts biallelic SNPs** with valid RSIDs (`rs*`) and exports them to clean TSVs.
4. **Converts the output** into `CHR:POS:REF:ALT` format for fast lookups.

---

### 📜 Script Overview

Save the following as `process_dbsnp.sh` and run it from your working directory:

```bash
#!/bin/bash
```

---

### 🔧 Step-by-Step Breakdown

#### 1. **Input and Rename Setup**

```bash
# Input VCF (bgzipped and tabix-indexed)
VCF_FILE="GCF_000001405.40.gz"

# Rename mapping file (RefSeq → UCSC-style names)
RENAME_FILE="RefSeq.seq.accession2UCSC.style.name.txt"

# Output renamed VCF path
RENAMED_VCF_TMP="GCF_000001405.40_UCSC.gz"
```

#### 2. **Rename Chromosomes Using bcftools**

```bash
bcftools annotate \
  --rename-chrs "$RENAME_FILE" \
  --threads 6 \
  "$VCF_FILE" -Oz -o "$RENAMED_VCF_TMP"

# Index the new renamed VCF
tabix -p vcf "$RENAMED_VCF_TMP"
```

📝 *This converts contig names like `NC_000001.11` to `chr1`, ensuring compatibility with UCSC-style tools.*

---

#### 3. **Define Target Chromosomes**

```bash
CHROMOSOMES=(chr{1..22} chrX)
OUT_DIR="dbsnp_by_chr"
mkdir -p "$OUT_DIR"
```

---

#### 4. **Extract Variants Per Chromosome**

```bash
for CHR in "${CHROMOSOMES[@]}"; do
  echo "🔍 Processing $CHR..."

  bcftools view \
    -v snps \                                # Only SNPs
    -i 'ID ~ "^rs"' \                        # Only RSIDs (skip . or non-rs)
    -r "$CHR" \                              # Filter by chromosome
    --threads 8 \
    "$RENAMED_VCF_TMP" -Ou | \

  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' | \

  # Handle multi-allelic sites → split into multiple lines
  awk -F'\t' '{
    split($4, alts, ",");
    for (i in alts) {
      print $1 "\t" $2 "\t" $3 "\t" alts[i] "\t" $5
    }
  }' > "${OUT_DIR}/${CHR}_dbsnp.tsv"
done
```

📂 Output:  
Each chromosome's data is saved in `dbsnp_by_chr/chr*_dbsnp.tsv`, with the format:
```
CHROM    POS    REF    ALT    RSID
```

---

#### 5. **Reformat to CHR:POS:REF:ALT → rsID Format**

```bash
INPUT_DIR="./dbsnp_by_chr"
OUTPUT_DIR="./cpra_rsID"
mkdir -p "$OUTPUT_DIR"
```

```bash
echo "🔄 Converting TSVs to CHR:POS:REF:ALT → rsID format..."

for infile in "$INPUT_DIR"/chr*_dbsnp.tsv; do
  chr_name=$(basename "$infile")
  outfile="$OUTPUT_DIR/$chr_name"

  echo "  → Processing $chr_name"

  awk '{
    gsub(/^chr/, "", $1);   # remove "chr" prefix
    print $1":"$2":"$3":"$4 "\t" $5
  }' "$infile" > "$outfile"
done
```

📂 Output:  
Formatted files like `cpra_rsID/chr1_dbsnp.tsv` will contain:
```
1:55516888:A:G    rs123456
```

---

### ✅ Summary

This script:
- Converts dbSNP to UCSC-style naming
- Filters for valid RSIDs
- Breaks up variants by chromosome
- Outputs clean TSVs ready for fast lookups

Ideal for use in:
- Variant annotation pipelines
- Indexing RSIDs for external datasets
- Cross-referencing with functional genomics files

## 📊 5. dbSNP Build 157 Summary

### Organism

- **Species**: *Homo sapiens*
- **Taxonomy ID**: 9606

### RefSNP Summary

- **Total RS IDs**: 1,206,053,617

| Chromosome | RS Count |   | Chromosome | RS Count |
|------------|----------|---|------------|----------|
| chr1       | 93,459,666 | | chr13      | 40,024,876 |
| chr2       | 98,025,545 | | chr14      | 35,970,751 |
| chr3       | 80,072,046 | | chr15      | 33,788,289 |
| chr4       | 77,039,720 | | chr16      | 36,346,514 |
| chr5       | 72,149,666 | | chr17      | 33,478,570 |
| chr6       | 67,614,181 | | chr18      | 31,556,825 |
| chr7       | 64,880,744 | | chr19      | 25,440,948 |
| chr8       | 60,512,771 | | chr20      | 25,938,627 |
| chr9       | 51,144,408 | | chr21      | 15,323,376 |
| chr10      | 54,273,273 | | chr22      | 16,138,508 |
| chr11      | 55,865,576 | | chrX       | 45,803,018 |
| chr12      | 54,080,082 | | chrY       | 3,048,013 |

**Other regions**:
- chrM: 9,229
- Alt Only: 35,320
- Patch: 4,679
- PAR: 1,427,972
- Not On: 88,222
- Unplaced: 279,151
- Unlocalized: 296,229

### RS Status:

- Live RS: 1,172,689,405  
- Unsupported RS: 4,922,330  
- Withdrawn RS: 6,854,924  
- Merged RS: 21,585,672  
- Locationless RS: 1,286  

---

## 📁 6. FTP Site Information

- FTP URL: [`ftp.ncbi.nih.gov/snp/archive/b157/VCF/`](https://ftp.ncbi.nih.gov/snp/archive/b157/VCF/)
- Data Source: NCBI dbSNP Build 157  
- Official Release Notes: Found in `/latest_release/release_notes.txt`

---

## 🧠 Notes on File Formats and Usage

- **VCF is bgzipped + tabix-indexed** for fast querying.
- **Regions (-R)** for `bcftools view` can use:
  - `.vcf`, `.tab`, or `.bed` formats
  - For `.tab`: use `CHROM [TAB] POS`
  - For `.bed`: use 0-based half-open intervals: `CHROM [TAB] START [TAB] END`

> See [GRC Assembly Definitions](https://www.ncbi.nlm.nih.gov/grc/help/definitions/) for ALT, PAR, and other region definitions.

---

## ✅ Final Output

You will get:
- `dbsnp_by_chr/chr*_dbsnp.tsv`: Variant rows in tab-separated format
- `cpra_rsID/chr*_dbsnp.tsv`: Converted format `CHR:POS:REF:ALT → rsID`

These files are ready for downstream variant lookups or annotation workflows.
