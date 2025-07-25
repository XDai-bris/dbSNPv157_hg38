# dbSNP Build 157 (GRCh38) Processing

This repository provides scripts and instructions to download, verify, and process dbSNP Build 157 VCF data for **Homo sapiens** (GRCh38 assembly), producing a BED file of variant positions.

---

## 1. Download dbSNP VCF File

Download the compressed dbSNP VCF file from the NCBI FTP server:

- [`GCF_000001405.40.gz`](https://ftp.ncbi.nih.gov/snp/archive/b157/VCF/GCF_000001405.40.gz)
- [`GCF_000001405.40.gz.tbi`](https://ftp.ncbi.nih.gov/snp/archive/b157/VCF/GCF_000001405.40.gz.tbi)

Place both files in your working directory.

---

## 2. Verify File Integrity

Run an MD5 checksum to verify the integrity of the downloaded file:

```bash
md5sum GCF_000001405.40.gz
```

Compare the result against the checksum provided from the FTP server.

       6a6f313e92a39c337571174dad12cfe1  GCF_000001405.40.gz
       ba10bcbae4f0ad9b01244efdd564d6e2  GCF_000001405.40.gz.tbi

---

## 3. Extract Variant only Biallelic with Valid RSid
```bash
bcftools view \
  -v snps \
  --min-alleles 2 --max-alleles 2 \
  -i 'ID ~ "^rs"' \
  GCF_000001405.40.gz -Oz -o dbsnp157_biallelic_rs_hg38.vcf.gz

```

## 4. Rename the CHROM from contigs to CHR*
to check the CHROM format in the reference vcf file
```bash
 bcftools query -f '%CHROM\n' GCF_000001405.40.gz |uniq
 bcftools query -f '%CHROM\n' GCF_000001405.40.gz |uniq|wc -l
```
found https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/ to download sequencing report (sequence_report.tsv) which contains 'UCSC.style.name', 'RefSeq.seq.accession'
following https://samtools.github.io/bcftools/bcftools.html and https://www.youtube.com/watch?v=LetWDG54hvg to make the rename.txt file which origin Seqname in the first colmun and new one in the second
sep=' ', seperate by whitespace
```r
data <- read.table("./sequence_report.tsv", sep = '\t', header = T)
out <- cbind(data$RefSeq.seq.accession, data$UCSC.style.name)
write.table(out, file = 'RefSeq.seq.accession2UCSC.style.name.txt', quote = F,
            sep = " ", row.names = F, col.names = F)
```
```bash
 bcftools annotate --rename-chrs RefSeq.seq.accession2UCSC.style.name.txt dbsnp157_biallelic_rs_hg38.vcf.gz -Oz -o dbsnp157_biallelic_rs_hg38_UCSC.vcf.gz

 bcftools query -f '%CHROM\n' dbsnp157_biallelic_rs_hg38_UCSC.vcf.gz|uniq
 
 bcftools query -f '%CHROM\t%POS\n' dbsnp157_biallelic_rs_hg38_UCSC.vcf.gz > dbsnp157_biallelic_rs_hg38_UCSC.tab 
```

# futher information when using bcftool -R with .bed and .tab file 
       Regions can be specified either on command line or in a VCF, BED, or tab-delimited file (the default). The columns of the tab-delimited file can contain either positions (two-column format: CHROM, POS) or intervals (three-column format: CHROM, BEG, END), but not both. Positions are 1-based and inclusive. The columns of the tab-delimited BED file are also CHROM, POS and END (trailing columns are ignored), but coordinates are 0-based, half-open. To indicate that a file be treated as BED rather than the 1-based tab-delimited file, the file must have the ".bed" or ".bed.gz" suffix (case-insensitive). Uncompressed files are stored in memory, while bgzip-compressed and tabix-indexed region files are streamed. Note that sequence names must match exactly, "chr20" is not the same as "20". Also note that chromosome ordering in FILE will be respected, the VCF will be processed in the order in which chromosomes first appear in FILE. However, within chromosomes, the VCF will always be processed in ascending genomic coordinate order no matter what order they appear in FILE. Note that overlapping regions in FILE can result in duplicated out of order positions in the output. This option requires indexed VCF/BCF files. Note that -R cannot be used in combination with -r.

## 5. dbSNP Build 157 Release Notes Summary

### Organism

- **Species**: Homo sapiens
- **Taxonomy ID**: 9606

### RefSNP (RS) Summary

- **Total RS IDs**: 1,206,053,617

#### RS Counts by Chromosome:

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
| chrM       | 9,229     |   | Alt Only   | 35,320    |
| Patch      | 4,679     |   | PAR        | 1,427,972 |
| Not On     | 88,222    |   | Unplaced   | 279,151   |
| Unlocalized| 296,229   |   |            |           |

> See [GRC Assembly Definitions](https://www.ncbi.nlm.nih.gov/grc/help/definitions/) for information on ALT, PAR, and other regions.

#### RS Counts by Status:

- **Live RS**: 1,172,689,405  
- **Unsupported RS**: 4,922,330  
- **Withdrawn RS**: 6,854,924  
- **Locationless RS**: 1,286  
- **Merged RS**: 21,585,672

---

### SubSNP (SS) Summary

- **Total SS IDs**: 4,849,775,973  
- **Unmapped SS IDs**: 153,387

---

## 6. dbSNP FTP Site Information

- FTP URL: [`ftp.ncbi.nih.gov/snp/archive/b157/VCF/`](https://ftp.ncbi.nih.gov/snp/archive/b157/VCF/)



- Data Source: NCBI dbSNP Build 157  
- Documentation adapted from official release notes and `README.txt`
***************************************
       CURRENT ANNOUNCEMENTS
***************************************
Revised: Sept 17, 2018
****************************************

This document describes the FTP repository of dbSNP data in the
sections listed below:

     I.   FTP SITE DESCRIPTION AND FTP UPDATE FREQUENCY
    II.   DIRECTORY STRUCTURE
    IV.   FTP SITE REVISION HISTORY

********************************************************************************

I. FTP SITE DESCRIPTION AND FTP UPDATE FREQUENCY

Access to the NCBI FTP site is available via the web or anonymous FTP. 
The URL/host addresses are:

World Wide Web:     https://ftp.ncbi.nlm.nih.gov/snp/
Anonymous FTP:      host ftp.ncbi.nlm.nih.gov
                    cd snp

Announcements of the release of new builds and notification of corrections
to existing data content will be posted on a public mail list. To receive
these notifications, you can subscribe to the dbSNP announcement list at
https://www.ncbi.nlm.nih.gov/mailman/listinfo/dbsnp-announce.

********************************************************************************

II. DIRECTORY STRUCTURE
DIRECTORIES: 

/bin             Contains demo software tools for using RefSNP JSON files.

/specs           Contains the RefSnp API schema used by RefSNP JSON files.
                         /refsnp_specification.yaml

/latest_release  Contains the most recent release of human SNP data, in VCF and
                 API JSON format, along with the release notes:
                         /JSON
                         /VCF
                         /release_notes.txt


SUBDIRECTORIES:

/JSON              RefSNP JSON files. Refer to JSON_README.txt for details.

/VCF               RefSNP VCF files for GRC (Genome Reference Consortium) human assembly
                   37 (GCF_000001405.25) and 38 (GCF_000001405.40). Files are compressed
                   by bgzip and with the tabix index.

IV. FTP SITE REVISION HISTORY:
Rev: Sept 17, 2018
Rewrite the readme file for redesigned dbSNP.
********************************************************************************