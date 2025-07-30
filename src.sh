bcftools view \
  -v snps \
  -i 'ID ~ "^rs"' \
  GCF_000001405.40.gz -Oz -o dbsnp157_rs_hg38.vcf.gz

  bcftools annotate --rename-chrs RefSeq.seq.accession2UCSC.style.name.txt dbsnp157_rs_hg38.vcf.gz \
   -Oz -o dbsnp157_rs_hg38_UCSC.vcf.gz

 bcftools query -f '%CHROM\t%POS\n' dbsnp157_rs_hg38_UCSC.vcf.gz > dbsnp157_rs_hg38_UCSC.tab 


bcftools view -v snps -i 'ID ~ "^rs"' GCF_000001405.40.gz -Ou | \
bcftools annotate --rename-chrs RefSeq.seq.accession2UCSC.style.name.txt -Ou | \
bcftools query -f '%CHROM\t%POS\n' > dbsnp157_rs_hg38_UCSC.tab