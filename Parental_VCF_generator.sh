#!/bin/bash

sample_name=$1
directory=$2
#vcf_sites_file=$3
vcf_sites_file=data_files/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz


################################################# PARENT 1 VCF Generation ####################################################
mkdir -p $directory

echo "Start SIMDrom"
#java -jar /softwares/simdrom-cli-0.0.2.jar -b /Genomes/data_files/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.sites.vcf.gz -bAF $ethnicity > $directory/genome$sample_number'_'$ethnicity.vcf
java -jar jars/simdrom-cli-0.0.2.jar -b $vcf_sites_file > $directory/$sample_name.vcf


echo "Convert to UCSC"
# # # Convert to UCSC naming
bcftools annotate --rename-chrs data_files/chromosome_name.txt $directory/$sample_name.vcf > $directory/$sample_name'_UCSC'.vcf

## REMOVE Chr X and Chr Y data (as we use only autosomes)
sed -i  -e '/chrX/d' -e '/chrY/d' $directory/$sample_name'_UCSC'.vcf

echo "Remove multi-allelic sites"
# # # Remove multi-allelic 
bcftools view -M2 $directory/$sample_name'_UCSC'.vcf -o $directory/$sample_name'_no_multiallelic_UCSC'.vcf
rm $directory/$sample_name'_UCSC'.vcf
rm $directory/$sample_name.vcf


echo "SNPs and INDELs"
# # # Select only SNPs and INDELs
bcftools view -v snps,indels $directory/$sample_name'_no_multiallelic_UCSC'.vcf -o $directory/$sample_name'_SNP_INDEL'.vcf
rm $directory/$sample_name'_no_multiallelic_UCSC'.vcf


echo "Normalize VCF"
# # # Normalize VCF
bcftools norm -d all $directory/$sample_name'_SNP_INDEL'.vcf -o $directory/$sample_name'_normalized'.vcf
rm $directory/$sample_name'_SNP_INDEL'.vcf

echo "Compress VCF file"
# # # Compress VCF file
bgzip -c $directory/$sample_name'_normalized'.vcf > $directory/$sample_name'_normalized'.vcf.gz
rm $directory/$sample_name'_normalized'.vcf

echo "Index compressed VCF file"
# # # Index compressed VCF file
tabix $directory/$sample_name'_normalized'.vcf.gz


exit