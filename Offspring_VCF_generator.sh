#!/bin/bash


directory=$1
father_vcf=$2
mother_vcf=$3
num_DNMs=$4
denovo_db_file=$5

# Merge mother and father VCF

mkdir -p $directory

echo "Merge mother and father VCF files.."
bcftools merge -m none --force-samples $father_vcf $mother_vcf -o $directory/merged_mother_father.vcf

## Remove ChrX and ChrY chromosomes SNPs if exists
sed -i  -e '/chrX/d' -e '/chrY/d' $directory/merged_mother_father.vcf

echo "Normalize merged VCF file.."
bcftools norm -d all $directory/merged_mother_father.vcf -o $directory/merged_mother_father_normalized.vcf
rm $directory/merged_mother_father.vcf


# Python program : generate Child VCF
echo "Generate Child VCF.."
python Offspring_VCF_creator.py $directory/merged_mother_father_normalized.vcf $directory
rm $directory/merged_mother_father_normalized.vcf


# Python program : Spike-in DNMs in Child VCF
echo "Spike-in DNMs.."
python spikeIn_denovo_mutations.py $directory/Child_VCF.vcf $directory $num_DNMs $denovo_db_file > $directory/Child_VCF_DNM.vcf 

echo "Sort Child VCF by chromosome.."
cat $directory/Child_VCF_DNM.vcf | vcf-sort -c > $directory/Child_VCF_DNM_sort.vcf 


rm $directory/Child_VCF.vcf
rm $directory/Child_VCF_DNM.vcf


exit
