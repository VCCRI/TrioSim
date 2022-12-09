#!/bin/bash

START=$(date "+%s")

trio_number=$1
target_directory=$2
num_DNMs=$3

#trio_number=2
#target_directory=/Genomes/family1_WGS_hg38_liftOver_V2/

mkdir -p $target_directory

sample_name_P1='trio_'$trio_number'_parent1'
sample_name_P2='trio_'$trio_number'_parent2'
sample_C1='trio_'$trio_number'_child'

directory_P1=$target_directory/'parent1'
directory_P2=$target_directory/'parent2'
directory_child=$target_directory/'child'



################################################# START VCF Generation for PARENTS ##############################

echo "Start VCF Generation Parent 1.."

./Parental_VCF_generator.sh $sample_name_P1 $directory_P1 

echo "Start VCF Generation Parent 2.."

./Parental_VCF_generator.sh $sample_name_P2 $directory_P2 

################################################# START VCF Generation for Child ##############################

vcfgz_file_P1=$directory_P1/$sample_name_P1'_normalized'.vcf.gz
vcfgz_file_P2=$directory_P2/$sample_name_P2'_normalized'.vcf.gz


##./Offspring_VCF_generator.sh $directory_child $vcfgz_file_P1 $vcfgz_file_P2 $num_DNMs

END=$(date "+%s")

echo $START
echo $END
echo "It takes $((END-START)) seconds to complete this script..."

cd $target_directory
#rm -rf *


exit








