#!/bin/bash


## Usage help function ##
function help()
{
    echo "Usage: ./TrioSim -s <trio_number>
	-d <output directory path>
	-n <Number of de novo mutations to be spiked-in>
	[-t <test_run>]
	[-h <help>]
	
	TrioSim is used to generated trio VCF files using variants from 1000G for parents.
	For offspring VCF, Mendelian ineritance laws is followed (using generated parental VCF file).  
	De novo mutations in offsppring VCF are spiked-in from denovo-db.
	
	Required arguments:
	-s <trio_number>                                  The sample number of the trio VCF generation.
	-d <output directory path>                        Path to the output directory.
	-n <Number of de novo mutations to be spiked-in>  Number of de novo mutations to be spiked-in into offspring VCF file from denovo-db.
	
	Options:
	-t <test_run>                                     The test run executes TrioSim and generates trio VCF based on chr22 variants from 1000G.
	-h <help>                                         This help message for TrioSim.
	
	
    "
	
	exit
}

## Test data run function ##
## Here, small test data based on chromosome 22 is run. Trio VCF files and DNMs are for chr22.
function run_test_data()
{
  echo "Run test data. It generates Trio VCF files for chromosome 22 including spiked-in de novo mutations!"
  echo "Output will be generated in test_run_output sub-directory in TrioSim installation direcctory"
  
  vcf_sites_file=data_files/chr22_hg38_sites.vcf
  denovo_db_file=data_files/denovo-db.non-ssc-samples.chr22.UCSC_SNP_INDEL.vcf
  trio_number=1
  target_directory=test_run_output
  num_DNMs=10
  
  #echo $test
  #echo $vcf_sites_file
  #echo $denovo_db_file
  
  trio_vcf_generation
  
  END=$(date "+%s")

  echo "It takes $((END-START)) seconds to complete TrioSim..."

  
}

function trio_vcf_generation()
{
	echo "Inside trio_vcf_generation function"
	#echo $test
	#echo $vcf_sites_file
	#echo $denovo_db_file
	
	mkdir -p $target_directory

	sample_name_P1='trio_'$trio_number'_parent1'
	sample_name_P2='trio_'$trio_number'_parent2'
	sample_C1='trio_'$trio_number'_child'

	directory_P1=$target_directory/'parent1'
	directory_P2=$target_directory/'parent2'
	directory_child=$target_directory/'child'
	
	
	## START VCF Generation for PARENTS ##

	echo "Start VCF Generation Parent 1.."

	./Parental_VCF_generator.sh $sample_name_P1 $directory_P1 $vcf_sites_file

	echo "Start VCF Generation Parent 2.."

	./Parental_VCF_generator.sh $sample_name_P2 $directory_P2 $vcf_sites_file

	## START VCF Generation for Child ##

	vcfgz_file_P1=$directory_P1/$sample_name_P1'_normalized'.vcf.gz
	vcfgz_file_P2=$directory_P2/$sample_name_P2'_normalized'.vcf.gz


	./Offspring_VCF_generator.sh $directory_child $vcfgz_file_P1 $vcfgz_file_P2 $num_DNMs $denovo_db_file

	
	#cd $target_directory
	#rm -rf *

	END=$(date "+%s")

	echo "It takes $((END-START)) seconds to complete TrioSim..."

}


################################################# User parameters options ##########################################

START=$(date "+%s")

test=false
vcf_sites_file=""
denovo_db_file=""

while getopts :s:d:n:ht option
do 
    case "${option}"
        in
		t)
			test=true
			run_test_data
			exit;;
		s)trio_number=${OPTARG};;
        d)target_directory=${OPTARG};;
		n)num_DNMs=${OPTARG};;
		h) 
			help ;;
		:) 
			echo "Error: -${OPTARG} requires an argument."
			help ;;
		*)
			echo "Unexpected option: ${OPTARG}"
			help ;;
    esac
done



################################################# Call trio vcf generation function from user arguments #####################
 
if [ "$trio_number" ] && [ "$target_directory" ] && [ "$num_DNMs" ]
then
	echo "Trio sample number : $trio_number"
	echo "Target Directory   : $target_directory"
	echo "Number of DNMs   : $num_DNMs"
	
	vcf_sites_file=data_files/ALL.WGS_GRCh38_sites.20170504.vcf.gz
	denovo_db_file=data_files/denovo-db.non-ssc-samples.variants_hg38LiftedOver_decoy_hla_UCSC_SNP_INDEL.vcf
	
	trio_vcf_generation
fi



exit
