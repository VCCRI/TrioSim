# TrioSim

TrioSim is linux shell script based pipeline to generate WGS trio data. It can be run in two steps :

**1. Generate trio VCF files**

* In this step, VCF files for the trio are generated based on Parental VCF Generator and Offspring VCF generator of TrioSim.
* It has features to spike-in de novo mutations(from denovo-db) in the VCF file of the offspring.
* VCF file from denovo-db is provided in data_files directory (from install_dependencies.sh) and will be automatically used by TrioSim script.

**2. Generate trio BAM files**

* In this step, independent simulated BAM files for the trio are generated. We have provided scripts for EAGLE and NEAT simulators for the same. 
* But the users are free to use VCF files generated from Step1 and use those as input to any other simulator of their choice.

**System requirements :**
These packages can also be installed as a part of install_dependencies.sh script.

* Python (version 2)
* Modules for python 2 : numpy, pandas, PyVCF
* bcftools (version 1.7)
* vcftools (version 0.1.15)
* htslib (version 1.16)
* Open JRE (version 11)

**Data files requirements :**

These files are required to be downloaded in following sub-directories of TrioSim :
Alternatively, they are installed as a part of install_dependencies.sh script.

a) jars
* https://github.com/visze/simdrom/releases/download/v0.0.2/simdrom-cli-0.0.2.jar

b) data_files : 
* s3://vccri-giannoulatou-lab-denovo-mutations/programs/TrioSim/data_files/chromosome_name.txt
* s3://vccri-giannoulatou-lab-denovo-mutations/programs/TrioSim/data_files/denovo-db.non-ssc-samples.variants_hg38LiftedOver_decoy_hla_UCSC_SNP_INDEL.vcf
* s3://vccri-giannoulatou-lab-denovo-mutations/programs/TrioSim/data_files/chr22_hg38_sites.vcf
* s3://vccri-giannoulatou-lab-denovo-mutations/programs/TrioSim/data_files/ALL.WGS_GRCh38_sites.20170504.vcf.gz
* s3://vccri-giannoulatou-lab-denovo-mutations/programs/TrioSim/data_files/ALL.WGS_GRCh38_sites.20170504.vcf.gz.tbi


**Step 1 : Generate trio VCF files**

   Usage :

    Usage: ./TrioSim -s <trio_number>
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

**Step 2 : Generate trio BAM files**

A. Usage for EAGLE genome simulator
    
    For further information on usage of EAGLE simualator, please refer to : https://github.com/sequencing/EAGLE
    
B. Usage for NEAT genome simulator
    
    For further information on usage of NEAT simualator, please refer to : https://github.com/zstephens/neat-genreads
    
    
