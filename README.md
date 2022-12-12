# TrioSim

TrioSim is linux shell script based pipeline to generate WGS trio data. It can be run in two steps :

** 1. Generate trio VCF files **

* In this step, VCF files for the trio are generated based on Parental VCF Generator and Offspring VCF generator of TrioSim.
* It has features to spike-in de novo mutations(from denovo-db) in the VCF file of the offspring.
* VCF file from denovo-db is provided in data_files directory and will be automatically used by TrioSim script.

** 2. Generate trio BAM files **

* In this step, independent simulated BAM files for the trio are generated. We have provided scripts for EAGLE and NEAT simulators for the same. 
* But the users are free to use VCF files generated from Step1 and use those as input to any other simulator of their choice.

** System requirements : **

* Python (version 2)
* Modules for python 2 : numpy, pandas, PyVCF
* bcftools (version 1.7)
* htslib (version 1.16)
* Open JRE (version 11)


** Step 1 : Generate trio VCF files **

   Usage :

    TrioSim.sh -i <sample number/name> \
    -d </output directory path> \
    -n </number of de novo mutations> \
    -b <hg19/hg38 for variant sites>

** Step 2 : Generate trio BAM files **

A. Usage for EAGLE genome simulator
    
    For further information on usage of EAGLE simualator, please refer to : https://github.com/sequencing/EAGLE
    
B. Usage for NEAT genome simulator
    
    For further information on usage of NEAT simualator, please refer to : https://github.com/zstephens/neat-genreads
    
    
