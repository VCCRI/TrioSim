#!/bin/bash

## Ubuntu based updates

sudo apt-get update

sudo apt install gcc

sudo apt-get install zlib1g-dev

sudo apt-get install libbz2-dev

sudo apt-get install liblzma-dev

sudo apt install make

###### JRE, bgzip  ############################

sudo apt install openjdk-11-jre-headless

sudo apt install tabix

######## bcftools (v 1.7), vcftools (v 0.1.15), Python2, numpy, pandas, pip and pyVCF ############

sudo apt install bcftools
sudo apt install vcftools

sudo add-apt-repository universe
sudo apt install python2
sudo apt install python-numpy
sudo apt install python-pandas

sudo apt install curl
curl https://bootstrap.pypa.io/pip/2.7/get-pip.py --output get-pip.py
sudo python2 get-pip.py
pip2 install PyVCF


########### SIMdrom ###############

mkdir jars
cd jars

wget https://github.com/visze/simdrom/releases/download/v0.0.2/simdrom-cli-0.0.2.jar

cd ..


########### Download files ###############

sudo apt  install awscli

aws s3 cp --recursive s3://vccri-giannoulatou-lab-denovo-mutations/programs/TrioSim/data_files data_files

#aws s3 cp s3://vccri-giannoulatou-lab-denovo-mutations/programs/TrioSim/data_files/chromosome_name.txt .

#aws s3 cp s3://vccri-giannoulatou-lab-denovo-mutations/programs/TrioSim/data_files/denovo-db.non-ssc-samples.variants_hg38LiftedOver_decoy_hla_UCSC_SNP_INDEL.vcf .

#aws s3 cp s3://vccri-giannoulatou-lab-denovo-mutations/programs/TrioSim/data_files/chr22_hg38_sites.vcf .

#aws s3 cp s3://vccri-giannoulatou-lab-denovo-mutations/programs/TrioSim/data_files/ALL.WGS_GRCh38_sites.20170504.vcf.gz .
#aws s3 cp s3://vccri-giannoulatou-lab-denovo-mutations/programs/TrioSim/data_files/ALL.WGS_GRCh38_sites.20170504.vcf.gz.tbi .


##### Executables ####

chmod +x *.sh *.py





