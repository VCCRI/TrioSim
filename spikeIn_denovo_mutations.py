#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import os
import io
import numpy as np
import pandas as pd


## Read de novo mutations VCF file
## Create a dictionary with Chr, position as key and REF and ALT as values

def read_vcf(path):
  
    f = open(path,'r')

    data = []
    header_lines = []

    for each_line in f:
    	if each_line.startswith('##'):
    		header_lines.append(each_line)
    	else:
    		data.append(each_line)

    data_lines = pd.read_table(io.StringIO(u''.join(data)),dtype={'#CHROM': str, 'POS': str, 'ID': str, 'REF': str, 'ALT': str,'QUAL': str, 'FILTER': str, 'INFO': str}).rename(columns={'#CHROM': 'CHROM'})

    f.close()

    return header_lines, data_lines


## Change the file according to hg19 or hg38
DNM_header, DNM_data = read_vcf('/home/ubuntu/softwares/data_files/denovo-db.non-ssc-samples.variants_hg38LiftedOver_decoy_hla_UCSC_SNP_INDEL.vcf')

childVcf_header, childVcf_data = read_vcf(sys.argv[1])
output_dir=sys.argv[2]

#print DNM_header
#print DNM_data.head()

#print childVcf_header
#print childVcf_data

## Find common SNPs between denovo-db and Child VCFs
#common_df = pd.merge(DNM_data, childVcf_data, on=['CHROM','POS'], how='inner')
common_df = pd.merge(DNM_data, childVcf_data, on=['CHROM','POS'])

#print "## Before filter ##"
#print DNM_data.shape
#print childVcf_data.shape
#print childVcf_data.head()
#print common_df.shape
#print common_df['POS'].head()

#common_df.to_csv('/Genomes/user_vcf/common_denonvo_dbSNPs_childSNPs.csv')

## Filter out SNPs from denovo-db which are common in Child VCF

#DNM_data = DNM_data[(~DNM_data.CHROM.isin(common_df.CHROM))&(~DNM_data.POS.isin(common_df.POS))]
# https://stackoverflow.com/questions/28901683/pandas-get-rows-which-are-not-in-other-dataframe?noredirect=1&lq=1
DNM_filter = DNM_data.merge(common_df.drop_duplicates(), on=['CHROM','POS'], how='left', indicator=True)

DNM_filter_leftonly = DNM_filter.loc[DNM_filter['_merge'] == 'left_only'][['CHROM','POS','ID','REF','ALT']]

#print DNM_filter_leftonly.head()
#print DNM_filter_leftonly.shape

## Subset partucular chromosome SNPs to be spiked-in (e.g. chr22)
#DNM_filter_leftonly = DNM_filter_leftonly.loc[DNM_filter_leftonly['CHROM'] == 'chr22']

## Exclude spiking-in X, Y chromosome & GL000209.1 SNPs
#DNM_filter_leftonly = DNM_filter_leftonly[(DNM_filter_leftonly.CHROM != 'chrX') & (DNM_filter_leftonly.CHROM != 'chrY') & (DNM_filter_leftonly.CHROM != 'GL000209.1')]
DNM_filter_leftonly = DNM_filter_leftonly[~DNM_filter_leftonly.CHROM.str.contains('KI|GL|chrX|chrY',case=False)]

## Randomly select "n" number of de novo mutations to be spiked-in
k = 100


## Randomly select "n" de novo mutations
DNM_sampled = DNM_filter_leftonly.sample(n=k)

## Insert column with default values
DNM_sampled.insert(5,"QUAL", "100") 
DNM_sampled.insert(6,"FILTER", "PASS") 
DNM_sampled.insert(7,"INFO", ".")
DNM_sampled.insert(8,"FORMAT", "GT")
DNM_sampled.insert(9,"Child_genotype", "0/1")

#print DNM_sampled.shape
#print DNM_sampled
DNM_sampled.to_csv(output_dir+'/DNM_sampled.csv')

## Spike-in de novo mutations in Child VCF (outer merge)
spikedIn_df = pd.merge(childVcf_data, DNM_sampled, on=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','Child_genotype'], how='outer')
spikedIn_df.rename(columns={'CHROM':'#CHROM'}, inplace=True)
spikedIn_df['ID'].fillna('.', inplace=True)
spikedIn_df.dropna()
#spikedIn_df.sort_values(by=['#CHROM'],inplace=True)


#spikedIn_df.to_csv('/Genomes/user_vcf/spikedIn_df.csv')

for x in childVcf_header:
    print x.rstrip('\n')

#print spikedIn_df.head()

print(spikedIn_df.to_csv(sep='\t', index=False))

del DNM_header
del DNM_data

del childVcf_header
del childVcf_data

del DNM_filter
del DNM_filter_leftonly
del spikedIn_df



