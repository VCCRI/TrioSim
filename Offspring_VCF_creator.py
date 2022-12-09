#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import vcf
import random
import json

# Usage : python generate_child_VCF.py > /Users/anusha/Documents/anushi/denovo_project/data/1000G_files/release_20130502/simulated_trios_genomes_v1/child_VCF.txt
# This script is run on Anushi's MAC machine through anaconda environment.

# TBD : Check for duplicate SNPs

# export PERL5LIB=/Users/anusha/Documents/anushi/software/vcftools/src/perl
# export PATH=$PATH:/Users/anusha/Documents/anushi/software/vcftools_anushi_install/bin
mother_father_merged_vcf = sys.argv[1]
directory = sys.argv[2]
child_VCF = directory+"/"+"Child_VCF.vcf"

childvcf_reader = vcf.Reader(open(mother_father_merged_vcf, 'r'))

#childvcf_first_record = next(childvcf_reader)

#print type(childvcf_reader)
#print childvcf_first_record
#print childvcf_first_record.samples

#print childvcf_reader.samples
vcf_writer = open(child_VCF, "wb")

vcf_writer.write("##fileformat=VCFv4.2\n")
vcf_writer.write("##contig=<ID=chr1,assembly=b37,length=249250621>\n")
vcf_writer.write("##contig=<ID=chr2,assembly=b37,length=243199373>\n")
vcf_writer.write("##contig=<ID=chr3,assembly=b37,length=198022430>\n")
vcf_writer.write("##contig=<ID=chr4,assembly=b37,length=191154276>\n")
vcf_writer.write("##contig=<ID=chr5,assembly=b37,length=180915260>\n")
vcf_writer.write("##contig=<ID=chr6,assembly=b37,length=171115067>\n")
vcf_writer.write("##contig=<ID=chr7,assembly=b37,length=159138663>\n")
vcf_writer.write("##contig=<ID=chr8,assembly=b37,length=146364022>\n")
vcf_writer.write("##contig=<ID=chr9,assembly=b37,length=141213431>\n")
vcf_writer.write("##contig=<ID=chr10,assembly=b37,length=135534747>\n")
vcf_writer.write("##contig=<ID=chr11,assembly=b37,length=135006516>\n")
vcf_writer.write("##contig=<ID=chr12,assembly=b37,length=133851895>\n")
vcf_writer.write("##contig=<ID=chr13,assembly=b37,length=115169878>\n")
vcf_writer.write("##contig=<ID=chr14,assembly=b37,length=107349540>\n")
vcf_writer.write("##contig=<ID=chr15,assembly=b37,length=102531392>\n")
vcf_writer.write("##contig=<ID=chr16,assembly=b37,length=90354753>\n")
vcf_writer.write("##contig=<ID=chr17,assembly=b37,length=81195210>\n")
vcf_writer.write("##contig=<ID=chr18,assembly=b37,length=78077248>\n")
vcf_writer.write("##contig=<ID=chr19,assembly=b37,length=59128983>\n")
vcf_writer.write("##contig=<ID=chr20,assembly=b37,length=63025520>\n")
vcf_writer.write("##contig=<ID=chr21,assembly=b37,length=48129895>\n")
vcf_writer.write("##contig=<ID=chr22,assembly=b37,length=51304566>\n")
vcf_writer.write("##contig=<ID=MT,assembly=b37,length=16569>\n")
vcf_writer.write("##contig=<ID=chrX,assembly=b37,length=155270560>\n")
vcf_writer.write("##contig=<ID=chrY,assembly=b37,length=59373566>\n")


vcf_writer.write("#CHROM"+"\t"+"POS"+"\t"+"ID"+"\t"+"REF"+"\t"+"ALT"+"\t"+"QUAL"+"\t"+"FILTER"+"\t"+"INFO"+"\t"+"FORMAT"+"\t"+"Child_genotype"+"\n")

for each_record in childvcf_reader:
	
	# print "-----------------------------------------------------------------"
	# print each_record
	# print each_record.QUAL
	# print each_record.FILTER
	# print type(each_record.INFO)
	# print each_record.INFO
	# print each_record.genotype

	mother_genotype = each_record.genotype(childvcf_reader.samples[0])['GT']
	father_genotype = each_record.genotype(childvcf_reader.samples[1])['GT']
	
	#print mother_genotype+"\t"+father_genotype
	if(mother_genotype == './.'):
		mother_genotype = '0/0'

	if(father_genotype == './.'):
		father_genotype = '0/0'
 	
 	mother_alleles = mother_genotype.split("/")
 	father_alleles = father_genotype.split("/")

 	#print str(mother_alleles)+"\t"+str(father_alleles)

 	
 	# Generate punnett square
 	#print str(mother_alleles)+"\t"+str(father_alleles)

 	child_genotype_list = []

 	for m in mother_alleles:
 		for f in father_alleles:                                                       
 				
 			if m == '1' and f == '0':
 				child_genotype_option = str(f)+"/"+str(m)
 				child_genotype_list.append(child_genotype_option)
 			else:
 				child_genotype_option = str(m)+"/"+str(f)
 				child_genotype_list.append(child_genotype_option)


 	#print child_genotype_list

 	#child_genotype = random.Random(500).choice(child_genotype_list)
 	child_genotype = random.choice(child_genotype_list)

 	#print "child_genotype : "+child_genotype
 	ALT =  ''.join([unicode(i) for i in each_record.ALT])
 	#FILTER = ''.join([unicode(i) for i in each_record.FILTER])
 	FILTER = "PASS"
 	
 	
 	## Extract INFO column
 	INFO_dict_list=[]
	
	for i,j in each_record.INFO.iteritems():
		key = i
		
		if(isinstance(j,list)):
			value = ','.join([unicode(k) for k in j])
		else:
			value = j

		item = key + "=" + str(value)
		
		INFO_dict_list.append (item)
	
	INFO = ';'.join([unicode(i) for i in INFO_dict_list])

	## ID field
	if each_record.ID is None:
		each_record.ID = "NA"

	
	# print each_record.CHROM
	# print str(each_record.POS)
	# print each_record.ID
	# print each_record.REF
	# print ALT
	# print str(each_record.QUAL)
	# print FILTER
	# print INFO
	# print str(child_genotype)

	if child_genotype <> '0/0':
		print_str = each_record.CHROM+"\t"+str(each_record.POS)+"\t"+each_record.ID+"\t"+each_record.REF+"\t"+ALT+"\t"+str(each_record.QUAL)+"\t"+FILTER+"\t"+INFO+"\t"+"GT"+"\t"+str(child_genotype)+"\n"
		vcf_writer.write(print_str)

vcf_writer.close() 		

 