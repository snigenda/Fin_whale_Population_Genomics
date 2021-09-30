# -*- coding: utf-8 -*-
'''
# Title:  extract the most deleterious snpEff Annotation_Impact depending on the 'ANN=' section
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Jan 29 14:43:14 2021
'''

###########################################################
## import packages
import os
import sys
import argparse
import csv
import gzip

###########################################################
## def functions
sys.path.append('/u/project/rwayne/meixilin/fin_whale/analyses/scripts/config')
from vcf_parser_finwhale import split_info

def parse_args():
	'''
	Parse command-line arguments
	'''
	parser = argparse.ArgumentParser(description='This script will traverse a vcf file and output snpEff Annotation_Impact')

	parser.add_argument(
			'--VCF', required=True,
			help='REQUIRED. Path to the VCF file. Should be gzipped.')

	parser.add_argument(
			'--filter', required=True,
			help='REQUIRED. Site filter to apply. Separate by comma. Supply entire strings.')

	parser.add_argument(
			'--outprefix', required=True,
			help='REQUIRED. Path and prefix to the output file.')

	args = parser.parse_args()
	return args

# input a list of ANN
# return the most deleterious effects for a list of annotations
def get_mostdel_impact(ANN, impactrank):
	# the third field is the annotation impact
	ANN_impacts = [ann.split('|')[2] for ann in ANN]
	if not set(ANN_impacts).issubset(set(impactrank)):
		raise ValueError('Wrong Annotation_Impact format!')
	ANN_order= [impactrank.index(ii) for ii in ANN_impacts]
	ANN_orderminid = [ii for ii, xx in enumerate(ANN_order) if xx == min(ANN_order)]
	# the snpEff should have already ranked the rankings by the sort order defined
	if 0 not in ANN_orderminid:
		raise ValueError('Wrong Annotation_Impact order: ' + ','.join(ANN))
	else:
		return ANN_impacts[0]

###########################################################
## def variables
header = '# CHROM\tSTART\tEND\tAnnotation_Impact\n'
# annotation fields extracted from vcf header: grep 'INFO=<ID=ANN' JointCalls_testfile.vcf
annfields = ['Allele','Annotation','Annotation_Impact','Gene_Name','Gene_ID','Feature_Type','Feature_ID','Transcript_BioType','Rank','HGVS.c','HGVS.p','cDNA.pos / cDNA.length','CDS.pos / CDS.length','AA.pos / AA.length','Distance','ERRORS / WARNINGS / INFO']
# impact ranks
impactrank=['HIGH', 'MODERATE', 'LOW', 'MODIFIER']

###########################################################
## main
def main():
	args = parse_args()
	outfile = open(args.outprefix + '.bed', 'w')
	# outfile.write(header)
	myfilter = str(args.filter).split(',')

	with gzip.open(args.VCF, 'rt') as VCF:
		for line0 in VCF:
			if not line0.startswith('#'):
				line = line0.strip().split('\t')
				# get lines you are interested in
				mychr=line[0]
				myendpos=line[1]
				mystartpos=str(int(myendpos) - 1)
				myinfo = line[7]
				# check if pass filter, i.e if the elements in myfilterinfo is a complete subset of myfilter
				# don't split by ';' check the entire string
				myfilterinfo = [line[6]]
				if not set(myfilterinfo).issubset(set(myfilter)):
					continue
				# split the info if passing filter
				infoFields=split_info(myinfo)
				# only work with lines that contain ANN (ANN annotations)
				if 'ANN' in infoFields.keys():
					# split by the ','
					myANN=infoFields['ANN'].split(',')
					myimpact=get_mostdel_impact(myANN, impactrank)
					output = [str(mychr),str(mystartpos),str(myendpos),str(myimpact)]
					outfile.write("\t".join(output) + "\n")
	outfile.close()


if __name__ == '__main__':
	sys.exit(main())

