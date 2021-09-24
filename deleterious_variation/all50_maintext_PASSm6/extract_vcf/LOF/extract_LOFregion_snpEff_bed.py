# -*- coding: utf-8 -*-
'''
# Title:  extract snpEff LOF regions depending on the 'LOF=' section
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Dec  4 10:27:46 2020
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
	parser = argparse.ArgumentParser(description='This script will extract bedfiles for given vcf files')

	parser.add_argument(
			'--VCF', required=True,
			help='REQUIRED. Path to the VCF file. Should be gzipped.')

	parser.add_argument(
			'--filter', required=True,
			help='Site filter to apply. Separate by comma')

	parser.add_argument(
			'--cutoff', required=True,
			help='Percentage affected cutoff (not including)')

	parser.add_argument(
			'--outprefix', required=True,
			help='REQUIRED. Path and prefix to the output file.')

	args = parser.parse_args()
	return args

###########################################################
## def variables
header = '# CHROM\tSTART\tEND\n'

###########################################################
## main

def main():
	args = parse_args()
	outfile = open(args.outprefix + '.bed', 'w')
	# outfile.write(header)
	myfilter = str(args.filter).split(',')
	# percent cutoff
	per_cutoff = float(args.cutoff)
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
				# only accept filters including: PASS,WARN_missing,WARN_excessHet
				myfilterinfo = [line[6]]
				if not set(myfilterinfo).issubset(set(myfilter)):
					continue
				# split the info if passing filter
				infoFields=split_info(myinfo)
				# only work with lines that contain LOF (LOF annotations)
				if 'LOF' in infoFields.keys():
					# split by the ','
					myLOF=infoFields['LOF'].strip().split(',')
					per_affected=[]
					for thislof in myLOF:
						thislof=thislof.strip('(').strip(')').split('|')
						# get the last line of the percent affected
						per_affected.append(float(thislof[3]))
					# as long as there is one that is bigger (works for the list)
					if per_affected > per_cutoff:
						output = [mychr, mystartpos, myendpos, infoFields['LOF']]
						outfile.write('\t'.join(output) + '\n')
					else:
						continue
				else:
					continue
			else:
				continue
	outfile.close()


if __name__ == '__main__':
	sys.exit(main())

