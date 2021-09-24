# -*- coding: utf-8 -*-
'''
# Title: Sliding window heterozygosity
# Author: jarobinson, snigenda
# Date: Sat Dec  5 23:22:27 2020
# Example usage:

Script to count number of called genotypes and number of heterozygotes per sample in sliding windows.
Input file is a single- or multi-sample VCF file that has been filtered (passing sites have "PASS" in the FILTER column) and compressed with gzip/bgzip.

Usage:
python ./SlidingWindowHet.py [vcf] [window size] [step size] [chromosome/scaffold name]

Windows will be non-overlapping if step size == window size.

Example:
python ./SlidingWindowHet.py input.vcf.gz 100000 10000 chr01
'''

###########################################################
## import packages
import sys
import pysam
import os
import os.path
import gzip
import csv
import pandas as pd
import argparse

###########################################################
## def functions
sys.path.append('/u/project/rwayne/meixilin/fin_whale/analyses/scripts/config')
from vcf_parser_finwhale import get_samples, read_minke_contig_len, read_minke_contig_idx

def parse_args():
	'''
	Parse command-line arguments
	'''
	parser = argparse.ArgumentParser(description='This script will generate heterozygosity in sliding windows')

	parser.add_argument(
			'--path', required=True,
			help='REQUIRED. Working directory')

	parser.add_argument(
			'--VCF', required=True,
			help='REQUIRED. Path to the VCF file. Should be gzipped.')

	parser.add_argument(
			'--outprefix', required=True,
			help='REQUIRED. Path to the output file and prefix.')

	parser.add_argument(
			'--windowsize', required=True,
			help='REQUIRED. Window size.')

	parser.add_argument(
			'--stepsize', required=True,
			help='REQUIRED. Step size.')

	parser.add_argument(
			'--idx', required=True,
			help='REQUIRED. Chromosome index used in the vf file formats. Should format like 01-96')

	args = parser.parse_args()
	return args

# Make sure the VCF file is indexed (if not, create index)
def index_invcf(filename):
	if not os.path.exists("%s.tbi" % filename):
		pysam.tabix_index(filename, preset="vcf")
		sys.stdout.write(str('Generated tbi for '+filename + '\n'))
	parsevcf = pysam.Tabixfile(filename)
	return parsevcf

# Get the start and end positions in the given contigs (allows for multiple contigs per vcf file), also check that the start position is 1.
def get_start_end(filename, chrom, chrom_size):
	with gzip.open(filename, 'rt') as VCF:
		for line0 in VCF:
			if line0[0] != '#':
				# split the vcf entry lines
				line=line0.strip().split()
				# get to the target chromosome (useful if you have multiple chromosomes)
				if line[0] == chrom:
					# if there are multiple chromosomes in one
					start_pos = int(line[1])
					# start_pos should always be 1 in the finwhale all50 dataset
					if (start_pos != 1):
						errmsg = str('\t'.join(line[:2]) + 'Wrong start pos')
						sys.exit(errmsg)
					end_pos = int(chrom_size[chrom])
				else:
					continue
				break
			else:
				continue
	return start_pos, end_pos

# Calculate heterozygosity
# Fetch a region, ignore sites that fail filters, tally genotype calls and heterozygotes
def snp_cal(parsevcf, samples, output, chrom, window_start, window_end, myfilter):
	print("%s:%s-%s" % (chrom, window_start, window_end))
	rows = tuple(parsevcf.fetch(region="%s:%s-%s" % (chrom, window_start, window_end), parser=pysam.asTuple()))
	sites_total=0
	calls=[0]*len(samples)
	hets=[0]*len(samples)
	for line in rows:
		if line[6]!=myfilter:
			continue
		sites_total+=1
		for i in range(0,len(samples)):
			if line[i+9][:1]=='.':
				continue
			calls[i]+=1
			GT=line[i+9].split(':')[0]
			if '/' in GT: sp='/'
			if '|' in GT: sp='|'
			if GT.split(sp)[0]!=GT.split(sp)[1]: hets[i]+=1
	output.write( '%s\t%s\t%s\t%s\t%s\t%s\n' % (chrom,window_start,window_end,sites_total,'\t'.join(map(str,calls)),'\t'.join(map(str,hets))) )
	return None

###########################################################
## def variables

###########################################################
## main
def main():
	# Set variables
	args = parse_args()
	os.chdir(args.path)
	filename = str(args.VCF) # should be a full path
	outprefix=str(args.outprefix)
	window_size = int(args.windowsize)
	step_size = int(args.stepsize)
	idx = str(args.idx).zfill(2) # should be 01-96

	# Filter have to be 'PASS'
	myfilter='PASS'
	contig_summary='/u/project/rwayne/snigenda/finwhale/scripts/config/minke_contig_summary.csv'
	# Read chromosome length dictionary
	chrom_size=read_minke_contig_len(contig_summary)
	# Read the contig list dictionary
	chrom_list=read_minke_contig_idx(contig_summary)
	# Use the one for this vcf file
	chrom_list=chrom_list[idx]
	# Check if the file is indexed and get the index
	parsevcf=index_invcf(filename)
	# Get sample names
	samples=get_samples(filename)
	# Write output header
	output = open(outprefix + '_het_%swin_%sstep.txt' % (window_size, step_size), 'w')
	# note that the `join` function does not join the first one, so you need to add 'calls_' before %s
	output.write('chrom\twindow_start\twindow_end\tsites_total\tcalls_%s\thets_%s\n' % ('\tcalls_'.join(samples), '\thets_'.join(samples)) )

	# Looping through the chromosomes in the chrom_list
	for chrom in chrom_list:
		# get start and end positions
		start_pos, end_pos = get_start_end(filename, chrom, chrom_size)
		# Initialize window start and end coordinates
		window_start = start_pos
		window_end = start_pos + window_size - 1
		# Calculate stats for window, update window start and end positions,
		# repeat to end of chromosome
		while window_start < end_pos:
			if window_end >= end_pos:
				window_end = end_pos
			snp_cal(parsevcf, samples, output, chrom, window_start, window_end, myfilter)
			window_start = window_start + step_size
			window_end = window_start + window_size - 1

	# Close files and exit
	parsevcf.close()
	output.close()

if __name__ == "__main__":
	sys.exit(main())
