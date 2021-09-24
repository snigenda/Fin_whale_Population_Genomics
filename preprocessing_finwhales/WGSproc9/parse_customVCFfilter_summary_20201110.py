# -*- coding: utf-8 -*-
'''
Title: parse customVCFfilter_summary output line by line
Author: Meixi Lin (meixilin@ucla.edu)
Date: Tue Nov 10 19:47:53 2020
'''

##########################################################
## import packages
import argparse
import os
import csv
import pandas
import sys
import gzip

##########################################################
## def functions

##########################################################
## def variables

##########################################################
## main
def main():
	# parse input
	parser=argparse.ArgumentParser(description='Count number of various filter reasons in a VCF summary files')
	parser.add_argument("--dir",required=True,help="Working directory")
	parser.add_argument("--prefix",required=True,help="Prefix to the file")
	parser.add_argument("--contig",required=True,help="Contig id to append")
	parser.add_argument("--outdir",required=True,help="Output directory")

	args=parser.parse_args()
	os.chdir(str(args.dir))
	myprefix=str(args.prefix)
	contigid=str(args.contig)
	outdir=str(args.outdir)
	# open the file
	infile=gzip.open(str(myprefix+".tsv.gz"), 'rt')

	# get sample names
	first_line=infile.readline().strip()
	samples=first_line.split('\t')[3:]

	# initialize
	filterdict={}
	# make empty list of dicts
	gtfilterdict=[{} for ii in range(len(samples))]
	# name the empty lists of dicts by the samples (dict of dicts)
	gtfilterdict=dict(zip(samples, gtfilterdict))

	# count filters from the second line
	for line0 in infile:
		line=line0.strip().split('\t')
		if line[2] not in filterdict:
			filterdict[line[2]]=1
		else:
			filterdict[line[2]]+=1
		# now start genotype filters
		for ii in range(len(samples)):
			if line[ii+3] not in gtfilterdict[samples[ii]]:
				gtfilterdict[samples[ii]][line[ii+3]]=1
			else:
				gtfilterdict[samples[ii]][line[ii+3]]+=1
	infile.close()

	# start writing to the text for the site filters
	with open(outdir+"/"+myprefix+"_FILTER_tally.csv", "w") as ff:
		outf_filter=csv.writer(ff)
		outf_filter.writerow(['SITE_FILTER', str('N_FILTER_'+contigid)])
		for key, val in filterdict.items():
			outf_filter.writerow([key, val])

	# start writing to the text for the gt filters
	gtfilterdf = pandas.DataFrame(gtfilterdict)
	gtfilterdf.to_csv(str(outdir+"/"+myprefix+"_GTFILTER_tally.csv"))


if __name__ == "__main__":
	sys.exit(main())
