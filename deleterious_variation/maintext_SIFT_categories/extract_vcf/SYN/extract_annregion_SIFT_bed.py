# -*- coding: utf-8 -*-
'''
# Title:  extract SIFT annoations regions depending on the annotation types
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Dec 18 00:56:25 2020
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
def parse_args():
	'''
	Parse command-line arguments
	'''
	parser = argparse.ArgumentParser(description='This script will extract bedfiles for given vcf files')

	parser.add_argument(
			'--VCF', required=True,
			help='REQUIRED. Path to the VCF file. Should be gzipped.')

	parser.add_argument(
			'--anntype_SIFT', required=True,
			help='REQUIRED. The SIFT annotations to be extracted. Only one type of annotation per time.')

	parser.add_argument(
			'--filter', required=True,
			help='Site filter to apply. Separate by comma')

	parser.add_argument(
			'--verbose', required=False, action="store_true",
			help='Should it be verbose and print the lines with multiple annotations.')

	parser.add_argument(
			'--outprefix', required=True,
			help='REQUIRED. Path and prefix to the output file.')

	args = parser.parse_args()
	return args

# load SIFT annotations sortorder and return a str list
def load_SIFT_sortorder(orderfile):
	with open(orderfile, 'r') as ff:
		SIFT_order = [line.strip() for line in ff]
	return SIFT_order

# input a list of ANN
# return a list with only one elemnt ANN that has the most deleterious outcome
# first ranked by annotation type second ranked by the SIFT score
def get_mostdel_ANN(line0, ANN, SIFT_order, verbose):
	# most deleterious one are adapted according to: VCFannotationformat_v1.0.pdf
	# for SIFT annotations with the same ranking, the most deleterious one is selected by the SIFT score
	ANN_Annotation = [ann.split('|')[5] for ann in ANN]
	ANN_Score = [ann.split('|')[8] for ann in ANN]
	# assign an impossible value of SIFT score if that field does not exist (NA or empty)
	ANN_Score = ['2.00' if score in ('', 'NA') else score for score in ANN_Score]
	# if there are multiple ANN_Annotation
	ANN_order = [SIFT_order.index(ii) for ii in ANN_Annotation]
	ANN_orderminid = [ii for ii, xx in enumerate(ANN_order) if xx == min(ANN_order)]
	if len(ANN_orderminid) == 1:
		ANN_minid = int(ANN_orderminid[0])
	else:
		# access the ANN_Score by the annotation of interest and zip it with the index in list ANN
		# sort the zipped object by the ANN_Score and access the most deleterious one in the first element of the zipped object
		ANN_minscore = sorted(zip([float(ANN_Score[ii]) for ii in ANN_orderminid], ANN_orderminid))
		# pick the most deleterious ANN from the first one in the ranking (or just choose the first out of all the same SIFT score and type of annotations)
		ANN_minid = int(ANN_minscore[0][1])
		# Add warning if can't distinguish the deleterious rankings
		if ANN_minscore[0][0] == ANN_minscore[1][0]:
			if verbose == True:
				errmsg=str('WARNING! Multiple SIFT annotations with same sort orders: ' + ','.join(ANN) + '\n' + 'WARNING! Picked ANN: Rank=' + str(ANN_minid) + ';Value=' + ANN[ANN_minid] + '\n')
				sys.stderr.write(errmsg)
	outANN = [ANN[ANN_minid]]
	if len(outANN) == 1:
		return outANN
	else:
		errmsg=str("ERROR! Something wrong getting most deleterious ANN at " + line0)
		sys.stderr.write(errmsg)
		sys.exit(1)
		return None



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
	# verbose
	if args.verbose:
		verbosity = True
	else:
		verbosity = False
	# load the annotation ranking file
	SIFT_order = load_SIFT_sortorder('/u/project/rwayne/meixilin/fin_whale/analyses/scripts/config/SIFT_annotation_sort_order.txt')
	with gzip.open(args.VCF, 'rt') as VCF:
		for line0 in VCF:
			if not line0.startswith('#'):
				line = line0.strip().split('\t')
				myinfo = line[7]
				# check if pass filter, i.e if the elements in myfilterinfo is a complete subset of myfilter
				# don't split by ';' check the entire string
				# only accept filters including: PASS,WARN_missing,WARN_excessHet
				myfilterinfo = [line[6]]
				if not set(myfilterinfo).issubset(set(myfilter)):
					continue
				# only work with lines that contain SIFTINFO (SIFT annotations)
				if 'SIFTINFO' in myinfo:
					# get lines you are interested in
					mychr=line[0]
					myendpos=line[1]
					mystartpos=str(int(myendpos) - 1)
					# instead of iterating through each one, make a dict.
					infoFields=dict(s.split('=') for s in myinfo.split(';'))
					# split by the ','
					myANN=infoFields['SIFTINFO'].split(',')
					# not allowing multiple mutations
					if len(myANN) != 1:
						if args.verbose:
							# write before and after
							errmsg=str('INFO! Multiple SIFT annotations: ' + line0)
							sys.stderr.write(errmsg)
						# start ranking the annotations and get the most deleterious one
						myANN = get_mostdel_ANN(line0, myANN, SIFT_order, verbosity)
						if args.verbose:
							errmsg=str('INFO! Picked annotation: ' + ','.join(myANN) + '\n')
							sys.stderr.write(errmsg)
					ANNline = myANN[0].split('|')
					ANN_ALT = str(ANNline[0]) # the alternative allele
					ANN_Annotation = str(ANNline[5]) # second field
					ANN_Type = str(ANNline[4]) # get the transcript biotype
					# if matches the target annotation type
					if ANN_Annotation == str(args.anntype_SIFT):
						# check that it matches the right ALT allele and biotype
						if ANN_ALT == str(line[4]) and ANN_Type == 'CDS':
							# output the bedfile
							output = [mychr, mystartpos, myendpos]
							outfile.write('\t'.join(output) + '\n')
						else:
							errmsg=str('ERROR! Wrong biotype or ALT allele in SIFTINFO: ' + line0)
							sys.stderr.write(errmsg)
							sys.exit(1)
	outfile.close()


if __name__ == '__main__':
	sys.exit(main())

