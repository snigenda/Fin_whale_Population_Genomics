# Title: Quickly extract total genotype depth and total sites from vcf files for variant filter
# Date: Sun Apr  5 15:12:08 2020
# @modification Tue Nov  3 21:45:20 2020
# @modification NOTE: bcftools stats should work better and also faster see "fin_whale/scripts_analyses/Summary_stats/count_sites/wrapper_bcftools_stats_20200829.sh". The output from bcftools stats and this script had been checked and they were the same. This script is still being use because it defines the genotype depth and total called sites explicitly and easier to be calculated across multiple vcf files.
# Example:
# python extract_vcf_gtdp.py --VCF Minke.chr01.ann.vcf.gz --outfile "variant_summary/summary_stats/total_dp"

###########################################################
## import packages
import os
import sys
import argparse
import csv
import gzip
import numpy

###########################################################
## def functions
def parse_args():
	"""
	Parse command-line arguments
	"""
	parser = argparse.ArgumentParser(description="This script will extract genotype depth at all sites")

	parser.add_argument(
			"--VCF", required=True,
			help="REQUIRED. Path to the VCF file. Should be gzipped.")
	parser.add_argument(
			"--contiglist", required=True,
			help="REQUIRED. contiglist used in this vcf file.")
	parser.add_argument(
			"--outfile", required=True,
			help="REQUIRED. Path to the output file.")
	args = parser.parse_args()
	return args

def tally_gtdp(GT_entry):
	"""
	Count up genotype depth for each iteration
	"""
	myDP = GT_entry.get("DP", "NA")
	badDP = ["NA", "."]
	# if genotype depth and allelic depth missing, skip this
	if myDP in badDP:
		gtDP = 0
		isgt = 0 # not counting as a genotype
	else:
		gtDP = int(myDP)
		isgt = 1
	return gtDP, isgt

###########################################################
## def variables
header_gt = "\tLISTID\n"

###########################################################
## main
def main():
	args = parse_args()
	out_dp = args.outfile + "_sumDP.tsv"
	outdp = open(out_dp, "w")
	out_gt = args.outfile + "_countgt.tsv"
	outgt = open(out_gt, "w")
	listid = str(args.contiglist)
	with gzip.open(args.VCF, "r") as VCF:
		# Get list of samples =========
		samples=[]
		for line in VCF:
			if line.startswith('##'):
				pass
			else:
				for i in line.split()[9:]: samples.append(i)
				break
		# write header for output
		outdp.write("\t".join(str(x) for x in (samples)))
		outdp.write(header_gt)
		outgt.write("\t".join(str(x) for x in (samples)))
		outgt.write(header_gt)

		# Get back to line one =========
		VCF.seek(0)
		sumDP=[0 for ii in range(0,len(samples))] # sum of genotype depth
		countgt=[0 for ii in range(0,len(samples))] # count of genotype that accounts to the sum
		for line in VCF:
			if not line.startswith("#"):
				line = line.rstrip("\n")
				line = line.split("\t")

				GT_format = line[8].split(':')
				thisDP=[] # for this line
				thisgt=[] # for this line
				for ii in range(0,len(samples)):
					GT_entry = dict(zip(GT_format, line[ii+9].split(':')))
					thisDP.append(tally_gtdp(GT_entry)[0])
					thisgt.append(tally_gtdp(GT_entry)[1])
				sumDP=numpy.add(sumDP, thisDP).tolist()
				countgt=numpy.add(countgt, thisgt).tolist()

		# write the output =========
		outdp.write('\t'.join(str(x) for x in (sumDP)))
		outdp.write('\t' + listid + '\n')
		outgt.write('\t'.join(str(x) for x in (countgt)))
		outgt.write('\t' + listid + '\n')
	VCF.close()

	# close other files
	outdp.close()
	outgt.close()

if __name__ == "__main__":
	sys.exit(main())