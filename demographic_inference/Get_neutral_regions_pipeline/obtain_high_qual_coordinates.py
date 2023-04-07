import os
import sys
import argparse
import gzip

def parse_args():
	"""
	Parse command-line arguments
	"""
	parser = argparse.ArgumentParser(description="This script takes in a VCF after step10_obtain_high_qual_sites and returns a list of coordinates that is high quality. The reason for this extra script is that GATK variant filtration will annotate a site as PASS when there is no information in the INFO column. When there is no information in the INFO column, it is simply a dot.")

	parser.add_argument(
            "--VCF", required=True,
            help="REQUIRED. VCF file. Should be gzipped")

	parser.add_argument("--outfile", required=True, 
			help="REQUIRED. Name of output file. End of .bed. Note that this bed file is 0-based.")

	args = parser.parse_args()
	return args

def main():
	args = parse_args()
	outfile = open(args.outfile, "w")

	with gzip.open(args.VCF, "r") as file:
		for line in file:
			if not line.startswith("#"):
				line = line.rstrip("\n")
				line = line.split("\t")
				if line[7] != ".":
					if line[6] == "PASS":
						zero_based_coord = int(line[1])-1
						toprint = [line[0], str(zero_based_coord), line[1]] # this is because vcf is 1-based while bed is 0-based
						print >>outfile, "\t".join(x for x in toprint)
main()
