# -*- coding: utf-8 -*-
'''
# Title: convert the gaps in transcripts to N
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Sun Jan 17 10:34:28 2021
# Example usage: python convert_gap2N.py <input> <output>
'''

###########################################################
## import packages
import sys
import argparse

###########################################################
## def functions

###########################################################
## def variables
def parse_args():
	'''
	Parse command-line arguments
	'''
	parser = argparse.ArgumentParser(description="Options for convert gaps to N")

	parser.add_argument(
			"--infasta", required=True,
			help="REQUIRED. Path to the fasta")

	parser.add_argument(
			"--outfasta", required=True,
			help="REQUIRED. Path to the output fasta")
	args = parser.parse_args()
	return args

###########################################################
## main
def main():
	args=parse_args()
	outFile=open(args.outfasta, 'w')
	with open(args.infasta, 'rt') as inFile:
		for line in inFile:
			if line[0]=='>':
				outFile.write(line)
			else:
				line=line.replace('-', 'N')
				outFile.write(line)
	outFile.close()

if __name__ == "__main__":
	sys.exit(main())
