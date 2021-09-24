'''
# Title: reformatting the INFO field for inputting GATK
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Tue Aug 18 20:39:05 2020
# @modification: Thu Nov  5 21:59:08 2020
# @modification: 1. Change to directly use vcf output from SIFT4G annotator
# @modification: 2. Remove the trailing '^M' from each line, otherwise not matching the gatk VCF specifications

'''

###########################################################
# import packages and input arguments
import sys

###########################################################
# define functions

###########################################################
# main
def main():
	vcf_file=sys.argv[1]
	# read input files
	with open(vcf_file, 'rt') as inVCF:
		for line0 in inVCF:
			# remove the before and trailing space, \r and \n
			line0=line0.strip()
			if line0.startswith('#'):
				if line0.startswith('##SIFT_Threshold'):
					line=line0.split(': ')
					sys.stdout.write('%s\n' % '='.join(line))
				else:
					sys.stdout.write('%s\n' % line0)
			else:
				line=line0.split('\t')
				# convert the INFO field
				if ' ' in line[7]:
					line[7]=line[7].replace(' ', '_')
					sys.stdout.write('%s\n' % '\t'.join(line))
				else:
					sys.stdout.write('%s\n' % line0)

if __name__ == "__main__":
	sys.exit(main())

