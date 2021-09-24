# -*- coding: utf-8 -*-
'''
# Title: generated parametric bootstrap files
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Sun Feb 21 23:34:31 2021
# Example usage: python fsc26_generate_bootpar.py 1D.1Epoch_maxL.par 1D.1Epoch_boot.par
# @modification: Sun Apr 25 18:32:34 2021
# @modification: Update the nchr to match the neutral regions better
'''

###########################################################
## import packages
import sys

###########################################################
## def functions

###########################################################
## def variables

###########################################################
## main
def main():
	inFile=str(sys.argv[1])
	outFile=str(sys.argv[2])
	nchr=str(sys.argv[3]) # number of chromosome: Total length (L) used in dadi = 392,707,916 (ENP), divided by the loci length of 100 bp. Note that different population will have different L value used in dadi because of differences in the counting of missing sites
	nloci='100' # loci length (default to 100 bp)
	prevline='FIRSTLINE' # some random characters

	with open(outFile, 'w') as outfile:
		with open(inFile, 'rt') as infile:
			for line in infile:
				if prevline.startswith('//Number of independent loci [chromosome]'):
					outline=' '.join([nchr, '0'])+'\n'
				elif line.startswith('FREQ '):
					linefields=line.split(' ')
					outfields=['DNA', nloci] + linefields[2:]
					outline=' '.join(outfields)
				else:
					outline=line
				# update previous line
				prevline=line
				# output new lines
				outfile.write(outline)


if __name__ == "__main__":
	sys.exit(main())