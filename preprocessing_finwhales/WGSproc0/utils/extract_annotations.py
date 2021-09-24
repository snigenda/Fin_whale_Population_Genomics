# Title:  extract snpEff annoation errors and combine with the previous filters
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Sat Mar 28 14:31:57 2020

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
    """
    Parse command-line arguments
    """
    parser = argparse.ArgumentParser(description="This script will extract annotations for vcf files")

    parser.add_argument(
            "--path", required=True,
            help="REQUIRED. Working directory")

    parser.add_argument(
            "--VCF", required=True,
            help="REQUIRED. Path to the VCF file. Should be gzipped.")

    # parser.add_argument(
    #         "--scaffold", required=True,
    #         help="REQUIRED. Full name of scaffold")

    parser.add_argument(
            "--outfile", required=True,
            help="REQUIRED. Path to the output file.")

    args = parser.parse_args()
    return args

###########################################################
## def variables 
header = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tANN_Allele\tANN_Annotation\tANN_Annotation_Impact\tANN_Gene_Name\tANN_Gene_ID\tANN_Feature_Type\tANN_Feature_ID\tANN_Transcript_BioType\tANN_Rank\tANN_HGVS.c\tANN_HGVS.p\tANN_cDNA.pos / cDNA.length\tANN_CDS.pos / CDS.length\tANN_AA.pos / AA.length\tANN_Distance\tANN_ERRORS / WARNINGS / INFO\n"

###########################################################
## main 

def main():
    args = parse_args()
    os.chdir(args.path)
    outfile = open(args.outfile, "w")
    # scaff=str(args.scaffold)
    outfile.write(header) 
    with gzip.open(args.VCF, "r") as VCF:
    #with gzip.open(filepath,"r") as VCF:
        for line in VCF:
            if not line.startswith("#"):
                line = line.rstrip("\n")
                line = line.split("\t")
                # scaffold= line[0]
                info_col = line[7]
                # # if the scaffold isn't the scaff you've chosen, skip rest of loop
                # if scaffold!=scaff:
                #    continue
                # only work with lines that contain ANN (annotations)         
                if "ANN" in info_col: 
                    myscaf=line[0:7] # get all the previous lines
                    myinfo=line[7]
                    # split info fields:
                    # instead of iterating through each one, make a dict.
                    infoFields=dict(s.split('=') for s in myinfo.split(";"))
                    myANN=infoFields["ANN"]
                    # first split by "," to get alternate alleles
                    myANN=myANN.split(",")
                    for ann in myANN: 
                        ANNline = ann.split("|")
                        # append the scaffold information 
                        output = myscaf + ANNline
                        outfile.write("\t".join(output))
                        outfile.write("\n")
                else:
                    continue
    VCF.close()
    print("closing vcf now")
    outfile.close()
    print("closing outfile now")
sys.exit(main())

