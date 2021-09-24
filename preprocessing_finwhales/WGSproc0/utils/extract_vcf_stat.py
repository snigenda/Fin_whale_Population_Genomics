# Title: Quickly extract information from vcf files for variant filter 
# Date: Sun Apr  5 15:12:08 2020
# @modification Wed Nov  4 10:45:09 2020
# @modification DEFUNCT, DO NOT USE. Can be replaced by the gatk VariantsToTable function
# Example:
# python extract_vcf_stat.py --VCF Minke.chr01.ann.vcf.gz --outfile "haha/chr01"

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
    parser = argparse.ArgumentParser(description="This script will extract CHROM, QUAL, DP, QD and genotpe DP for all variant sites")

    parser.add_argument(
            "--VCF", required=True,
            help="REQUIRED. Path to the VCF file. Should be gzipped.")

    parser.add_argument(
            "--outfile", required=True,
            help="REQUIRED. Path to the output file.")

    args = parser.parse_args()
    return args

# get genotype Allelic balance
# input should be a dictionary
def get_gtAB(GT_entry):
    myDP = GT_entry.get("DP", "NA")
    myAD = GT_entry.get("AD", "NA")
    badDP = ["0", "NA", "."]
    badAD = ["NA", "."]
    # if genotype depth and allelic depth missing, skip this  
    if myDP in badDP or myAD in badAD: 
        AB = "NA"
    else:
        REF=float(myAD.split(',')[0])
        DP=float(myDP)
        AB=float(REF/DP)
        AB = "{:.3f}".format(AB)
    return AB

# get genotype GQ or reference RGQ
# input should be a dictionary
def get_gtGQ(GT_entry):
    myGQ = GT_entry.get("GQ", "NA")
    myRGQ = GT_entry.get("RGQ", "NA")
    if myGQ == "NA":
        return myRGQ # NA or the actual "RGQ"
    else:
        return myGQ


###########################################################
## def variables 
header_info = "CHROM\tPOS\tQUAL\tDP\tQD\n"
header_gt = "\tCHROM\tPOS\n"

###########################################################
## main 
def main():
    args = parse_args()
    out_info = args.outfile + "_QUAL_DP_QD.tsv"
    out_gt = args.outfile + "_gtVAL.tsv"
    out_dp = args.outfile + "_gtDP.tsv"
    out_gq = args.outfile + "_gtGQ.tsv"
    out_ab = args.outfile + "_gtAB.tsv"
    outinfo = open(out_info, "w")
    outgt = open(out_gt, "w")
    outdp = open(out_dp, "w")
    outgq = open(out_gq, "w")
    outab = open(out_ab, "w")

    with gzip.open(args.VCF, "r") as VCF:
        # Get list of samples
        samples=[]
        for line in VCF:
            if line.startswith('##'):
                pass
            else:
                for i in line.split()[9:]: samples.append(i)
                break
        # write header for output 
        outinfo.write(header_info)
        outgt.write("\t".join(str(x) for x in (samples)))
        outgt.write(header_gt)
        outdp.write("\t".join(str(x) for x in (samples)))
        outdp.write(header_gt)
        outgq.write("\t".join(str(x) for x in (samples)))
        outgq.write(header_gt)
        outab.write("\t".join(str(x) for x in (samples)))
        outab.write(header_gt)

        # Get back to line one 
        VCF.seek(0)
        for line in VCF:
            if not line.startswith("#"):
                line = line.rstrip("\n")
                line = line.split("\t")
                myqual=line[5]
                info_col = line[7]

                # split info fields:
                # instead of iterating through each one, make a dict.
                if info_col == '.':
                    myDP = "NA"
                    myQD = "NA"
                else: 
                    infoFields=dict(s.split('=') for s in info_col.split(";"))
                    myDP=infoFields.get("DP", "NA")
                    myQD=infoFields.get("QD", "NA")
                    
                outinfo.write("\t".join(str(x) for x in (line[0], line[1], myqual, myDP, myQD)))
                outinfo.write("\n")

                # now read for genotype, genotype depth, genotype quality and allelic depth 
                GT_format = line[8].split(':')
                gtVAL=[]
                gtDP=[]
                gtGQ=[]
                gtAB=[]
                for ii in range(0,len(samples)):
                    GT_entry = dict(zip(GT_format, line[ii+9].split(':')))
                    gtVAL.append(GT_entry.get("GT", "NA"))
                    gtDP.append(GT_entry.get("DP", "NA"))
                    gtGQ.append(get_gtGQ(GT_entry))
                    gtAB.append(get_gtAB(GT_entry))
                gtVAL.extend(line[:2])
                gtDP.extend(line[:2])
                gtGQ.extend(line[:2])
                gtAB.extend(line[:2])
                # some of the DP could contain "."
                # write the output 
                outgt.write("\t".join(str(x) for x in (gtVAL)))
                outgt.write("\n")
                outdp.write("\t".join(str(x) for x in (gtDP)))
                outdp.write("\n")
                outgq.write("\t".join(str(x) for x in (gtGQ)))
                outgq.write("\n")
                outab.write("\t".join(str(x) for x in (gtAB)))
                outab.write("\n")
    VCF.close()

    # close other files 
    outinfo.close()
    outgt.close()
    outdp.close()
    outgq.close()
    outab.close()

sys.exit(main())

