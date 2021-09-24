# -*- coding: utf-8 -*-
'''
# Title:  extract snpEff annoations from a vcf file
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Fri Apr 16 00:25:24 2021
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
    parser = argparse.ArgumentParser(description='This script will extract annotations for given vcf files')

    parser.add_argument(
            '--path', required=True,
            help='REQUIRED. Working directory')

    parser.add_argument(
            '--VCF', required=True,
            help='REQUIRED. Path to the VCF file. Should be gzipped.')

    parser.add_argument(
            '--filter', required=True,
            help='Site filter to apply. Separate by comma')

    parser.add_argument(
            '--verbose', required=False, action="store_true",
            help='Should it be verbose and print the lines with multiple annotations.')

    parser.add_argument(
            '--local', required=False, action="store_true",
            help='Should the script be run locally.')

    parser.add_argument(
            '--outfile', required=True,
            help='REQUIRED. Path to the output file.')

    args = parser.parse_args()
    return args

# load snpEff annotations sortorder and return a str list
def load_snpEff_sortorder(orderfile):
    with open(orderfile, 'r') as ff:
        snpEff_order = [line.strip() for line in ff]
    return snpEff_order

# input a list of ANN
# return a list with only one elemnt ANN that has the most deleterious outcome
# usually the first annotation should be the most deleterious outcome already
def get_mostdel_ANN(line0, ANN, snpEff_order, verbose):
    # most deleterious one should be reported according to: VCFannotationformat_v1.0.pdf
    # does not work with annotations with the same ranking
    ANN_Annotation = [ann.split('|')[1] for ann in ANN]
    # if there are multiple ANN_Annotation
    ANN_order = []
    for ii in ANN_Annotation:
        if '&' in ii:
            ii = ii.split('&')
            iiorder = min([snpEff_order.index(jj) for jj in ii])
        else:
            iiorder = snpEff_order.index(ii)
        ANN_order.append(iiorder)
    ANN_orderminid = [ii for ii, xx in enumerate(ANN_order) if xx == min(ANN_order)]
    ANN_minid = ANN_orderminid[0]
    # if there are multiple minimum
    if len(ANN_orderminid) != 1:
        # pick the first one anyways but generate a warning
        if verbose == True:
            errmsg=str('WARNING! Multiple snpEff annotations with same sort orders: ' + ','.join(ANN) + '\n' + 'WARNING! Picked ANN: Rank=' + str(ANN_minid) + ';Value=' + ANN[ANN_minid] + '\n')
            sys.stderr.write(errmsg)
    # the snpEff should have already ranked the rankings by the sort order defined
    if ANN_minid != 0:
        errmsg=str('ERROR! Multiple snpEff annotations and the first annotation is not the most deleterious: ' + line0)
        sys.stderr.write(errmsg)
        # sys.exit(1)
        # return None
    # pick the most deleterious ANN from the first one in the ranking
    outANN = [ANN[0]]
    return outANN



###########################################################
## def variables
header = '# CHROM\tSTART\tEND\n'

###########################################################
## main

def main():
    args = parse_args()
    os.chdir(args.path)
    outfile = open(args.outfile + '.txt', 'w')
    # outfile.write(header)
    myfilter = str(args.filter).split(',')
    # verbose
    if args.verbose:
        verbosity = True
    else:
        verbosity = False
    # load the annotation ranking file
    if args.local:
        snpEff_order = load_snpEff_sortorder('/Users/linmeixi/Lab/fin_whale/scripts_analyses/config/snpEff_annotation_sort_order.txt')
    else:
        snpEff_order = load_snpEff_sortorder('/u/project/rwayne/meixilin/fin_whale/analyses/scripts/config/snpEff_annotation_sort_order.txt')
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
                # only work with lines that contain ANN (annotations)
                if 'ANN' in myinfo:
                    # get lines you are interested in
                    mychr=line[0]
                    myendpos=line[1]
                    mystartpos=str(int(myendpos) - 1)
                    # instead of iterating through each one, make a dict.
                    infoFields=dict(s.split('=') for s in myinfo.split(';'))
                    # split by the ','
                    myANN=infoFields['ANN'].split(',')
                    # not allowing multiple mutations
                    if len(myANN) != 1:
                        if args.verbose:
                            # write before and after
                            errmsg=str('INFO! Multiple snpEff annotations: ' + line0)
                            sys.stderr.write(errmsg)
                        # start ranking the annotations and get the most deleterious one
                        myANN = get_mostdel_ANN(line0, myANN, snpEff_order, verbosity)
                        if args.verbose:
                            errmsg=str('INFO! Picked annotation: ' + ','.join(myANN) + '\n')
                            sys.stderr.write(errmsg)
                    ANNline = myANN[0].split('|')
                    ANN_ALT = str(ANNline[0]) # the alternative allele
                    ANN_Annotation = str(ANNline[1]) # second field
                    ANN_Type = str(ANNline[7]) # get the transcript biotype
                    # check that it matches the right ALT allele not restricting biotype
                    if ANN_ALT == str(line[4]):
                        # output only the annotation type
                        output = [mychr, myendpos, ANN_Annotation, ANN_Type]
                        outfile.write('\t'.join(output) + '\n')
                    else:
                        errmsg=str('ERROR! Wrong ALT allele in ANN: ' + line0)
                        sys.stderr.write(errmsg)
                        sys.exit(1)
    outfile.close()


if __name__ == '__main__':
    sys.exit(main())

