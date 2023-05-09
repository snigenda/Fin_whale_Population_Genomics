# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 10:13:40 2018
@author: annabelbeichman
"""
import sys
import os
import gzip
from collections import OrderedDict
import argparse
import datetime
from itertools import combinations
################# This function will count the number of monomorphci sites passing a projection level for minimum calls ################
# sites must be all 0/0 
# must have at least (projection / 2 ) genotype calls (opposed to missing data)
# note this script assumes unphased data with 0/0 calls rather than 0|0 calls
# This script modifies code from Overcast's EasySFS
# and will give monomorphic site counts passing given projections for each individual pop
# and for each pair of pops (passing sites meet projection levels specific to each population for bot pops in the pair)
# will output both sets of counts.
############### Parse input arguments ########################
parser = argparse.ArgumentParser(description='Count the number of 0/0 sites that have have called genotypes in at least [your projection value  / 2 ] or more individuals. Note that easy SFS projection values are in *haploids* Not Diploids')
parser.add_argument("--vcf",required=True,help="path to vcf file")
parser.add_argument("--popMap",required=True,help="file mapping individuals to populations")
parser.add_argument("--proj",required=True,help="**haploid** projection values for each population, separated by commas (for 8 diploids you would put in 16 as your projection value (consistent with easysfs)). e.g. 16,14,13,12")
parser.add_argument("--popIDs",required=True,help="the populations in the SAME order as you put the projection values, separated by commas. (tells us which proj goes to which population; these should be the same popIDs that are in your popMap. e.g. KUR,AK,CA,AL,COM")
parser.add_argument("--outdir",required=True,help="path to output directory")
parser.add_argument("--outPREFIX",required=False,help="output file prefix (optional)",default="")

args = parser.parse_args()
vcfFile=args.vcf
popMap=str(args.popMap)
outdir=str(args.outdir)
prefix=str(args.outPREFIX)
projectionValues= [int(x) for x in args.proj.split(",")]
popIDs = [str(x) for x in args.popIDs.split(",")]



if len(popIDs)==len(projectionValues):
    projDict=dict(zip(popIDs,projectionValues))
else:
    print("you must supply the same number of projection values as populations that you want to project!")
    sys.exit()

# write out projection values:
proj_outputFile=open(str(outdir)+"/projectionValues.txt","w")
proj_outputFile.write("population\tProjectionValueHaploids\n")
for key,value in projDict.items():
    proj_outputFile.write('{0}\t{1}\n'.format(key,value))
proj_outputFile.close()


# get sample names
def getSamplesFromVCF(filepath):
    inVCF = gzip.open(filepath, 'r')
    samples=[]
    for line in inVCF:
        if line.startswith('##'):
            pass
        else:
            for i in line.split()[9:]: samples.append(i)
            break
    inVCF.close()
    return samples
samples  = getSamplesFromVCF(vcfFile) # these are what are in the VCF file
######## function from easysfs:  ###############
#### get the population assignments for the samples from a pops_file  ####
def get_populations(pops_file, verbose=False):
    # Here we need to read in the individual population
    # assignments file and do this:
    # - populate the locs dictionary with each incoming population name
    # - populate another dictionary with individuals assigned to populations
    # Add the 'U' to handle opening files in universal mode, squashes the
    # windows/mac/linux newline issue.

    try:
        with open(pops_file, 'rU') as popsfile:
            ind2pop = {}
            pops = OrderedDict()
        
            lines = popsfile.readlines()
            ## Get all the populations
            for line in lines:
                pops.setdefault(line.split()[1], [])
        
            for line in lines:
                ind = line.split()[0]
                pop = line.split()[1]
                ind2pop[ind] = pop
                pops[pop].append(ind)

        print("Processing {} populations - {}".format(len( pops ), pops.keys()))
        if(verbose):
            for pop,ls in pops.items():
                print(pop, ls)

    except Exception as inst:
        msg = """
    Problem reading populations file. The file should be plain text with one
    individual name and one population name per line, separated by any amount of
    white space. There should be no header line in this file. 
    An example looks like this:
        ind1    pop1
        ind2    pop1
        ind3    pop2
        ind4    pop2"""
        print(msg)
        print("    File you specified is: ".format(pops_file))
        print("    Error - {}".format(inst))
        raise

    return ind2pop, pops # note this returns TWO things, so you need to give it two variables when you run it : see below
    
######## another function from easysfs to check the validity of the inputs and make the popMap and vcf sample list match ##################
def check_inputs(ind2pop, indnames, pops):
    ## Make sure all samples are present in both pops file and VCF, give the user the option
    ## to bail out if something is goofy
    pop_set = set(ind2pop.keys())
    vcf_set = set(indnames)
    
    if not pop_set == vcf_set:
        print("\nSamples in pops file not present in VCF: {}\n"\
            .format(", ".join(pop_set.difference(vcf_set))))
        ## Remove the offending samples from ind2pop
            # in this case "pop" is popping off the samples that are in the diff bet popset and vcfset
        map(ind2pop.pop, pop_set.difference(vcf_set))
        ## want to rewrite this to be safer.
        # **** BAD CODING **** don't iterate over list as you change it
        #### rewrite to be safer. and check the rest of his code!!! #####
        print("Samples in VCF not present in pops file: {}\n"\
            .format(", ".join(vcf_set.difference(pop_set))))

        ## Remove the offending individuals from the pops dict
        for k,v in pops.items():
            for ind in pop_set.difference(vcf_set):
                # try a safer change:
                if ind in v:
                    v.remove(ind) # doesn't return anything, just removes the entry from v; ab modifie to have this happen in two steps (correct)
                    pops[k] = v
        for k, v in pops.items():
            if not v:
                print("Empty population, removing - {}\n".format(k))
                pops.pop(k)
        for key,value in pops.items():
            print(str(len(value)) + ' Surviving individuals for {0}: {1}\n'.format(key,value))
        # AB: 20181121: removing this part of the script so that there isn't an interactive portion
        #cont = raw_input("\nContinue, excluding samples not in both pops file and VCF? (yes/no)\n")
        #while not cont in ["yes", "no"]:
            #cont = raw_input("\nContinue, excluding samples not in both pops file and VCF? (yes/no)\n")
        #if cont == "no":
            #sys.exit()
######## something is off with this code -- why does it work in easysfs? or does it? ############
    return ind2pop, indnames, pops
######################################################
    
    

ind2pop, pops = get_populations(popMap) # this works! gives you
# ind2pop which is 1:1 mapping of inds : pop and pops which is collections of each pop and is a dictionary (useful)

ind2pop, samples, pops = check_inputs(ind2pop, samples, pops)
############# THIS ISN"T WORKING -- WHY????? ############
# set up empty dictionary to keep track of sites that pass for each population

################# read through VCF #########
################ SET UP PER SAMPLE DICTIONARY OF GENOTYPE CALLS ########
# read through vcf, skip the header lines
# input: pops file resulting from check_inputs; your projection value; your vcf file
def count_PassingMonomorphicSites(pops,projDict,VCF):
    print("beginning counts")
    inVCF = gzip.open(VCF, 'r')
    inVCF.seek(0)
    # set up dictionary:
    countDict=dict()
    #TotalDict=dict()
    # set up 2D dictionary for monomorphic site that pass each population in the pair's projection level
    pairDict=dict()
    # combos of every pair of populations (for 2D)
    combos=list(combinations(pops.keys(),2)) # every pair of populations
    # populate dicts with zeroes:    
    for population in pops.keys():
        countDict[population]=0
        #TotalDict[population]=0
    for popPair in combos:
        pairDict[popPair]=0
    # skip header lines
    for line0 in inVCF:
        if line0.startswith('#'):
            continue
        ### For all other non-header lines, set up a dicionary of genotype calls
        line=line0.strip().split('\t') # this splits line by tabs
        #CHROM0    POS1    ID2    REF3    ALT4    QUAL5    FILTER    INFO    FORMAT    [indivudals]
        mygenoinfo=line[9:]
        allCalls=[i.split(":")[0] for i in mygenoinfo] # get genotype calls
        # want to exclude variable sites AT THIS STAGE because any variable site 
        # made it into neutral.snps file and is taken care of by EasySFs even if it's monomorphic within a population
        if "0/1" in allCalls or "1/0" in allCalls or "1|0" in allCalls or "0|1" in allCalls or "1|1" in allCalls or "1/1" in allCalls:
            continue
        else: # it's entirely monomorphic
            callDict = dict(zip(samples,allCalls))
            # set up a dict that keeps track of if a particular population passes or fails at this site:
            # reset passDict each site: 
            passDict=dict()
            #for population in pops.keys():
            #    passDict[population]="NA"
            for population in pops.keys(): # go through each population 
                    pop_gts = [ callDict[x] for x in pops[population] ] 
                    popProjValue=projDict[population]
                    # this is saying to go through each individual of the population and get that call eg. callDict["145_Elut_CA_145"] is 0/0 (fake example) and pops[pop] gives you all individuals from the pops dictionary (from the popmap file) 
                    # check if the number of 0/0 for that population is greater than the projection value / 2 (divided by 2 because projection of 16 represents 8 diploids and these are diploid genotypes)
                    # have already selected for only monomorphic sites above
                    # passes: 
                    if pop_gts.count("0/0") >= (float(popProjValue)/2):
                        countDict[population] += 1
                        passDict[population]="PASS"
                    # fails:
                    else:
                        passDict[population]="FAIL"
                        continue
            # Get sites that pass projection values for each pair of populations
            for popPair in combos:
                pop1=popPair[0]
                pop2=popPair[1]
                # if both pops passed their projection at this site (even if at different projection levels, that's okay) then add one to their shared monomorphic site bin. if not, then continue.
                if passDict[pop1]=="PASS" and passDict[pop2]=="PASS":
                    pairDict[popPair]+=1
                else:
                    continue
                
    #print(str(TotalDict))
    inVCF.close()
    return countDict, pairDict
    
singlePopcounts, popPaircounts, = count_PassingMonomorphicSites(pops,projDict,vcfFile)

############ write out counts 
#todaysdate=datetime.datetime.today().strftime('%Y%m%d')

# write out single pop counts:
singlePop_outputFile=open(str(outdir)+"/countsOfMonomorphicSites.perPop."+str(prefix)+".txt","w")
singlePop_outputFile.write("population\tHomREFcountPassingProjThreshold\tProjectionValue\n")
for key,value in singlePopcounts.items():
    singlePop_outputFile.write('{0}\t{1}\t'.format(key,value))
    singlePop_outputFile.write(str(projDict[key])+"\n")
singlePop_outputFile.close()

# write out dual pop counts
dualPop_outputFile=open(str(outdir)+"/countsOfMonomorphicPassingProjectionThresholds.perPair."+str(prefix)+".txt","w")
dualPop_outputFile.write("population1\tpopulation2\tHomREFcountPassingBothProjThresholds\tProjectionValue1\tProjectionValue2\n")
for key,value in popPaircounts.items():
    # split key into the two populations that make it up (key[0] and key[1])
    dualPop_outputFile.write('{0}\t{1}\t{2}\t'.format(key[0],key[1],value))
    # projection value of first population: 
    dualPop_outputFile.write(str(projDict[key[0]])+"\t")
    dualPop_outputFile.write(str(projDict[key[1]])+"\n")
dualPop_outputFile.close()
#[outputFile.write(('{0}\t{1}\t'.format(key,value)+str(projDict[key])) for key,value in counts.items())]
#outputFile.close()

sys.exit()
