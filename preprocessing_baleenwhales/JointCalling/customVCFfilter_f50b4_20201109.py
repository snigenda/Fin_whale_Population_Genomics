# -*- coding: utf-8 -*-
'''
Custom script for site- and genotype-level filtering of a VCF file
Author: Jacqueline Robinson
Modified by Meixi Lin.
Modification: 2020/11/09
1. Fix the maxD[sample] bug in gtfiltertype function. Now change the maxD to more specific maxDPDict and maxD=maxDPDict[sample]
2. Fix the maxD float datatype bug in reading customVCFfilter files
3. Change site filter variable from `filter` to `sitefilter` to avoid overlapping python builtin function filter()
4. Trim trailing space
5. Add functions `check_samples` and `get_samples`
6. Add functions `split_info` and `combine_info` to tolerate INFO fields of Type=Flag

Input_1 = VCF file (unphased)
Input_2 = Table with maximum depth for each individuals
Input_3 = Output filter statistics filename (should be in scratch, large text file)

Output_1 = filtered VCF file (pipe to stdout)
Output_2 = site and genotype filters applied (writes to Input_3 specified output file)
- Filtered sites are marked as FAIL_[?] in the 7th (FILTER) column
- Sites that pass go on to genotype filtering
- Genotypes that fail filters are changed to './.'

Example usage:

python customVCFfilter.py myJointVCFfile.vcf.gz maxD_file out_gtfilter | \
bgzip > myJointVCFfile_filtered.vcf.gz
tabix -p vcf myJointVCFfile_filtered.vcf.gz

NOTE: Will recalculate AC, AF, AN in INFO field.

NOTE on the script's activity:
1. get list of samples
2. add new headers for filter names in the VCF file
3. loop through each line of record and look for:
	1. sitefilters already applied in GATKfilter --> append sitefilters
	2. reference N --> FAIL_refN --> continue;
	3. reference not A/T/C/G --> FAIL_badRef
	4. alternate not A/T/C/G/. --> FAIL_badAlt
	5. no valid QUAL --> FAIL_noQUAL
	6. no valid INFO --> FAIL_noINFO
	7. fail if not monomorphic or simple SNP or no VariantType annotation --> FAIL_mutType
	8. In the GT field:
		1. No genotype_AD information --> FAIL_noADi
		2. No genotype_DP information --> FAIL_noDPi
		3. No genotype_GQ/RGQ information --> FAIL_noGQi
	If any sitefilter exists --> continue;
	9. Perform individual sample GT filter (using function GTfilter):
		Passing genotypes have to suffice ALL of these conditions
		1. Biallelic allele: ('0/0','0/1','1/1')
		2. GQ is not '.' and GQ >=20
		3. DP is not '.' and 8 <= DP <= maxD
		4. For genotype '0/1', 0.2<= ABHet <=0.8; for genotype '0/0' ABHet > 0.9; for genotype '1/1' ABHet < 0.1
		Otherwise, convert GT field to missing './.' --> next record
	10. Finish filtering, recalculate the fields.
		1. If REF + ALT == 0 after filtering --> FAIL_noGT --> continue;
		2. Recalculate: AC, AN, AF fields
	11. Perform all sample further site filter based on the modified GT:
		1. Filter out sites with more than 75% heterozygous genotypes --> WARN_excessHet
		2. Filter out sites with more than 20% missing genotypes --> WARN_missing


About continue: if the filters suffice before continue command, this script will finish up the analyses on this record, and continue to the next record (e.g. maybe it has FAIL_noGQi and FAIL_refN, but because FAIL_refN supercedes the FAIL_noGQi, the latter sitefilter will not show; the GT filters and AC, AN, AF statistics will not be analyzed if any of the FAIL_[?] sitefilters are applied earlier). if there is no continue, the script will keep analyzing (e.g. site filters will all add up before GTfilter analyses)
'''

###########################################################
# import packages and input arguments
import sys
import gzip
import re
import csv

vcf_file=sys.argv[1]
maxD_file=sys.argv[2]
out_file=sys.argv[3]

###########################################################
# define functions

# Filter to be applied to individual genotypes
### sample is the sample name
### GT_entry is the entire genotype entry for that individual
### minD is the minimum genotype depth in float
### maxDPDict is the maximum genotype depth dictionary in float
### maxDPDict[sample] is the maximum genotype depth (maxD) in float
### ADpos is the position of the AD in FORMAT (typically GT:AD:DP:GQ)
### DPpos is the position of the DP in FORMAT
### GQpos is the position of the GQ/RGQ in FORMAT
def GTfilter(sample, GT_entry, minD, maxDPDict, ADpos, DPpos, GQpos):
	if GT_entry[:1]=='.' : return GT_entry
	else:
		gt=GT_entry.split(':')
		GT_fail=str('./.:' + ':'.join(gt[1:]))
		if gt[0] in ('0/0','0/1','1/1') and gt[GQpos]!='.' and gt[DPpos]!='.':
			DP=int(gt[DPpos])
			GQ=float(gt[GQpos])
			minD=float(minD)
			maxD=float(maxDPDict[sample])
			if GQ>=20.0 and minD<=DP<=maxD:
				REF=float(gt[ADpos].split(',')[0])
				AB=float(REF/DP)
				if gt[0]=='0/0':
					if AB>=0.9: return GT_entry
					else: return GT_fail
				elif gt[0]=='0/1':
					if 0.2<=AB<=0.8: return GT_entry
					else: return GT_fail
				elif gt[0]=='1/1':
					if AB<=0.1: return GT_entry
					else: return GT_fail
				else: return GT_fail
			else: return GT_fail
		else: return GT_fail

# Based on the GTfilter function above, get the type of filters applied
def GTfiltertype(sample, GT_entry, minD, maxDPDict, ADpos, DPpos, GQpos):
	gtfilter=[]
	if GT_entry[:1]=='.' :
		gtfilter.append('gtmiss_before')
	else:
		gt=GT_entry.split(':')
		if not gt[0] in ('0/0','0/1','1/1'):
			gtfilter.append('gttype_badGT')
		if gt[GQpos]=='.':
			gtfilter.append('gttype_missingGQ')
		if gt[DPpos]=='.':
			gtfilter.append('gttype_missingDP')
		if gt[0] in ('0/0','0/1','1/1') and gt[GQpos]!='.' and gt[DPpos]!='.':
			DP=int(gt[DPpos])
			GQ=float(gt[GQpos])
			minD=float(minD)
			maxD=float(maxDPDict[sample])
			if DP < minD:
				gtfilter.append('dpqd_lowminDP')
			if DP > maxD:
				gtfilter.append('dpqd_highmaxDP')
			if GQ < 20.0:
				gtfilter.append('dpqd_lowGQ20')
			if GQ>=20.0 and minD<=DP<=maxD:
				REF=float(gt[ADpos].split(',')[0])
				AB=float(REF/DP)
				if gt[0]=='0/0' and AB < 0.9:
					gtfilter.append('ab_homrefAB09')
				if gt[0]=='0/1' and AB < 0.2:
					gtfilter.append('ab_hetAB02')
				if gt[0]=='0/1' and AB > 0.8:
					gtfilter.append('ab_hetAB08')
				if gt[0]=='1/1' and AB > 0.1:
					gtfilter.append('ab_homaltAB01')
	# get all the filters
	if gtfilter==[]:
		gtfilter.append('PASS')
	outgtfilter=str(';'.join(gtfilter))
	return outgtfilter

# confirm the two gt are right
# GTbefore = before GTfilter() function
# GTafter = after GTfilter() function
# filtertype = GTfiltertype() output
def check_gtfilter(GTbefore, GTafter, filtertype):
	if filtertype=='PASS':
		if GTbefore[:3]==GTafter[:3]:
			return 0
		else:
			return 1
	else:
		if GTafter[:3]=='./.':
			return 0
		else:
			return 1

# get sample names from vcf file
def get_samples(vcf_file):
	samples=[]
	with gzip.open(vcf_file, 'rt') as inVCF:
		for line in inVCF:
			if line.startswith('##'):
				pass
			else:
				for ii in line.split()[9:]: samples.append(ii)
				break
	return samples

# check the DP file has the same samples as the samples in the VCF
def check_samples(samples, maxDPDict):
	dictsamples=maxDPDict.keys()
	dictsamples.sort()
	if samples==dictsamples:
		return 0
	else:
		return 1

# Read csv based filter files
def readDPfile(DP_file):
	with open(DP_file) as f:
		ff=filter(None,csv.reader(f))
		dd={str(rows[0]):float(rows[1]) for rows in ff}
	return dd

# clear variables start of each line
def clearvars():
	varlist=['line', 'INFO', 'infodict', 'ADpos', 'DPpos', 'GQpos', 'formatfields', 'REF', 'ALT', 'outfilter', 'missing', 'excesshet', 'GT_list', 'GT_filter','GT_before', 'GT', 'gtfilter', 'i', 'x']
	for var in varlist:
		if var in globals():
			del globals()[var]
	return None

# split INFO field with more tolerance to INFO fields with Type=Flag
# @line7: should be called from line[7], a str in the INFO field. CANNOT BE A str with '.'
def split_info(line7):
	infodict={}
	if line7!='.':
		INFO=line7.split(';')
		for x in INFO:
			xx=x.split('=')
			if len(xx)==1:
				xx=[str(xx[0]), True]
			else:
				xx=[str(xx[0]), str(xx[1])]
			infodict[xx[0]]=xx[1]
	return infodict

# combine the infodict modified during the genotype filtering
# @infodict: should be derived from split_info(line[7]) function
def combine_info(infodict):
	infolist=[]
	for key,val in sorted(infodict.items()):
		if type(val) is bool:
			infolist.append(str(key))
		else:
			infolist.append('{0}={1}'.format(key, val))
	return infolist

###########################################################
# processing preparations

# Note that these depths includes possibly bad reads
# Set minimum genotype depth (1/3x mean coverage lower than 8)
minD=8.0

# Set maximum genotype depth for each sample (2.5x mean coverage)
maxDPDict=readDPfile(maxD_file) # output a dictionary with float values

# Get list of samples from given vcf files
samples=get_samples(vcf_file)

# check if the maxDPDict has the same values as the samples
if check_samples(samples, maxDPDict)!=0:
	errmsg='inVCF sample names not matching maxDPDict'
	sys.exit(errmsg)

# load data in text mode
inVCF=gzip.open(vcf_file, 'rt')

# Output file (for output filter statistics)
outfile=open(out_file, 'w')
# Write output header
outheader=['CHROM', 'POS', 'FILTER'] + samples
outfile.write('\t'.join(str(x) for x in (outheader)) + '\n')

###########################################################
# main

# Add new header lines for filters being added - for GATK compatibility
for line0 in inVCF:
	if line0.startswith('#'):
		if line0.startswith('##FORMAT'):
			sys.stdout.write('##FILTER=<ID=FAIL_refN,Description="Reference sequence N">\n')
			sys.stdout.write('##FILTER=<ID=FAIL_badRef,Description="Bad reference alleles Not A or T or C or G">\n')
			sys.stdout.write('##FILTER=<ID=FAIL_badAlt,Description="Bad alternate alleles Not A or T or C or G or dot">\n')
			sys.stdout.write('##FILTER=<ID=FAIL_noQUAL,Description="No QUAL in site">\n')
			sys.stdout.write('##FILTER=<ID=FAIL_noINFO,Description="No INFO in site">\n')
			sys.stdout.write('##FILTER=<ID=FAIL_mutType,Description="Bad VariantType">\n')
			sys.stdout.write('##FILTER=<ID=FAIL_noADi,Description="No AD in genotype">\n')
			sys.stdout.write('##FILTER=<ID=FAIL_noDPi,Description="No DP in genotype">\n')
			sys.stdout.write('##FILTER=<ID=FAIL_noGQi,Description="No GQ in genotype">\n')
			sys.stdout.write('##FILTER=<ID=FAIL_noGT,Description="No genotype after filter">\n')
			sys.stdout.write('##FILTER=<ID=WARN_excessHet,Description="More than 75 percent heterozygous genotype">\n')
			sys.stdout.write('##FILTER=<ID=WARN_missing,Description="More than 20 percent missing genotype">\n')
			sys.stdout.write(line0)
			break
		else: sys.stdout.write(line0)

# Main: perform filtering line by line (no restart to seek(0))
for line0 in inVCF:
	if line0.startswith('#'):
		sys.stdout.write(line0); continue

### For all other lines:
# clear variable buffer first
	clearvars()
	line=line0.strip().split('\t')

### Site filtering:
# Keep any filters that have already been applied (like CpGRep)
	sitefilter=[]
	if line[6] not in ('.', 'PASS'): sitefilter.append(line[6])

### Reference must not be N
	if line[3]=='N':
		sitefilter.append('FAIL_refN')
		outfilter=[line[0], line[1], str(';'.join(sitefilter))] + ['site_fail']*len(samples)
		outfile.write('\t'.join(str(x) for x in (outfilter)) + '\n')
		sys.stdout.write('%s\t%s\t%s\n' % ('\t'.join(line[0:6]), ';'.join(sitefilter), '\t'.join(line[7:]))); continue

### Check REF allele
	if line[3] not in ('A','C','G','T'):
		sitefilter.append('FAIL_badRef')

### Check ALT allele
	if line[4] not in ('A','C','G','T','.'):
		sitefilter.append('FAIL_badAlt')

### Must have a valid QUAL
	if line[5]=='.':
		sitefilter.append('FAIL_noQUAL')

### Access INFO field annotations
	infodict={}
	if line[7]=='.':
		sitefilter.append('FAIL_noINFO')
	else:
		infodict=split_info(line[7])

### Only accept sites that are monomorphic or simple SNPs
### Only check the variables with info sections
	if infodict!={}:
		if 'VariantType' not in infodict or infodict['VariantType'] not in ('NO_VARIATION', 'SNP'):
			sitefilter.append('FAIL_mutType')

### Get the position of AD value in genotype fields
	if 'AD' in line[8]:
		ADpos=line[8].split(':').index('AD')
	else:
		sitefilter.append('FAIL_noADi')

### Get the position of DP value in genotype fields
	if 'DP' in line[8]:
		DPpos=line[8].split(':').index('DP')
	else:
		sitefilter.append('FAIL_noDPi')

### Get the position of GQ/RGQ value in genotype fields use the for loop to search for either
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=RGQ,Number=1,Type=Integer,Description="Unconditional reference genotype confidence, encoded as a phred quality -10*log10 p(genotype call is wrong)">
	if 'GQ' in line[8]:
		formatfields=line[8].split(':')
		GQpos=[formatfields.index(x) for x in formatfields if 'GQ' in x][0]
	else:
		sitefilter.append('FAIL_noGQi')

### If any filters failed, write out line and continue
	if sitefilter!=[]:
		outfilter=[line[0], line[1], str(';'.join(sitefilter))] + ['site_fail']*len(samples)
		outfile.write('\t'.join(str(x) for x in (outfilter)) + '\n')
		sys.stdout.write('%s\t%s\t%s\n' % ('\t'.join(line[0:6]), ';'.join(sitefilter), '\t'.join(line[7:])))
		continue

### Genotype filtering:
	missing,excesshet=0,0
	GT_list=[]
	GT_filter=[] # genotype filter annotation
	for i in range(0,len(samples)):
		GT_before=line[i+9]
		GT=GTfilter(samples[i],GT_before,minD,maxDPDict,ADpos,DPpos,GQpos)
		gtfilter=GTfiltertype(samples[i],GT_before,minD,maxDPDict,ADpos,DPpos,GQpos)
		# check if the gt filter annotation is correct
		if check_gtfilter(GT_before, GT, gtfilter)!=0:
			inVCF.close()
			outfile.close()
			errmsg=str(';'.join([line[0], line[1], samples[i], 'GT filters not matching']))
			sys.exit(errmsg)
		if GT[:3]=='./.':
			missing+=1
		if GT[:3]=='0/1':
			excesshet+=1
		GT_list.append(GT)
		GT_filter.append(gtfilter)

### Recalculate INFO fields (count REF and ALT genotypes)
	REF=2*[x[:3] for x in GT_list].count('0/0') + [x[:3] for x in GT_list].count('0/1')
	ALT=2*[x[:3] for x in GT_list].count('1/1') + [x[:3] for x in GT_list].count('0/1')

	if REF+ALT==0:
		sitefilter.append('FAIL_noGT')
		outfilter=[line[0], line[1], str(';'.join(sitefilter))] + GT_filter
		outfile.write('\t'.join(str(x) for x in (outfilter)) + '\n')
		sys.stdout.write('%s\t%s\t%s\t%s\n' % ('\t'.join(line[0:6]), ';'.join(sitefilter), '\t'.join(line[7:9]), '\t'.join(GT_list)))
		continue

### Apply warnings for sites with more than 75% heterozygous ('0/1') genotypes
# Change depending on experimental design!
	if excesshet>0.75*len(samples):
		sitefilter.append('WARN_excessHet')

### Apply warnings for sites with more than 20% missing ('./.') genotypes
	if missing>0.2*len(samples):
		sitefilter.append('WARN_missing')

### Recalculate INFO fields (count REF and ALT genotypes)
	# AC: allele count in genotypes, for each ALT allele, in the same order as listed
	# AN: total number of alleles in called genotypes
	infodict['AC']=ALT
	infodict['AN']=REF+ALT
	infodict['AF']=round(float(ALT)/(float(REF)+float(ALT)), 4)

	if sitefilter==[]:
		sitefilter.append('PASS')

### Write out new line
	outfilter=[line[0], line[1], str(';'.join(sitefilter))] + GT_filter
	outfile.write('\t'.join(str(x) for x in (outfilter)) + '\n')
	sys.stdout.write('%s\t%s\t%s\t%s\t%s\n' % ('\t'.join(line[0:6]), ';'.join(sitefilter), ';'.join(combine_info(infodict)), line[8], '\t'.join(GT_list)))

###########################################################
# cleanup
outfile.close()
inVCF.close()
exit()

