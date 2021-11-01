# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 16:46:27 2018
@author: annabelbeichman
Wed Aug  5 15:27:39 2020
@modification: Meixi Lin
@modification: Paulina Nunez-Valencia
"""
import matplotlib
matplotlib.use('Agg') # so graphics show up on hoffman
import sys
import argparse
import dadi
from dadi import Numerics, PhiManip, Integration, Spectrum
from numpy import array # don't comment this out
import datetime
todaysdate=datetime.datetime.today().strftime('%Y%m%d')


modelName="1D.4Epoch"

############### Parse input arguments ########################
parser = argparse.ArgumentParser(description='Infer a '+ modelName +' model from a 1D folded SFS in dadi')
parser.add_argument("--runNum",required=True,help="iteration number (e.g. 1-50)")
parser.add_argument("--pop",required=True,help="population identifier, e.g. 'ENP/GOC'")
parser.add_argument("--mu",required=True,help="supply mutation rate in mutation/bp/gen")
parser.add_argument("--L",required=True,help="number of called neutral/synonymous sites that went into making SFS (monomorphic+polymorphic)")
parser.add_argument("--sfs",required=True,help="path to FOLDED SFS in dadi format from easySFS (mask optional)")
parser.add_argument("--outdir",required=True,help="path to output directory")
# usage:
# python 1D.Bottleneck.dadi.py --runNum $i --pop CA --mu 8.64411385098638e-09 --L 4193488 --sfs [path to sfs] --outdir [path to outdir]
args = parser.parse_args()
runNum=str(args.runNum)
pop=str(args.pop)
mu=float(args.mu)
L=float(args.L)
outdir=str(args.outdir)
sfs=str(args.sfs)
maxiter=100
############### Input data ####################################
fs=dadi.Spectrum.from_file(sfs) # this is folded from easy SFS

if fs.folded==False:
    fs=fs.fold()
else:
    fs=fs
############### Set up General Dadi Parameters ########################
ns = fs.sample_sizes # get sample size from SFS (in haploids)
pts_l = [ns[0]+5,ns[0]+15,ns[0]+25] # this should be slightly larger (+5) than sample size and increase by 10
############### Set up Specific Model -- this will change from script to script ########################
# Not fixing the bottleneck duration for now

def double_bottleneck(params, ns, pts):
    nuB,nuR,nuC,TB,TR,TC = params
    xx = Numerics.default_grid(pts) # sets up grid
    phi = PhiManip.phi_1D(xx) # sets up initial phi for population
    phi = Integration.one_pop(phi, xx, TB, nuB)  
    phi = Integration.one_pop(phi, xx, TR, nuR) 
    phi = Integration.one_pop(phi, xx, TC, nuC) 
    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs
param_names=("nuB","nuR","nuC","TB","TF","TC")

# for each lower bound, should be larger than 0
# for the upper bound, depending on specific conditions
upper_bound = [5, 5, 5, 3, 3, 3]
lower_bound = [1e-5, 1e-5, 1e-7, 1e-10, 1e-8, 1e-10]
p0 = [0.01,0.001,0.06,0.09,0.1,0.1] # initial parameters


func=double_bottleneck # set the function

############### Carry out optimization (same for any model) ########################
# Make extrapolation function:
func_ex = dadi.Numerics.make_extrap_log_func(func)
# perturb parameters
p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,lower_bound=lower_bound)
# optimize:
print('Beginning optimization ************************************************')
popt = dadi.Inference.optimize_log(p0, fs, func_ex, pts_l,
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(p0), maxiter=maxiter)
print('Finshed optimization **************************************************')

# Calculate the best-fit model AFS.
model = func_ex(popt, ns, pts_l)
# Likelihood of the data given the model AFS.
ll_model = dadi.Inference.ll_multinom(model, fs)

# also get the LL of the data to itself (best possible ll)
ll_data=dadi.Inference.ll_multinom(fs, fs)

# calculate best fit theta
theta = dadi.Inference.optimal_sfs_scaling(model, fs)

###### model specific scaling of parameters (will depend on mu and L that you supply) #######

Nanc=theta / (4*mu*L) # based on the definition of theta and units of mutation rate
nuB_scaled_dip=popt[0]*Nanc
nuR_scaled_dip=popt[1]*Nanc
nuC_scaled_dip=popt[2]*Nanc
TB_scaled_gen=popt[3]*2*Nanc # t was scaled by 2N
TF_scaled_gen=popt[4]*2*Nanc
TC_scaled_gen=popt[5]*2*Nanc
scaled_param_names=("Nanc_FromTheta_scaled_dip","nuB_scaled_dip","nuR_scaled_dip","nuC_scaled_dip","TB_scaled_gen","TF_scaled_gen","TC_scaled_gen")
scaled_popt=(Nanc,nuB_scaled_dip,nuR_scaled_dip,nuC_scaled_dip,TB_scaled_gen,TF_scaled_gen,TC_scaled_gen)
############### Write out output (same for any model) ########################
print('Writing out parameters **************************************************')

outputFile=open(str(outdir)+"/"+str(pop)+".dadi.inference."+str(modelName)+".runNum."+str(runNum)+".output","w")
# get all param names:
param_names_str='\t'.join(str(x) for x in param_names)
scaled_param_names_str='\t'.join(str(x) for x in scaled_param_names)
header=param_names_str+"\t"+scaled_param_names_str+"\ttheta\tLL\tLL_data\tmodelFunction\tmu\tL\tmaxiter\trunNumber\trundate\tinitialParameters\tupper_bound\tlower_bound" # add additional parameters theta, log-likelihood, model name, run number and rundate
popt_str='\t'.join(str(x) for x in popt) # get opt'd parameters as a tab-delim string
scaled_popt_str='\t'.join(str(x) for x in scaled_popt)
# joint together all the output fields, tab-separated:
output=[popt_str,scaled_popt_str,theta,ll_model,ll_data,modelName,mu,L,maxiter,runNum,todaysdate,p0,upper_bound,lower_bound] # put all the output terms together
output='\t'.join(str(x) for x in output) # write out all the output fields
# this should result in a 2 row table that could be input into R / concatenated with other runs
outputFile.write(('{0}\n{1}\n').format(header,output))
outputFile.close()

############### Output SFS ########################
print('Writing out SFS **************************************************')

outputSFS=str(outdir)+"/"+str(pop)+".dadi.inference."+str(modelName)+".runNum."+str(runNum)+"."+str(todaysdate)+".expSFS"

model.to_file(outputSFS)


############### Output plot ########################
print('Making plots **************************************************')


import matplotlib.pyplot as plt
fig=plt.figure(1)

outputFigure=str(str(outdir)+"/"+str(pop)+".dadi.inference."+str(modelName)+".runNum."+str(runNum)+"."+str(todaysdate)+".figure.png")
dadi.Plotting.plot_1d_comp_multinom(model, fs)

plt.savefig(outputFigure)



###### exit #######
sys.exit()
