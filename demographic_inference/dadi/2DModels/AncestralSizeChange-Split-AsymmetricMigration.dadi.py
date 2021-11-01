# -*- coding: utf-8 -*-
"""
Created on Tue Oct 16 16:46:27 2018
@author: annabelbeichman/paulinanunezvalencia
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


modelName="2D.4Epochs.ASiMig.Ancestral"

############### Parse input arguments ########################
parser = argparse.ArgumentParser(description='Infer a '+ modelName +' model from a 2D folded SFS in dadi')
parser.add_argument("--runNum",required=True,help="iteration number (e.g. 1-50)")
parser.add_argument("--pop",required=True,help="population identifier, e.g. 'CA'")
parser.add_argument("--mu",required=True,help="supply mutation rate in mutation/bp/gen")
parser.add_argument("--L",required=True,help="number of called neutral sites that went into making SFS (monomorphic+polymorphic)")
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

# check if it's folded, if not folded, fold it
if fs.folded==False:
    fs=fs.fold()
else:
    fs=fs
############### Set up General Dadi Parameters ########################
ns = fs.sample_sizes # get sample size from SFS (in haploids)
pts_l = [ns[0]+5,ns[0]+15,ns[0]+25] # need 6 points because two populations
############### Set up Specific Model -- this will change from script to script ########################

def asym_mig_size(params, ns, pts):
    
     
    Ta, nua, Td, nu1, nu2, m12, m21 = params
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)

    phi = Integration.one_pop(phi, xx, T=Ta, nu=nua) #Ancetral size change
   
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, Td, nu1=nu1, nu2=nu2,m12=m12, m21=m21) #Divergence time

        
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs    


param_names=("Ta","nua","Td","nu1","nu2","m12","m21")

lower_bound = [1e-4, 1e-9, 1e-5, 1e-5, 1e-5, 1e-9, 1e-5]
upper_bound = [10, 20, 10, 10,5,10,700] # 20 as upper bound on mig rec by dadi
p0 = [0.1133,1.542,0.0188,1.063,0.00696,3.020,111.69] # initial parameters

func=asym_mig_size # set the function

############### Carry out optimization (same for any model) ########################
# Make extrapolation function:
func_ex = dadi.Numerics.make_extrap_log_func(func)
# perturb parameters
p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound, lower_bound=lower_bound)
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
# calculate best fit theta
theta = dadi.Inference.optimal_sfs_scaling(model, fs)

###### model specific scaling of parameters (will depend on mu and L that you supply) #######

Nanc=theta / (4*mu*L)
nua_scaled_dip=popt[1]*Nanc
nu1_scaled_dip=popt[3]*Nanc
nu2_scaled_dip=popt[4]*Nanc

Ta_scaled_gen=popt[0]*2*Nanc
Td_scaled_gen=popt[2]*2*Nanc

m12_fraction=popt[5]/(2*Nanc)
m21_fraction=popt[6]/(2*Nanc)


scaled_param_names=("Nanc","Ta","nua","Td","nu1", "nu2", "migrationFraction12", "migrationFraction21")

scaled_popt=(Nanc,Ta_scaled_gen,nua_scaled_dip,Td_scaled_gen,nu1_scaled_dip,nu2_scaled_dip,m12_fraction,m21_fraction)

############### Write out output (same for any model) ########################
print('Writing out parameters **************************************************')

outputFile=open(str(outdir)+"/"+str(pop)+".dadi.inference."+str(modelName)+".runNum."+str(runNum)+".output","w")
# get all param names:
param_names_str='\t'.join(str(x) for x in param_names)
scaled_param_names_str='\t'.join(str(x) for x in scaled_param_names)
header=param_names_str+"\t"+scaled_param_names_str+"\ttheta\tLL\tmu\tL\tmaxiter\trunNumber\trundate\tinitialParameters\tupper_bound\tlower_bound" # add additional parameters theta, log-likelihood, model name, run number and rundate
popt_str='\t'.join(str(x) for x in popt) # get opt'd parameters as a tab-delim string
scaled_popt_str='\t'.join(str(x) for x in scaled_popt)
# joint together all the output fields, tab-separated:
output=[popt_str,scaled_popt_str,theta,ll_model,mu,L,maxiter,runNum,todaysdate,p0,upper_bound,lower_bound] # put all the output terms together
output='\t'.join(str(x) for x in output) # write out all the output fields
# this should result in a 2 row table that could be input into R / concatenated with other runs
outputFile.write(('{0}\n{1}\n').format(header,output))
outputFile.close()

############### Output SFS ########################
print('Writing out SFS **************************************************')

outputSFS=str(outdir)+"/"+str(pop)+".dadi.inference."+str(modelName)+".runNum."+str(runNum)+".expSFS"

model.to_file(outputSFS)


############### Output plot ########################
print('Making plots **************************************************')


import matplotlib.pyplot as plt
fig=plt.figure(1)

outputFigure=str(str(outdir)+"/"+str(pop)+".dadi.inference."+str(modelName)+".runNum."+str(runNum)+".figure.png")
dadi.Plotting.plot_2d_comp_multinom(model, fs)

plt.savefig(outputFigure)



###### exit #######
sys.exit()
