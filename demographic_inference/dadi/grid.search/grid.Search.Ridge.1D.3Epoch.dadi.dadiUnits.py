# -*- coding: utf-8 -*-
'''
Title: Grid search for ENP population, 1D3Epoch (whaling bottleneck). Only do the search within the given ridge regions.
Notes: 1. TF+TB and nuB are fixed according to the best estimation.
       2. Bottleneck function used pnunez version.
Author: annabelbeichman, pnunez
Modified by Meixi Lin (meixilin@ucla.edu)
Date: Wed Aug  4 18:17:13 2021
Example usage:
python grid.Search.Ridge.1D.3Epoch.dadi.dadiUnits.py --pop 'ENP' --sfs SFS/ENP-44.sfs --TBplusTF_Fix 0.134 ... --outdir <output directory>

Notes on the dadi units:
nu: Ratio of contemporary to ancient population size
T: Time in the past at which size change happened (in units of 2*Na)
Nanc is an input parameter, and are set externally from past runs for each population
nu and T are searched through a grid
but should be relative to Nanc
'''

###########################################################
# %% import packages
import math # to get log10
import numpy as np
import sys
import os
import argparse
import dadi
from dadi import Numerics, PhiManip, Integration, Spectrum
# note: module  Demographics1D has the following models file:///Users/annabelbeichman/Documents/UCLA/dadi/dadi/doc/api/dadi.Demographics1D-module.html (not importing Demographics1D, since using the custom model)
# has 2 epoch, and 3 epoch (bottleneck)

from numpy import array # don't comment this out
import datetime

# https://stackoverflow.com/questions/13370570/elegant-grid-search-in-python-numpy
from sklearn.model_selection import ParameterGrid

###########################################################
# %% def functions
#### input parsers
def gridSearchParser(modelName):
    parser = argparse.ArgumentParser(description='Carry out a grid search for a '+ modelName +' model ')
    parser.add_argument("--pop",type=str,required=True,help="population identifier, e.g. ENP or GOC")
    parser.add_argument("--sfs",type=str,required=True,help="path to DATA FOLDED SFS in dadi format from easysfs (mask optional)")
    parser.add_argument("--L",type=int,required=True,help="number of called neutral sites that went into making SFS (monomorphic+polymorphic)")
    parser.add_argument("--TBplusTF_Fix",type=float,required=True,help="input FIXED duration of bottleneck plus whaling time. Provide in dadi units")
    parser.add_argument("--nuB_Fix",type=float,required=True,help="input FIXED pre-whaling population size. Provide in dadi units")
    parser.add_argument("--nuF_Low",type=float,required=True,help="input lower bound on nuF (whaling size) in dadi units -- distance between low and high will be searched evenly along a log scale")
    parser.add_argument("--nuF_High",type=float,required=True,help="input upper bound on nuF (whaling size) in dadi units -- distance points between nu_Low and nu_High will be searched evenly along a log scale")
    parser.add_argument("--TF_Low",type=float,required=True,help="input lower bound on TF (duration of whaling period; TF+TB is fixed at given TBplusTF in dadi units; distance between TF_Low and TF_High will be searched evenly along a log scale")
    parser.add_argument("--TF_High",type=float,required=True,help="input upper bound on TF (duration of whaling period; TF+TB is fixed at given TBplusTF in dadi units; distance between TF_Low and TF_High will be searched evenly along a log scale")
    parser.add_argument("--intercept_Low",type=float,required=True,help="The intercept of the bounding line. In log10 space.")
    parser.add_argument("--intercept_High",type=float,required=True,help="The intercept of the bounding line. In log10 space.")
    parser.add_argument("--slope",type=float,required=True,help="The slope of the bounding line. In log10 space.")
    parser.add_argument("--numGridPoints",type=int,required=True,help="number of grid points per parameter you want. Keep in mind this can drastically affect run time (e.g. 10 --> 100 calculations; 25 -->  625 calculations")
    parser.add_argument("--outdir",type=str,required=True,help="path to output directory")
    return parser

#### the bottleneck model used for inference
# source file: /u/project/rwayne/pnunez/FinWhale/dadi/1D.1Bottleneck.dadi.py (Last modified Feb 16, 2021)
def bottleneck(params, ns, pts):
    nuB,nuF,TB,TF = params
    xx = Numerics.default_grid(pts) # sets up grid
    phi = PhiManip.phi_1D(xx) # sets up initial phi for population
    phi = Integration.one_pop(phi, xx, TB, nuB)  # bottleneck
    phi = Integration.one_pop(phi, xx, TF, nuF) # recovery
    fs = Spectrum.from_phi(phi, ns, (xx,))
    return fs

# function to get the exp SFS and LL based on parameters YOU provide
def provideParams_3EpochWhalingModel(fs,nuB,nuF,TB,TF,func,func_ex,mu,L):
    # get sample size from sfs
    ns= fs.sample_sizes
    # get pts from sample size
    pts_l = [ns[0]+5,ns[0]+15,ns[0]+25]
    # get the expected sfs for this set of parameters:
    # parameters in order: (nuB,nuF,TB,TF)
    model=func_ex([nuB,nuF,TB,TF],ns,pts_l)
    # get the LL of model
    ll_model = dadi.Inference.ll_multinom(model, fs)
    # also get the LL of the data to itself (best possible ll)
    ll_data=dadi.Inference.ll_multinom(fs, fs)
    # HINT: Different from AB's script. Here Nanc is calculated from the SFS. In the follow up plotting R script, I added an option to set it externally (using previously inferred best theta0)
    theta = dadi.Inference.optimal_sfs_scaling(model, fs)
    Nanc=theta / (4*mu*L)
    nuB_scaled_dip=nuB*Nanc
    nuF_scaled_dip=nuF*Nanc
    TB_scaled_gen=TB*2*Nanc
    TF_scaled_gen=TF*2*Nanc

    # fold exp sfs?
    model_folded=model.fold()
    # note use of sfs.tolist() has to have () at the end
    # otherwise you get weird newlines in the mix
    output=[nuB,nuF,TB,TF,theta,
            Nanc,nuB_scaled_dip,nuF_scaled_dip,TB_scaled_gen,TF_scaled_gen,
            ll_model,ll_data,func.func_name,ns[0],model_folded.tolist(),fs.tolist()] # put all the output terms together
    return(output)

###########################################################
# %% def variables
todaysdate=datetime.datetime.today().strftime('%Y%m%d')
modelName="1D.3Epoch"
mu=2.77e-8

# testing variables and justifications MAKE SURE THEY ARE IN DADI UNITS
# SFS (Dec 20, 2020): /u/project/rwayne/snigenda/finwhale/SFS/Neutral/SFS_projection_fsc/ENP-44.sfs
# ENP run: 1D3Epoch, version 1 rundate: 20201223
# TBplusTF_Fix (April 12, 2020, ENP 1D3Epoch, version 1): TB + TF in dadi units for the best run. The top 15 ranking runs fall within (0.127,0.140)
# nuB_Fix (April 12, 2020, ENP 1D3Epoch, version 1): nuB in dadi units for the best run. The top 15 ranking runs fall within (1.430,1.482)
# nuF and TF low and high: upper and lower bound in dadi settings. TF_High has to be smaller than TBplusTF. Realistically, nuF is not larger than nuB_Fix a lot

#args={'pop': 'ENP',
#     'sfs': '/Users/linmeixi/Lab/fin_whale/scripts_analyses/dadi/grid.search/SFS/ENP-44.sfs',
#     'L': 392707916,
#     'TBplusTF_Fix': 0.134235802901809,
#     'nuB_Fix': 1.45113705136459,
#     'nuF_Low': 0.0001,
#     'nuF_High': 2.,
#     'TF_Low': 1e-05,
#     'TF_High': 1.,
#     'intercept_Low': -2.904,
#     'intercept_High': -2.004,
#     'slope':1.041,
#     'numGridPoints': 100,
#     'outdir': 'test'}

###########################################################
# %% main
def main():
    #### parse variables
    parser = gridSearchParser(modelName)
    prog = parser.prog
    args = vars(parser.parse_args())
    # print some logs
    sys.stdout.write('Beginning execution of {0} in directory {1}\n'.format(prog, os.getcwd()))
    sys.stdout.write('Parsed the following arguments:\n{0}\n'.format(
            '\n'.join(['\t{0} = {1}'.format(*tup) for tup in args.items()])))

    pop=args['pop']
    sfs=args['sfs']
    L=args['L']
    TBplusTF_Fix=args['TBplusTF_Fix']
    nuB_Fix=args['nuB_Fix']
    nuF_Low=args['nuF_Low']
    nuF_High=args['nuF_High']
    TF_Low=args['TF_Low']
    TF_High=args['TF_High']
    if TF_High > TBplusTF_Fix:
        TF_High = TBplusTF_Fix
    intercept_Low=args['intercept_Low']
    intercept_High=args['intercept_High']
    slope=args['slope']

    numGridPoints=args['numGridPoints']
    outdir=args['outdir']

    #### input data
    fs=dadi.Spectrum.from_file(sfs) # this is folded if from easy SFS

    # check if it's folded, if not folded, fold it
    if fs.folded==False:
        fs=fs.fold()
    else:
        fs=fs

    #### set up demography function
    func=bottleneck
    param_names= ("nuB","nuF","TB","TF")
    # Make extrapolation function:
    func_ex = dadi.Numerics.make_extrap_log_func(func) # this will give you the exp SFS

    #### set up grids
    # from numpy use linspace where you tell it how many points you want between two numbers
    # note that log10(nu_Low_rescaled) gives you the exponent to put in as the start
    # for example if you wanted to search 10 points between 0.0001 and 0.001
    # you would do np.logspace(-4,-3,10) where -4 and -3 are the exponents
    # so to get those exponents you take the base10 log of 0.0001 to get -4 and put that in: *make sure to use log10!!! (not natural log) ***
    nuFs= np.logspace(math.log10(nuF_Low),math.log10(nuF_High),numGridPoints) # return values evenly spaced along log scale,from base**start to base**stop base =10, get 10 data points
    TFs= np.logspace(math.log10(TF_Low),math.log10(TF_High),numGridPoints)
    # set up your set of parameters
    param_grid = {'nuF' : nuFs, 'TF' : TFs}
    # use the sklearn module to make a list with every pair of parameters in it:
    # similar to R expand
    grid = ParameterGrid(param_grid)

    #%%### set up output file
    outputFile=open(str(outdir)+"/"+"dadi.grid.search."+str(pop)+"."+str(modelName)+".LL.output.txt","w")
    other_output=("theta", "Nanc","nuB_scaled_dip","nuF_scaled_dip","TB_scaled_gen","TF_scaled_gen",
                  "LL_model", "LL_data", "modelFunction", "sampleSize","expectedSFS_fold_Theta1", "observedSFS_folded",
                  "pop", "rundate")
    header='\t'.join(str(x) for x in param_names+other_output) + '\n'
    outputFile.write(header)

    #%%### run the grid search
    # run the function on the grid of all parameter pairs:
    gridcount=0
    for params in grid:
        # test if the parameter fall within the given range
        min_TF=10.**(slope*(math.log10(params['nuF'])) + intercept_Low)
        max_TF=10.**(slope*(math.log10(params['nuF'])) + intercept_High)
        if params['TF'] < max_TF and params['TF'] > min_TF:
            sys.stdout.write('Working  {0}\n'.format(params))
            gridcount+=1
            # calculate TB = TB+TF-TF
            TB_temp=TBplusTF_Fix - params['TF']
            output=provideParams_3EpochWhalingModel(fs=fs,
                                             nuB=nuB_Fix,
                                             nuF=params['nuF'],
                                             TB=TB_temp,
                                             TF=params['TF'],
                                             func=func,
                                             func_ex=func_ex,
                                             mu=mu,
                                             L=L)
            output=output+[pop,todaysdate]
            output='\t'.join(str(x) for x in output)
            outputFile.write(output)
            outputFile.write("\n")
        else:
            # sys.stdout.write('Skipping {0}\n'.format(params))
            pass

    # print the final number of grids
    sys.stdout.write('Final parameter pairs tested = {0}\n'.format(gridcount))

    outputFile.close()


#%%
if __name__ == "__main__":
    sys.exit(main())
