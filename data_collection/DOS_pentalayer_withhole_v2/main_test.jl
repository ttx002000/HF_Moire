
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2
using Plots
using CSV,DataFrames

include(joinpath(@__DIR__,"../../src/operators_DOS_pentalayer_withhole_v2.jl"))

args=parse.(Float64,ARGS)
#args=[10.0,10*10^5,-0.02,0.02,1,120,0.01,1]
uD=args[1]
numsample=Int(args[2])
nstart=args[3]
nend=args[4]
perturb=Int(args[5])
DOS_n_binnum=Int(args[6])
Einterval=args[7]
file_pos=args[8]

perturb_Ham=get_perturb_Ham(perturb)
E_lower,E_upper, kradius=find_E_cut(nstart,nend,uD,perturb_Ham)
valuesset, N_valence, DOS_E_binnum=sample_value(uD, numsample,E_lower,E_upper,perturb_Ham,kradius,Einterval)
nE_bin_centers,nE_bin_means,bin_centers,nE,Nstates=process_data(valuesset,numsample,kradius,DOS_n_binnum,DOS_E_binnum,nstart,nend)

scratch_dir = ENV["SCRATCH"]

savepath=joinpath(scratch_dir, "pentalayer_DOS_v2/data_output$(Int(args[8]))/$(args[1])uD$(args[2])sample$(args[3])nstart$(args[4])nend$(args[5])perturb$(args[6])DosEbin$(args[7])Eint.jld2")
jldsave(savepath,nE_bin_centers=nE_bin_centers,
   nE_bin_means=nE_bin_means,
   bin_centers=bin_centers,nE=nE,
   Nstates=Nstates,E_lower=E_lower,E_upper=E_upper,kradius=kradius,
   final_num_sample=length(valuesset),DOS_E_binnum=DOS_E_binnum,N_valence=N_valence)
   

  