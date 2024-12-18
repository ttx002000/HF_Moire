
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2
using Plots
using CSV,DataFrames

include(joinpath(@__DIR__,"../../src/operators_DOS_pentalayer_withmoire.jl"))

args=parse.(Float64,ARGS)
#args=[-50.0,10*10^5,1.5,50,10^4,0.2/4,4/4,200]
uD=args[1]
numsample=Int(args[2])
θ=args[3]/180*π
DOS_n_binnum=Int(args[4])
DOS_E_binnum=Int(args[5])
Density_start=args[6]
Density_end=args[7]
Ecutoff=args[8]

valuesset,ns,gm,conduction_bandmin,conduction_bandmax,valence_bandmin,valence_bandmax=sample_value(uD, numsample,θ,Ecutoff)
CNP_point=(valence_bandmax+conduction_bandmin)/2
nE_bin_centers,nE_bin_means,bin_centers,nE,Nstates=process_data(valuesset,ns,numsample,DOS_n_binnum,DOS_E_binnum,Density_start,Density_end,gm,CNP_point)

scratch_dir = ENV["SCRATCH"]

savepath=joinpath(scratch_dir, "pentalayer_DOS/data_output/with_moire/withmoire_$(args[1])uD$(args[2])sample$(args[3])angle$(args[4])Dosnbin$(args[5])DosEbin$(args[6])denstart$(args[7])denend$(args[8])Ecut.jld2")
jldsave(savepath,nE_bin_centers=nE_bin_centers,
   ns=ns,nE_bin_means=nE_bin_means,
   bin_centers=bin_centers,nE=nE,
   Nstates=Nstates,conduction_bandmin=conduction_bandmin,conduction_bandmax=conduction_bandmax,
   valence_bandmin=valence_bandmin,valence_bandmax=valence_bandmax,
   final_num_sample=length(valuesset))
   

  