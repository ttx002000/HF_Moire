
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2
using Plots
using CSV,DataFrames

include(joinpath(@__DIR__,"../../src/operators_DOS_pentalayer.jl"))

args=parse.(Float64,ARGS)
#args=[30.0,10*10^5,1.5,1.0,50,10^4,0.3/4,5/4]
uD=args[1]
numsample=Int(args[2])
θ=args[3]/180*π
rad=args[4]
DOS_n_binnum=Int(args[5])
DOS_E_binnum=Int(args[6])
Density_start=args[7]
Density_end=args[8]
Ecutoff=args[9]

valuesset,ns,gm,bandmin,bandmax=sample_value(uD, numsample,θ,rad,Ecutoff)
nE_bin_centers,nE_bin_means,bin_centers,nE,Nstates=process_data(valuesset,ns,numsample,rad,DOS_n_binnum,DOS_E_binnum,Density_start,Density_end,gm)
savepath=joinpath(@__DIR__, "data_output/$(args[1])uD$(args[2])sample$(args[3])angle$(args[4])radius$(args[5])Dosnbin$(args[6])DosEbin$(args[7])denstart$(args[8])denend$(args[9])Ecut.jld2")
jldsave(savepath,nE_bin_centers=nE_bin_centers,ns=ns,nE_bin_means=nE_bin_means,bin_centers=bin_centers,nE=nE,Nstates=Nstates,bandmin=bandmin,bandmax=bandmax,final_num_sample=length(valuesset))