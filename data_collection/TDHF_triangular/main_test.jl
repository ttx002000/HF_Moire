using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames


include("../../src/operators_TDHF.jl")

args=parse.(Float64,ARGS)
flux=args[1]*π
V0=args[2]
ϕ=args[3]/180*π
Nq=Int(args[4]);
scale=args[5];
constq=args[6]/Nq^2
trytimes=Int(args[7])


overlapmatrix, wave, single_MoirePo, single_Ham,allowedq, T1, T2=single_part(flux,V0,ϕ,scale,Nq)

 


