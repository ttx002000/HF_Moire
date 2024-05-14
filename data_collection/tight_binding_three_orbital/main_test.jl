using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames
using LinearAlgebra

args=parse.(Int64,ARGS)
trytime=args[1]
include(joinpath(@__DIR__,"../../src/operators_threeorbital.jl"))

observable,parameters,energymatrix,eoutmatrix,single_eigenvalue,HF_eigenvalue=excecute_loop()

jldsave(joinpath(@__DIR__, "data_output/try$(trytime)_tprime=minus0.02.jld2"),observable=observable,parameters=parameters,energymatrix=energymatrix,eoutmatrix=eoutmatrix,single_eigenvalue=single_eigenvalue,HF_eigenvalue=HF_eigenvalue)




