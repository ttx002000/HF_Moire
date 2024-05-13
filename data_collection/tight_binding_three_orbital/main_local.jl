using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames
using LinearAlgebra

include(joinpath(@__DIR__,"../../src/operators_threeorbital.jl"))

observable,parameters,energymatrix,eoutmatrix=excecute_loop()
jldsave(joinpath(@__DIR__, "data_output/try1.jld2"),observable=observable,parameters=parameters,energymatrix=energymatrix,eoutmatrix=eoutmatrix)



