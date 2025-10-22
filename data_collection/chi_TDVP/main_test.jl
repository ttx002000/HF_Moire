using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames
using LinearAlgebra

include("../../src/operators_chi_TDVP.jl")


args=parse.(Float64,ARGS)

Nq=Int(args[1])
powerindex=Int(args[2])
bigQindex=Int(args[3])
file_pos=Int(args[4])

 MHmatrix, Bmatrix,MMmatrix=big_func(powerindex,bigQindex,file_pos)


scratch_dir = ENV["SCRATCH"]
output_path=joinpath(scratch_dir, "chi_TDVP/data_output$(Int(args[4]))/$(args[1])Nq$(args[2])powerindex$(args[3])bigQ")


 jldsave(output_path,
         Amatrix=MHmatrix,Bmatrix=Bmatrix,Gmatrix=MMmatrix)


