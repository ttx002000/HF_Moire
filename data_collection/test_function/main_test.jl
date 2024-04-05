using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames
using LinearAlgebra


include("../../src/operators_TDHF.jl")


args=parse.(Float64,ARGS)
flux=args[1]*π
V0=args[2]
ϕ=args[3]/180*π
Nq=Int(args[4]);
scale=args[5];
constq=args[6]/Nq^2
trytimes=Int(args[7])
bigQindex=Int(args[8])


Amatrix,Avec=test_function()
println(Amatrix)
println(Avec)
jldsave(joinpath(@__DIR__, "data_output/spectrum$(args[4])Nq$(args[1])flux$(args[2])V0$(args[3])phi$(args[5])scale$(args[6])constq$(args[7])try$(args[8])bigQ.jld2"),Amatrix=Amatrix)


