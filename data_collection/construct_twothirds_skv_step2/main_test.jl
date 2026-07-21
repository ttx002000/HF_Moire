using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2


include("../../src/construct_twothirds_skv_step2.jl")

args=parse.(Float64,ARGS)

main_func(args)

