using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

using LinearAlgebra
using JLD2
using Plots

BLAS.set_num_threads(1)

include("../../src/operators_R5G_contactinteraction_interpolation_withhBN_skyrmionexcitation_v6_time_evolution_Strang.jl")

args = parse.(Float64, ARGS)

run_tdhf_from_args!(args)


