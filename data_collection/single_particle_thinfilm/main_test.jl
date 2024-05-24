using JLD2
include("../../src/operators_thinfilm.jl")
using LinearAlgebra


V0space,am,L1,trace_condition,chern,uniform,gapup,gapdown,direct_gapup,direct_gapdown,bandwidth=main()


jldsave(joinpath(@__DIR__, "data_output/SbTe_L1$(L1)_am$(am).jld2"),V0space=V0space,am=am,L1=L1,trace_condition=trace_condition,chern=chern,uniform=uniform,gapup=gapup,gapdown=gapdown,direct_gapup=direct_gapup,direct_gapdown=direct_gapdown,bandwidth=bandwidth)
