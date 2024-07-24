using JLD2
include("../../src/operators_thinfilm.jl")
using LinearAlgebra

args=parse.(Int64,ARGS)
Nb=args[1]
ϕ=args[2]/180*π
V0space,am,L1,trace_condition,chern,uniform,gapup,gapdown,direct_gapup,direct_gapdown,bandwidth=main(Nb,ϕ)
#V0space,am,L1,gapup,gapdown,direct_gapup,direct_gapdown,bandwidth=main_SbTe_zdependence(Nb,ϕ)


jldsave(joinpath(@__DIR__, "data_output/CdAs_L1$(L1)_am$(am)_valence_zdep_$(args[1])band$(args[2])phi.jld2"),V0space=V0space,am=am,L1=L1,trace_condition=trace_condition,chern=chern,uniform=uniform,gapup=gapup,gapdown=gapdown,direct_gapup=direct_gapup,direct_gapdown=direct_gapdown,bandwidth=bandwidth)
#jldsave(joinpath(@__DIR__, "data_output/SbTe_L1$(L1)_am$(am)_conduction_zdep_$(args[1])band$(args[2])phi.jld2"),V0space=V0space,am=am,L1=L1,gapup=gapup,gapdown=gapdown,direct_gapup=direct_gapup,direct_gapdown=direct_gapdown,bandwidth=bandwidth)
