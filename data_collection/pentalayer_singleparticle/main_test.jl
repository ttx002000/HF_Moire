using JLD2
include("../../src/operators_pentalayerperiodicpo.jl")
using LinearAlgebra

args=parse.(Float64,ARGS)

θ=args[1]
Nq=15
Vperiod=args[2]
uD=args[3]
Nband=7


eigenvector,eigenvalue,wave,allowedq,allowedq_dic,T1,T2,BW,direct_gapup,direct_gapdown,indirect_gapup, indirect_gapdown=single_particle_periodicpo(θ,Nq,uD,Nband,Vperiod)
trace_condition,tra,Flink,chern,uniform=get_chern(Nq,wave,allowedq,allowedq_dic,T1,T2,eigenvector)
jldsave(joinpath(@__DIR__, "data_output/$(args[1])theta$(args[2])Vperiod$(args[3])uD.jld2"),chern=chern,Flink=Flink,TC=trace_condition,TC_kresolved=tra,uniform=uniform,BW=BW,direct_gapup=direct_gapup,direct_gapdown=direct_gapdown,indirect_gapup=indirect_gapup,indirect_gapdown=indirect_gapdown)
