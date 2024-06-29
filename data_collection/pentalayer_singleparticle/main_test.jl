using JLD2
include("../../src/operators_pentalayerperiodicpo.jl")
using LinearAlgebra

args=parse.(Float64,ARGS)

θ=args[1]/180*π
Nq=30
Vperiod=args[2]
uD=args[3]
period=args[4]
Nband=7
phaseangle=args[5]/180*π

eigenvector,eigenvalue,wave,allowedq,allowedq_dic,T1,T2,BW,direct_gapup,direct_gapdown,indirect_gapup, indirect_gapdown=single_particle_periodicpo(θ,Nq,uD,Nband,Vperiod,phaseangle,period)
trace_condition,tra,Flink,chern,uniform=get_chern(Nq,wave,allowedq,allowedq_dic,T1,T2,eigenvector)
jldsave(joinpath(@__DIR__, "data_output/$(args[1])theta$(args[2])Vperiod$(args[3])uD$(args[4])period.jld2"),chern=chern,Flink=Flink,TC=trace_condition,TC_kresolved=tra,uniform=uniform,BW=BW,direct_gapup=direct_gapup,direct_gapdown=direct_gapdown,indirect_gapup=indirect_gapup,indirect_gapdown=indirect_gapdown)
