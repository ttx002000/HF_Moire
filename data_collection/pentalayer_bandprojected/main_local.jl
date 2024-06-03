using JLD2
include("../../src/operators_pentalayer.jl")
using LinearAlgebra

args=[5.0,0.77,6,50.0,7.0]
ϵr=args[1]
θ=args[2]/180*π
Nq=Int(args[3]);
uD=args[4]
Nband=Int(args[5])




eigenvector,eigenvalue,wave,wave_diff,allowedq,allowedq_dic,T1,T2,form_factors,constq,ϵ=single_particle(ϵr,θ,Nq,uD,Nband)
println("Finished single particle")
initial_projector,band_Ham=get_initial_projector(Nq,Nband,eigenvalue)

DIIS_input_projector,energy,HF_eigenvalue,HF_eigenvector,bound=iteration(Nq,Nband,initial_projector,form_factors,constq,wave_diff,allowedq,allowedq_dic,T1,T2,band_Ham)
trace_condition,tra,Flink,chern,uniform=get_chern(Nq,wave,allowedq,allowedq_dic,T1,T2,eigenvector,HF_eigenvector)

jldsave(joinpath(@__DIR__, "data_output/tryDIIS.jld2"),chern=chern,Flink=Flink,energy=energy,TC=trace_condition,TC_kresolved=tra,HF_eigenvalue=HF_eigenvalue,bound=bound,uniform=uniform,args=args)
