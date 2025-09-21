using JLD2
include("../../src/operators_RMG_bandprojected_interpolation.jl")
using LinearAlgebra

args=parse.(Float64,ARGS)
#args=[5.0,0.77,6,50.0,7.0,1.0,0.0,3.01,5,1,1]
ϵr=args[1]
θ=args[2]/180*π
Nq=Int(args[3]);
uD=args[4]
Nband=Int(args[5])
λ=args[6]
contact_strength=args[7]
gcutoff=args[8]
NL=Int(args[9])
coupling_ratio=args[10]
trytime=Int(args[11])
filepos=Int(args[12])



eigenvector,eigenvalue,wave,wave_diff,allowedq,allowedq_dic,T1,T2,form_factors,Area,ϵ=single_particle(λ,θ,NL,Nq,uD,Nband,gcutoff,coupling_ratio)
push!(args,ϵ)
println("Finished single particle")
initial_projector,band_Ham=get_initial_projector(Nq,Nband,eigenvalue)

DIIS_input_projector,energy,HF_eigenvalue,HF_eigenvector,bound=iteration(Nq,Nband,initial_projector,form_factors,Area,ϵr,contact_strength,wave_diff,allowedq,allowedq_dic,T1,T2,band_Ham)
trace_condition,tra,Flink,chern,uniform=get_chern(NL,Nq,wave,allowedq,allowedq_dic,T1,T2,eigenvector,HF_eigenvector)



scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "triangle_R5G_bandproject_interpolation_withhBN/data_output$(Int(args[12]))/$(args[1])er$(args[2])theta$(args[3])Nq$(args[4])uD$(args[5])Nband$(args[6])lambda$(args[7])contact$(args[8])cutoff$(args[9])NL$(args[10])couplingratio$(args[11])trytime.jld2")


jldsave(savepath,
     chern=chern,Flink=Flink,energy=energy,TC=trace_condition,TC_kresolved=tra,
     HF_eigenvalue=HF_eigenvalue,bound=bound,uniform=uniform,args=args,eigenvector=eigenvector,eigenvalue=eigenvalue,
     final_projector=DIIS_input_projector[1],HF_eigenvector=HF_eigenvector,allowedq=allowedq,wave_diff=wave_diff,T1=T1,T2=T2)
