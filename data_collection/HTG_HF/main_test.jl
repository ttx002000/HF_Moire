using JLD2
include("../../src/operators_HTG_HF.jl")
using LinearAlgebra

args=parse.(Float64,ARGS)


wAA=75.0
wAB=110.0
#vF=579.2265 #unit meV*nm
vF=580.25

θ=1.8/180*π
ϵr=8.0
Nband=Int(2)
geonum=Int(2)
filling=7

eigenvector,eigenvalue,wave,wave_diff,wave_dic,allowedq,allowedq_dic,T1,T2,constq,Minv,g1mT,g2mT=single_particle(geonum,θ,wAA,wAB,vF,ϵr,Nband)
Npa=length(allowedq)*filling
formfactors=get_formfactors(allowedq,wave,wave_diff,wave_dic,Minv,Nband,eigenvector)
initial_projector, bg_projector, single_Ham=get_initial_proj(allowedq,eigenvalue,Nband)


HF_eigenvalue,HF_eigenvector,energy, DIIS_input_projector,bound=iteration(formfactors,initial_projector,bg_projector,constq,Nband,wave_diff,allowedq,allowedq_dic,T1,T2,single_Ham,Npa)
println(energy)
jldsave(joinpath(@__DIR__,"data.jld2"),HF_eigenvalue=HF_eigenvalue,eigenvalue=eigenvalue,energy=energy)