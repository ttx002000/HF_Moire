using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames
using LinearAlgebra


include(joinpath(@__DIR__,"../../src/operators_ansatz_physicalunits.jl"))

args=parse.(Float64,ARGS)

#args=[1.0,3.0,0.0,0.0,0.1,1.0,40*π^2/√3]


flux=args[1]*π
Nq=Int(args[2])
V0=args[3]
ϕ=args[4]
chistart=args[5]
chiend=args[6]
Vc=args[7]


β=4*π/(√3)
scale=sqrt(flux/π)
am=4π/(scale*√3)
Auc=√3/2*am^2
constq=Vc/(Nq^2*Auc)


chispace=collect(LinRange(chistart,chiend,200))#Notice the subtlety in the definition of chi

single_Ham, eigenvector, wave, wave_diff, allowedq,T1,T2=get_wavefunction(scale,flux,Nq,V0,ϕ,chispace,constq)

form_overlapmatrix=get_formoverlap(Nq,allowedq,wave,wave_diff,flux,T1,T2)
   
gkpqmap,gkmqmap=get_gkpgmap(Nq,allowedq,wave_diff)


Energy,kinetic,Fock=get_energy(chispace,constq,allowedq,form_overlapmatrix,single_Ham,eigenvector,wave,wave_diff,T1,T2,gkpqmap,gkmqmap)
      
jldsave(joinpath(@__DIR__, "data_output/physicalunits_$(args[1])flux$(args[2])Nq$(args[3])V0$(args[4])phi$(args[5])chist$(args[6])chien$(args[7])Vc.jld2"),chispace=chispace,Vc=Vc,Nq=Nq,scale=scale,Energy=Energy,Fock=Fock,kinectic=kinetic)