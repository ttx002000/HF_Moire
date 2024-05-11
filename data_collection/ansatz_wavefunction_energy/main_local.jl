using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames
using LinearAlgebra

include(joinpath(@__DIR__,"../../src/operators_ansatz.jl"))

args=[1.0,1.0,3.0,0.0,0.0,0.0,1.0,1.0]
scale=args[1]
flux=args[2]*π
Nq=Int(args[3])
V0=args[4]
ϕ=args[5]
chistart=args[6]
chiend=args[7]
constq=args[8]/Nq^2

chispace=collect(LinRange(0.2, 2.0, 30)) 

single_Ham, eigenvector, wave, wave_diff, allowedq,T1,T2=get_wavefunction(scale,flux,Nq,V0,ϕ,chispace,constq)

form_overlapmatrix=get_formoverlap(Nq,allowedq,wave,wave_diff,flux, T1,T2)
   
gkpqmap,gkmqmap=get_gkpgmap(Nq,allowedq,wave_diff)


Energy,kinetic,Fock=get_energy(chispace,constq,allowedq,form_overlapmatrix,single_Ham,eigenvector,wave,wave_diff,T1,T2,gkpqmap,gkmqmap)
   
jldsave(joinpath(@__DIR__, "data_output/$(args[1])scale$(args[2])flux$(args[3])Nq$(args[4])V0$(args[5])phi$(args[6])chist$(args[7])chien$(args[8])Cq.jld2"),chispace=chispace,constq=constq,Nq=Nq,scale=scale,Energy=Energy,Fock=Fock,kinectic=kinetic)