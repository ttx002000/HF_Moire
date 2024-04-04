using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames
using LinearAlgebra


include("../../src/operators_TDHF.jl")
input=load(joinpath(@__DIR__, "data_input/$(args[4])Nq$(args[1])flux$(args[2])V0$(args[3])phi$(args[5])scale$(args[6])constq_$(args[7])try.jld2"))

flux=0.0
V0=0.0
ϕ=0.0
Nq=3;
scale=1.0;
constq=1.0/Nq^2
DIIS_input_DensityMatrix=input["DIIS_input_DensityMatrix"][1]
loop_dic=input["loop_dic"]
bigQ=[0,0]
num_bandup=25
num_bandbelow=1
tot_bd=num_bandup+num_bandbelow

overlapmatrix, wave,wave_diff,single_MoirePo, single_Ham,allowedq,T1,T2=single_part(flux,V0,ϕ,scale,Nq)
form_overlapmatrix=get_foverlap(wave_diff,wave,allowedq,Nq,T1,T2,flux)
Flevel,HF_eigenvector,HF_eigenvalue=get_HFeigenvectors(loop_dic,allowedq,T1,T2,Nq,wave,DIIS_input_DensityMatrix,single_Ham,single_MoirePo,constq,overlapmatrix)
Bandvector,Aindexset,AmQindexset,B2indexset,gkpqmap,gkmqmap=get_indexset(Flevel,HF_eigenvalue,wave,wave_diff,num_bandbelow,num_bandup,allowedq,Nq,bigQ)
Fmatrix=get_Fmatrix(Bandvector,HF_eigenvector,wave,wave_diff,form_overlapmatrix)
Amatrix,AmQmatrix,Bmatrix=Construct_Amatrix(Aindexset,AmQindexset,B2indexset,allowedq,Fmatrix,gkpqmap,gkmqmap,T1,T2,wave_diff,HF_eigenvalue)
Totalmatrix=vcat(hcat(Amatrix,Bmatrix),hcat(-Bmatrix',-conj(AmQmatrix)))
ω=eigvals(Totalmatrix)





