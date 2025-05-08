using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames
using LinearAlgebra


include("../../src/operators_TDHF_exciton.jl")


args=parse.(Float64,ARGS)
#args=[1,1,1,1,1]

seed_num=Int(args[1])
bigQ_index=Int(args[2])
num_bandup=Int(args[3])
num_bandbelow=Int(args[4])
filepos=Int(args[5])


scratch_dir = ENV["SCRATCH"]

loadpath=joinpath(scratch_dir, "TDHF_exciton/data_output$(Int(args[5]))/data_input/seed$(Int(args[1])).jld2")
#loadpath=joinpath(@__DIR__,"seed1.jld2")

input=load(loadpath)

HF_eigenvalues=input["HF_eigenvalues"]
HF_eigenvectors=input["HF_eigenvectors"]

ϵr=input["ϵr"]
TDHF_Area=input["TDHF_Area"]
bigQ_index_set=input["bigQ_index_set"]


ildis=input["ildis"]
NL=Int(input["NL"])
bigQ=bigQ_index_set[bigQ_index]
println(bigQ)


TDHF_k_set=input["TDHF_k_set"]
TDHF_k_index=input["TDHF_k_index"]
TDHF_k_pos=input["TDHF_k_pos"]

k_index=input["k_index"]
k_set=input["k_set"]
JH=input["JH"]




tot_bd=num_bandup+num_bandbelow



Bandvector,Aindexset,AmQindexset,B2indexset=get_indexset(num_bandbelow,num_bandup,bigQ,TDHF_k_set,TDHF_k_index,TDHF_k_pos,k_index,k_set)
println(length(Aindexset))
println("finish index set")
flush(stdout)

#fcmatrix=get_formfactors(NL,k_set,ildis)
#println("I am here")
#flush(stdout)
Amatrix,AmQmatrix,Bmatrix=Construct_Amatrix(Aindexset,AmQindexset,B2indexset,HF_eigenvalues,HF_eigenvectors,TDHF_Area,ϵr,ildis,JH,NL,k_set)

 
 
 Totalmatrix=vcat(hcat(Amatrix,Bmatrix),hcat(-Bmatrix',-conj(AmQmatrix)))
 Smatrix=vcat(hcat(Amatrix,Bmatrix),hcat(Bmatrix',conj(AmQmatrix)))
 ω=eigvals(Totalmatrix)
 Sspectrum=eigvals(Smatrix)
 Aspectrum=eigvals(Amatrix)

savepath=joinpath(scratch_dir, "TDHF_exciton/data_output$(Int(args[5]))/seed$(args[1])bigQ$(args[2])nup$(args[3])ndown$(args[4]).jld2")
#savepath="test.jld2"

 jldsave(savepath,omegaspectrum=ω,Sspectrum=Sspectrum,Aspectrum=Aspectrum)


