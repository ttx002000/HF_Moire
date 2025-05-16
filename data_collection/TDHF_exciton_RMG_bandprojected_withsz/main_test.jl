using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames
using LinearAlgebra


include("../../src/operators_TDHF_exciton_withsz.jl")


args=parse.(Float64,ARGS)


seed_num=Int(args[1])
bigQ_index=Int(args[2])
num_bandup=Int(args[3])
num_bandbelow=Int(args[4])
mode=Int(args[5])
filepos=Int(args[6])


scratch_dir = ENV["SCRATCH"]

loadpath=joinpath(scratch_dir, "TDHF_exciton_withsz/data_output$(Int(args[6]))/data_input/seed$(Int(args[1])).jld2")
#loadpath=joinpath(@__DIR__,"seed1.jld2")

input=load(loadpath)

ϵr=input["ϵr"]
TDHF_Area=input["TDHF_Area"]
bigQ_index_set=input["bigQ_index_set"]

ildis=input["ildis"]

NL=input["NL"]
TDHF_k_set=input["TDHF_k_set"]
TDHF_k_index=input["TDHF_k_index"]
TDHF_k_pos=input["TDHF_k_pos"]
JH=input["JH"]
k_index=input["k_index"]
k_set=input["k_set"]

HF_eigenvalues_TDHF=input["HF_eigenvalues_TDHF"]

HF_eigenvectors_sublatticebasis_TDHF=input["HF_eigenvectors_sublatticebasis_TDHF"]
valley_num=2
layer_num=2
spin_num=2
sublattice_num=2*NL
band_num=valley_num*layer_num







tot_bd=num_bandup+num_bandbelow

z_pos=zeros(Float64,layer_num,sublattice_num)
z_pos[1,:]=0.335*[i for i in 0:NL-1 for _ in 1:2]
z_pos[2,:]=0.335*[i for i in -NL+1:0 for _ in 1:2].-ildis


bigQ=bigQ_index_set[bigQ_index]

println(bigQ)





Bandvector,Aindexset,AmQindexset,B2indexset=get_indexset(num_bandbelow,num_bandup,mode,bigQ,TDHF_k_set,TDHF_k_index,TDHF_k_pos,k_index,k_set)
println(length(Aindexset))
println("finish index set")
flush(stdout)


Amatrix,AmQmatrix,Bmatrix=Construct_Amatrix(Aindexset,AmQindexset,B2indexset,HF_eigenvalues_TDHF,HF_eigenvectors_sublatticebasis_TDHF,TDHF_Area,ϵr,NL,k_set,JH,z_pos)

 
 
 Totalmatrix=vcat(hcat(Amatrix,Bmatrix),hcat(-Bmatrix',-conj(AmQmatrix)))
 Smatrix=vcat(hcat(Amatrix,Bmatrix),hcat(Bmatrix',conj(AmQmatrix)))
 FFF=eigen(Totalmatrix)
 ω=FFF.values
 XYvectors=FFF.vectors
 Sspectrum=eigvals(Smatrix)
 #Aspectrum=eigvals(Amatrix)

savepath=joinpath(scratch_dir, "TDHF_exciton_withsz/data_output$(Int(args[6]))/seed$(args[1])bigQ$(args[2])nup$(args[3])ndown$(args[4])mode$(args[5]).jld2")
#savepath="test.jld2"

 jldsave(savepath,omegaspectrum=ω,Sspectrum=Sspectrum,Aindexset=Aindexset,Bandvector=Bandvector,
         AmQindexset=AmQindexset,B2indexset=B2indexset,Amatrix=Amatrix,AmQmatrix=AmQmatrix,Bmatrix=Bmatrix,XYvectors=XYvectors)


