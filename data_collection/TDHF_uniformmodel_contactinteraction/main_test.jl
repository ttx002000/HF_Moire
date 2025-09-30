using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames
using LinearAlgebra

#This is designed for triangular_uniformmodel_contactinteraction
include("../../src/operators_TDHF_uniformmodel_contactinteraction.jl")


args=parse.(Float64,ARGS)

flux=(args[1])*π
V0=args[2]
ϕ=args[3]/180*π
Nq=Int(args[4]);
scale=args[5];
contact_strength=args[6]
filling=Int(args[7])
gcutoff=args[8]
bigQindex=Int(args[9])
filepos=Int(args[10])



scratch_dir = ENV["SCRATCH"]
input_path=joinpath(scratch_dir, "triangle_uniformmodel_contact/data_output$(Int(args[10]))/TDHF_seed/$(args[1])flux$(args[2])V0$(args[3])phi$(args[4])Nq$(args[5])scale$(args[6])constq$(args[7])filling$(args[8])cutoff_seed.jld2")
#input_path=joinpath(@__DIR__,"data.jld2")
println("found it")

st=load(input_path)
T1=st["T1"]
T2=st["T2"]
allowedq=st["allowedq"]
wave=st["wave"]
bigQ=allowedq[bigQindex]
num_bandup=length(wave)-1
num_bandbelow=1
tot_bd=length(wave)
HF_eigenvalue=st["HFeigenvalue"]
HF_eigenvector=st["HF_eigenvector"]




wave_diff=get_wavediff(gcutoff,T1,T2,scale)


 am=4π/(√3*scale)
  Area=(√3/2*am^2*Nq^2)


Bandvector,Aindexset,AmQindexset,B2indexset,gkpqmap,gkmqmap=get_indexset(wave_diff,num_bandbelow,num_bandup,allowedq,Nq,bigQ)

println("I am here")
flush(stdout)
form_overlapmatrix=get_foverlap(wave_diff,wave,allowedq,Nq,T1,T2,flux,scale)
println("I am here")
flush(stdout)

Fmatrix=get_Fmatrix(Bandvector,HF_eigenvector,wave,wave_diff,T1,T2,flux,scale,form_overlapmatrix)
 

println("I am here")
flush(stdout)
Amatrix,AmQmatrix,Bmatrix=Construct_Amatrix(Aindexset,AmQindexset,B2indexset,
                  allowedq,Fmatrix,gkpqmap,gkmqmap,
                  T1,T2,wave_diff,
                  HF_eigenvalue,Nq,Area,contact_strength)
   println("I am here")
   flush(stdout)


 
 Totalmatrix=vcat(hcat(Amatrix,Bmatrix),hcat(-Bmatrix',-conj(AmQmatrix)))
 Smatrix=vcat(hcat(Amatrix,Bmatrix),hcat(Bmatrix',conj(AmQmatrix)))
 FFF=eigen(Totalmatrix)
 ω=FFF.values
 XYvecotrs=FFF.vectors
 Sspectrum=eigvals(Smatrix)
 Aspectrum=eigvals(Amatrix)

output_path=joinpath(scratch_dir, "triangle_uniformmodel_contact/data_output$(Int(args[10]))/TDHF_result/$(args[1])flux$(args[2])V0$(args[3])phi$(args[4])Nq$(args[5])scale$(args[6])constq$(args[7])filling$(args[8])cutoff$(args[9])bigQ.jld2")


 jldsave(output_path,
         omegaspectrum=ω,
         Sspectrum=Sspectrum,XYvecotrs=XYvecotrs,
         Aspectrum=Aspectrum,Amatrix=Amatrix,Bmatrix=Bmatrix,AmQmatrix=AmQmatrix,
         Bandvector=Bandvector,Aindexset=Aindexset,AmQindexset=AmQindexset,B2indexset=B2indexset)


