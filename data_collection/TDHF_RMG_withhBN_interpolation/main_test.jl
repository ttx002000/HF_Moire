using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames
using LinearAlgebra

#This is designed for triangular_R5G_contactinteraction_interpolation_withhBN
include("../../src/operators_TDHF_RMG_withhBN_interpolation.jl")


args=parse.(Float64,ARGS)
NL=Int(args[1])

Nq=Int(args[2]);
θ=args[3]/180*pi;
contact_strength=args[4]
ϵr=args[5]
uD=args[6]
filling=Int(args[7])
gcutoff=args[8]
λ=args[9]
bigQindex=Int(args[10])
filepos=Int(args[11])




scratch_dir = ENV["SCRATCH"]
input_path=joinpath(scratch_dir, "triangle_R5G_contact_interpolation_withhBN/data_output$(Int(args[11]))/TDHF_seed/$(args[1])NL$(args[2])Nq$(args[3])theta$(args[4])constq$(args[5])ϵr$(args[6])uD$(args[7])filling$(args[8])cutoff$(args[9])lambda_seed.jld2")
#input_path=joinpath(@__DIR__,"data.jld2")
println("found it")

st=load(input_path)
T1=st["T1"]
T2=st["T2"]
b1=st["b1"]
b2=st["b2"]
allowedq=st["allowedq"]
wave=st["wave"]
bigQ=allowedq[bigQindex]
num_bandup=length(wave)-1
num_bandbelow=1
tot_bd=length(wave)
HF_eigenvalue=st["HFeigenvalue"]
HF_eigenvector=st["HF_eigenvector"]

spinor_set=st["spinor_set"]
am=norm(st["a1m"])
wave_diff=get_wavediff(gcutoff,T1,T2,b1,b2)



  Area=(√3/2*am^2*Nq^2)


Bandvector,Aindexset,AmQindexset,B2indexset,gkpqmap,gkmqmap=get_indexset(wave_diff,num_bandbelow,num_bandup,allowedq,Nq,bigQ)

println("I am here")
Fmatrix=get_Fmatrix(Bandvector,spinor_set,
                          HF_eigenvector,
                          wave,wave_diff,NL)

println("I am here")
flush(stdout)
Amatrix,AmQmatrix,Bmatrix=Construct_Amatrix(Aindexset,AmQindexset,B2indexset,
                  allowedq,Fmatrix,gkpqmap,gkmqmap,
                  T1,T2,wave_diff,
                  HF_eigenvalue,Nq,Area,contact_strength,ϵr)
   println("I am here")
   flush(stdout)


 
 Totalmatrix=vcat(hcat(Amatrix,Bmatrix),hcat(-Bmatrix',-conj(AmQmatrix)))
 Smatrix=vcat(hcat(Amatrix,Bmatrix),hcat(Bmatrix',conj(AmQmatrix)))
 FFF=eigen(Totalmatrix)
 ω=FFF.values
 XYvector=FFF.vectors
 Sspectrum=eigvals(Smatrix)
 Aspectrum=eigvals(Amatrix)

output_path=joinpath(scratch_dir, "triangle_R5G_contact_interpolation_withhBN/data_output$(Int(args[11]))/TDHF_result/$(args[1])NL$(args[2])Nq$(args[3])theta$(args[4])constq$(args[5])ϵr$(args[6])uD$(args[7])filling$(args[8])cutoff$(args[9])lambda$(args[10])bigQ.jld2")


 jldsave(output_path,
         omegaspectrum=ω,
         Sspectrum=Sspectrum,
         Aspectrum=Aspectrum,Amatrix=Amatrix,Bmatrix=Bmatrix,AmQmatrix=AmQmatrix,
         Bandvector=Bandvector,Aindexset=Aindexset,AmQindexset=AmQindexset,B2indexset=B2indexset,XYvector= XYvector)


