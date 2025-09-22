using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames
using LinearAlgebra

#This is designed for triangular_R5G_contactinteraction_v2
include("../../src/operators_TDHF_RMG_contactinteraction_v2.jl")


args=parse.(Float64,ARGS)
#args=[5.0,0.77,6,50.0,7.0,1.0,0.0,3.01,5,1,1,1]
args=[5.0,0.0,0.0,6.0,2.0,944060.8762859226,1.0,3.51,1.0,1.0]
NL=Int(args[1])
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
input_path=joinpath(scratch_dir, "triangle_R5G_contact_v2/data_output$(Int(args[10]))/TDHF_seed/$(args[1])NL$(args[2])V0$(args[3])phi$(args[4])Nq$(args[5])scale$(args[6])constq$(args[7])filling$(args[8])cutoff_seed.jld2")
#input_path=joinpath(@__DIR__,"data.jld2")
println("found it")

st=load(input_path)
T1=st["T1"]
T2=st["T2"]
allowedq=st["allowedq"]
wave=st["wave"]
#wave_diff=st["wave_diff"]
bigQ=allowedq[bigQindex]
num_bandup=length(wave)-1
num_bandbelow=1
tot_bd=length(wave)
HF_eigenvalue=st["HFeigenvalue"]
HF_eigenvector=st["HF_eigenvector"]
#moire_eigenvector_single=st["eigenvector"]
spinor_set=st["spinor_set"]



wave_diff=get_wavediff(gcutoff,T1,T2,scale)


 am=4π/(√3*scale)
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
                  HF_eigenvalue,Nq,Area,contact_strength)
   println("I am here")
   flush(stdout)


 
 Totalmatrix=vcat(hcat(Amatrix,Bmatrix),hcat(-Bmatrix',-conj(AmQmatrix)))
 Smatrix=vcat(hcat(Amatrix,Bmatrix),hcat(Bmatrix',conj(AmQmatrix)))
 ω=eigvals(Totalmatrix)
 Sspectrum=eigvals(Smatrix)
 Aspectrum=eigvals(Amatrix)

output_path=joinpath(scratch_dir, "triangle_R5G_contact_v2/data_output$(Int(args[10]))/TDHF_result/$(args[1])NL$(args[2])V0$(args[3])phi$(args[4])Nq$(args[5])scale$(args[6])constq$(args[7])filling$(args[8])cutoff$(args[9])bigQ.jld2")


 jldsave(output_path,
         omegaspectrum=ω,Sspectrum=Sspectrum,Aspectrum=Aspectrum)


