using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames
using LinearAlgebra


include("../../src/operators_TDHF_RMG_bandprojected.jl")


args=parse.(Float64,ARGS)
#args=[5.0,0.77,6,50.0,7.0,1.0,0.0,3.01,5,1,1,1]
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
bigQindex=Int(args[11])
filepos=Int(args[12])



scratch_dir = ENV["SCRATCH"]
input_path=joinpath(scratch_dir, "triangle_R5G_bandproject_interpolation_withhBN/data_output$(Int(args[12]))/TDHF_seed/$(args[1])er$(args[2])theta$(args[3])Nq$(args[4])uD$(args[5])Nband$(args[6])lambda$(args[7])contact$(args[8])cutoff$(args[9])NL$(args[10])couplingratio_seed.jld2")


st=load(input_path)
T1=st["T1"]
T2=st["T2"]
allowedq=st["allowedq"]
wave=st["wave"]
wave_diff=st["wave_diff"]
bigQ=allowedq[bigQindex]
num_bandup=Nband-1
num_bandbelow=1
tot_bd=Nband
HF_eigenvalue=st["HF_eigenvalue"]
HF_eigenvector=st["HF_eigenvector"]
moire_eigenvector_single=st["eigenvector"]


 g1=T1*Nq
 g2=T2*Nq
 a1=inv([g1';g2'])*[2π,0]
 am=norm(a1)
 Area=(√3/2*am^2*Nq^2)


Bandvector,Aindexset,AmQindexset,B2indexset,gkpqmap,gkmqmap=get_indexset(wave_diff,num_bandbelow,num_bandup,allowedq,Nq,bigQ)
 Fmatrix=get_Fmatrix(Bandvector,moire_eigenvector_single,
                          HF_eigenvector,
                          wave,wave_diff,NL)


 Amatrix,AmQmatrix,Bmatrix=Construct_Amatrix(Aindexset,AmQindexset,B2indexset,
                  allowedq,Fmatrix,gkpqmap,gkmqmap,
                  T1,T2,wave_diff,
                  HF_eigenvalue,Nq,ϵr,Area,contact_strength)
   

 
 Totalmatrix=vcat(hcat(Amatrix,Bmatrix),hcat(-Bmatrix',-conj(AmQmatrix)))
 Smatrix=vcat(hcat(Amatrix,Bmatrix),hcat(Bmatrix',conj(AmQmatrix)))
 ω=eigvals(Totalmatrix)
 Sspectrum=eigvals(Smatrix)
 Aspectrum=eigvals(Amatrix)

output_path=joinpath(scratch_dir, "triangle_R5G_bandproject_interpolation_withhBN/data_output$(Int(args[12]))/TDHF_result/$(args[1])er$(args[2])theta$(args[3])Nq$(args[4])uD$(args[5])Nband$(args[6])lambda$(args[7])contact$(args[8])cutoff$(args[9])NL$(args[10])couplingratio$(args[11])bigQ.jld2")


 jldsave(output_path,
         omegaspectrum=ω,Sspectrum=Sspectrum,Aspectrum=Aspectrum)


