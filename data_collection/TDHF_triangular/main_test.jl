using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames


include("../../src/operators_TDHF.jl")

args=parse.(Float64,ARGS)
flux=args[1]*π
V0=args[2]
ϕ=args[3]/180*π
Nq=Int(args[4]);
scale=args[5];
constq=args[6]/Nq^2
trytimes=Int(args[7])
DIIS_input_DensityMatrix=
loop_dic=
num_bandup=25
num_bandbelow=1
tot_bd=num_bandup+num_bandbelow

overlapmatrix, wave,wave_diff,single_MoirePo, single_Ham,allowedq, T1, T2=single_part(flux,V0,ϕ,scale,Nq)
form_overlapmatrix=get_foverlap(wave_diff,wave,allowedq,Nq,T1,T2)
Flevel,HF_eigenvector,HF_eigenvalue=get_HFeigenvectors(loop_dic,allowedq,T1,T2,Nq,wave,DIIS_input_DensityMatrix,single_Ham,single_MoirePo,constq,overlapmatrix)
Bandvector,Aindexset,AmQindexset,B2indexset=
 


