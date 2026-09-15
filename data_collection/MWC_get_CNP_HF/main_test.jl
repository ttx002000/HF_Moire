using Pkg
Pkg.activate(joinpath(@__DIR__,"../.."))
using JLD2
include("../../src/MWC_get_CNP_HF.jl")
using LinearAlgebra
args=parse.(Float64,ARGS)
radius=args[1]
num_points=Int(args[2])
uD=args[3]
ϵr=args[4]
NL=Int(args[5])
gatedis=args[6]
coherence=Int(args[7])
trytime=args[8]
file_pos=args[9]


eig_set,k_set,k_index,eig_vec_set,Area,single_matrix=get_single_particle(radius,num_points,uD,NL)

println("finish1")


initial_density_matrix, BG_density_matrix=get_initial_proj(k_set,eig_vec_set,NL)
HF_eigenvalues,HF_eigenvectors,energy, DIIS_input_density_matrix,fermi_level,Hartree_matrix,Fock_matrix,eout=iteration(initial_density_matrix,BG_density_matrix,ϵr,k_set,single_matrix,Area,NL,gatedis,coherence)
scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "MWC_CNP/data_output$(Int(args[9]))/$(args[1])radius$(args[2])num_points$(args[3])uD$(args[4])er$(args[5])NL$(args[6])gatedis$(args[7])coh$(args[8])trytime.jld2")
mkpath(dirname(savepath))
jldsave(savepath,
             final_density_matrix=DIIS_input_density_matrix[1],HF_eigenvalues=HF_eigenvalues,
             k_set=k_set,eig_vec_set=eig_vec_set,HF_eigenvectors=HF_eigenvectors,fermi_level=fermi_level,energy=energy,eig_set=eig_set,
             single_matrix=single_matrix,Hartree_matrix=Hartree_matrix,Fock_matrix=Fock_matrix,k_index=k_index,eout=eout)
