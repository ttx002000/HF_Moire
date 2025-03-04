using JLD2
include("../../src/operators_single_RMG_full_HF_noscreen.jl")
using LinearAlgebra,Plots
args=parse.(Float64,ARGS)
#args=[1.4,15.0,10.0,0.0,8.0,1.0,3.0,1.0,1.0]
radius=args[1]
num_points=Int(args[2])
uD=args[3]
ϵr=args[4]
NL=Int(args[5])
target_density=args[6]
temp=args[7]
trytime=args[8]
file_pos=args[9]


eig_set,k_set,k_index,eig_vec_set,Area,single_matrix=get_single_particle(radius,num_points,uD,NL)

println("finish1")



initial_density_matrix, BG_density_matrix=get_initial_proj_withansatz(k_set,eig_vec_set,NL)
HF_eigenvalues,HF_eigenvectors,energy, DIIS_input_density_matrix,fermi_level,Hartree_matrix,Fock_matrix,eout,renormalized_density=iteration(initial_density_matrix,BG_density_matrix,ϵr,k_set,single_matrix,Area,NL,target_density,temp)
scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "single_RMG_full_HF_noscreen/data_output$(Int(args[9]))/$(args[1])radius$(args[2])num_points$(args[3])uD$(args[4])er$(args[5])NL$(args[6])tgden$(args[7])temp$(args[8])trytime.jld2")

  
jldsave(savepath,
             final_density_matrix=DIIS_input_density_matrix[1],HF_eigenvalues=HF_eigenvalues,
             k_set=k_set,eig_vec_set=eig_vec_set,HF_eigenvectors=HF_eigenvectors,fermi_level=fermi_level,energy=energy,eig_set=eig_set,
             single_matrix=single_matrix,Hartree_matrix=Hartree_matrix,Fock_matrix=Fock_matrix,k_index=k_index,eout=eout,renormalized_density=renormalized_density)
