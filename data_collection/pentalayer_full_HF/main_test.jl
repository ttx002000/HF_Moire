using JLD2
include("../../src/operators_pentalayer_full_HF.jl")
using LinearAlgebra,Plots
args=parse.(Float64,ARGS)
#args=[1,1,1.4,20.0,30.0,10.0,1.0,5.0,1.0,-40.0,1.0]
#args=[1.4,20.0,10.0,0.02,8.0,1.0,1.0]
radius=args[1]
num_points=Int(args[2])
uD=args[3]
target_density=args[4]
ϵr=args[5]
temp=args[6]
trytime=args[7]



eig_set,k_set,k_index,eig_vec_set,Area,single_matrix=get_single_particle(radius,num_points,uD)

println("finish1")



initial_density_matrix, BG_density_matrix=get_initial_proj(k_set,eig_vec_set)
HF_eigenvalues,HF_eigenvectors,energy, DIIS_input_density_matrix,fermi_level,renormalized_density,Hartree_matrix,Fock_matrix=iteration(initial_density_matrix,BG_density_matrix,ϵr,k_set,single_matrix,target_density,temp)

scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "pentalayer_full_HF/data_output1/$(args[1])radius$(args[2])num_points$(args[3])uD$(args[4])tg_density$(args[5])er$(args[6])temp$(args[7])trytime.jld2")

  
jldsave(savepath,
             final_density_matrix=DIIS_input_density_matrix[1],HF_eigenvalues=HF_eigenvalues,
             k_set=k_set,eig_vec_set=eig_vec_set,HF_eigenvectors=HF_eigenvectors,fermi_level=fermi_level,energy=energy,eig_set=eig_set,
             single_matrix=single_matrix,Hartree_matrix=Hartree_matrix,Fock_matrix=Fock_matrix,k_index=k_index,renormalized_density=renormalized_density)
