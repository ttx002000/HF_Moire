using JLD2
include("../../src/operators_double_RMG_full_HF_allowcoherence.jl")
using LinearAlgebra,Plots
args=parse.(Float64,ARGS)
#args=[1.4,15.0,10.0,0.0,8.0,1.0,3.0,1.0,1.0]
radius=args[1]
num_points=Int(args[2])
uD=args[3]
CNP=args[4]
ϵr=args[5]
ildis=args[6]
NL=Int(args[7])
target_density=args[8]
temp=args[9]
trytime=args[10]
file_pos=args[11]


eig_set,k_set,k_index,eig_vec_set,Area,single_matrix=get_single_particle(radius,num_points,uD,CNP,NL)

println("finish1")



_, BG_density_matrix=get_initial_proj(k_set,eig_vec_set,NL)
scratch_dir = ENV["SCRATCH"]
seed_num=mod(Int(trytime),5)+1
savepath=joinpath(scratch_dir, "double_RMG_full_HF_allowcoherence/data_output$(Int(args[11]))/seeds/seeds$(seed_num).jld2")
st=load(savepath)

initial_density_matrix=st["DS_seed"]

HF_eigenvalues,HF_eigenvectors,energy, DIIS_input_density_matrix,fermi_level,Hartree_matrix,Fock_matrix,eout,renormalized_density=iteration_noDIIS(initial_density_matrix,BG_density_matrix,ϵr,k_set,single_matrix,ildis,Area,NL,target_density,temp)
scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "double_RMG_full_HF_allowcoherence/data_output$(Int(args[11]))/$(args[1])radius$(args[2])num_points$(args[3])uD$(args[4])CNP$(args[5])er$(args[6])ildis$(args[7])NL$(args[8])tgden$(args[9])temp$(args[10])trytime.jld2")

  
jldsave(savepath,
             final_density_matrix=DIIS_input_density_matrix[1],HF_eigenvalues=HF_eigenvalues,
             k_set=k_set,eig_vec_set=eig_vec_set,HF_eigenvectors=HF_eigenvectors,fermi_level=fermi_level,energy=energy,eig_set=eig_set,
             single_matrix=single_matrix,Hartree_matrix=Hartree_matrix,Fock_matrix=Fock_matrix,k_index=k_index,eout=eout,renormalized_density=renormalized_density,seed_num=seed_num)
