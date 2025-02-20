using JLD2
include("../../src/operators_double_RMG_full_HF_withmagneticfield.jl")
using LinearAlgebra,Plots
args=parse.(Float64,ARGS)

radius=args[1]
num_points=Int(args[2])
uD=args[3]
CNP=args[4]
ϵr=args[5]
ildis=args[6]
NL=Int(args[7])
target_density=args[8]
temp=args[9]
Bfield=args[10]
trytime=args[11]
file_pos=args[12]


eig_set,k_set,k_index,eig_vec_set,Area,single_matrix=get_single_particle(radius,num_points,uD,CNP,NL)

println("finish1")

magnetic_single_matrix=get_magnetic_single(k_set,eig_vec_set,single_matrix,Bfield,uD,NL)


initial_density_matrix, BG_density_matrix=get_magnetic_initial_proj(k_set,magnetic_single_matrix,NL)

HF_eigenvalues,HF_eigenvectors,energy, DIIS_input_density_matrix,fermi_level,Hartree_matrix,Fock_matrix,eout,renormalized_density=iteration(initial_density_matrix,BG_density_matrix,ϵr,k_set,magnetic_single_matrix,ildis,Area,NL,target_density,temp)
scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "double_RMG_full_HF_withmagneticfield_v2/data_output$(Int(args[12]))/$(args[1])radius$(args[2])num_points$(args[3])uD$(args[4])CNP$(args[5])er$(args[6])ildis$(args[7])NL$(args[8])tgden$(args[9])temp$(args[10])B$(args[11])trytime.jld2")

  
jldsave(savepath,
             final_density_matrix=DIIS_input_density_matrix[1],HF_eigenvalues=HF_eigenvalues,
             k_set=k_set,eig_vec_set=eig_vec_set,HF_eigenvectors=HF_eigenvectors,fermi_level=fermi_level,energy=energy,eig_set=eig_set,
             single_matrix=single_matrix,Hartree_matrix=Hartree_matrix,Fock_matrix=Fock_matrix,k_index=k_index,eout=eout,
             renormalized_density=renormalized_density,magnetic_single_matrix=magnetic_single_matrix)
