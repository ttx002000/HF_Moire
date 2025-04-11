using JLD2
include("../../src/operators_bandprojected_RMG.jl")
using LinearAlgebra,Plots
args=parse.(Float64,ARGS)
#args=[0.8,35.0,15.0,5.0,3.0,0.01,0.1,1.0,0.1,1.0,1000.0,1.0]
radius=args[1]
num_points=Int(args[2])
uD=args[3]
ϵr=args[4]
NL=Int(args[5])
target_density=args[6]
temp=args[7]
JH=args[8]
SOCcoef=args[9]
whichside=Int(args[10])

trytime=args[11]
file_pos=args[12]


eig_set,k_set,k_index,eig_vec_set,Area,formfactors,Coulombmatrix,single_matrix=get_single_particle(radius,num_points,uD,NL,whichside,SOCcoef)

println("finish1")




initial_density_matrix, BG_density_matrix=get_initial_proj(k_set,whichside)
HF_eigenvalues,HF_eigenvectors,energy, DIIS_input_density_matrix,fermi_level,Hartree_matrix,Fock_matrix,eout,renormalized_density=iteration(initial_density_matrix,BG_density_matrix,
                                                                                                                                  ϵr,k_set,single_matrix,Area,
                                                                                                                               target_density,temp,Coulombmatrix,formfactors,JH,whichside)


scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "bandprojected_RMG/data_output$(Int(args[12]))/$(args[1])radius$(args[2])num_points$(args[3])uD$(args[4])er$(args[5])NL$(args[6])tgden$(args[7])temp$(args[8])JH$(args[9])SOC$(args[10])whichside$(args[11])trytime.jld2")

  
jldsave(savepath,
             final_density_matrix=DIIS_input_density_matrix[1],HF_eigenvalues=HF_eigenvalues,
             k_set=k_set,eig_vec_set=eig_vec_set,HF_eigenvectors=HF_eigenvectors,fermi_level=fermi_level,energy=energy,eig_set=eig_set,
             single_matrix=single_matrix,Hartree_matrix=Hartree_matrix,Fock_matrix=Fock_matrix,k_index=k_index,eout=eout,renormalized_density=renormalized_density)
