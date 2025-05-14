using JLD2
include("../../src/operators_double_RMG_bandprojected_withscreening_withHunds.jl")
using LinearAlgebra,Plots
args=parse.(Float64,ARGS)
#args=[1.4, 61.0, 10.0, 24.0, 5.0, 0.002, 0.01, -166.0, 0.0, 1.0, 2.0, 2.0]
#args=[0.8,25.0,15.0,10.0,3.0,0.0,0.01,-500.0,5.0,1.5,0.0,1.0,1.0]
radius=args[1]
num_points=Int(args[2])
uD=args[3]
ϵr=args[4]
NL=Int(args[5])
target_density=args[6]
temp=args[7]
JH=args[8]
CNP=args[9]
ildis=args[10]
shift=Int(args[11])
trytime=args[12]
file_pos=args[13]


eig_set,k_set,k_index,eig_vec_set,Area,Fock_formfactors,Hartree_formfactors,Hunds_formfactors,single_matrix=get_single_particle(radius,num_points,uD,NL,ildis,shift)


println("finish1")




initial_density_matrix, BG_density_matrix=get_initial_proj(k_set)
HF_eigenvalues,HF_eigenvectors,energy, DIIS_input_density_matrix,fermi_level,Hartree_matrix,Fock_matrix,eout,renormalized_density=iteration(initial_density_matrix,BG_density_matrix,
                                                                                                                              ϵr,k_set,single_matrix,Area,
                                                                                                                               target_density,temp,Fock_formfactors,Hartree_formfactors
                                                                                                                               ,Hunds_formfactors,JH,0)




final_density_matrix=DIIS_input_density_matrix[1]
second_input=deepcopy(final_density_matrix)

DIIS_input_density_matrix=nothing
initial_density_matrix=nothing
GC.gc()
println("start second iterations")

parent_HF_eigenvalues,parent_HF_eigenvectors,parent_energy, parent_DIIS_input_density_matrix,parent_fermi_level,parent_Hartree_matrix,parent_Fock_matrix,_,_=iteration(second_input,BG_density_matrix,
                                                                                                                              ϵr,k_set,single_matrix,Area,
                                                                                                                               target_density,temp,Fock_formfactors,Hartree_formfactors
                                                                                                                               ,Hunds_formfactors,JH,-1)
println("finish second iterations")




 
#scratch_dir = ENV["SCRATCH"]
#savepath=joinpath(scratch_dir, "double_RMG_bandprojected_screen_withHunds/data_output$(Int(args[13]))/$(args[1])radius$(args[2])num_points$(args[3])uD$(args[4])er$(args[5])NL$(args[6])tgden$(args[7])temp$(args[8])JH$(args[9])CNP$(args[10])ildis$(args[11])shift$(args[12])trytime.jld2")

savepath=joinpath(@__DIR__,"test.jld2")
jldsave(savepath,
             final_density_matrix=final_density_matrix,HF_eigenvalues=HF_eigenvalues,
             k_set=k_set,eig_vec_set=eig_vec_set,HF_eigenvectors=HF_eigenvectors,fermi_level=fermi_level,energy=energy,eig_set=eig_set,
             single_matrix=single_matrix,Hartree_matrix=Hartree_matrix,Fock_matrix=Fock_matrix,k_index=k_index,eout=eout,renormalized_density=renormalized_density,
             parent_HF_eigenvalues=parent_HF_eigenvalues,parent_HF_eigenvectors=parent_HF_eigenvectors,parent_DS=parent_DIIS_input_density_matrix[1],
             parent_Hartree_matrix=parent_Hartree_matrix,parent_Fock_matrix=parent_Fock_matrix,parent_energy=parent_energy)
