using JLD2
include("../../src/operators_double_RMG_full_HF_allowcoherence_constrainedpairing.jl")
using LinearAlgebra,Plots
args=parse.(Float64,ARGS)
println("here is my parameters",args)
flush(stdout)
#args=[1.4,15.0,10.0,20.0,8.0,1.0,5.0,0.0,0.1,1.0,1.0,1.0]
radius=args[1]
num_points=Int(args[2])
uD=args[3]
CNP=args[4]
ϵr=args[5]
ildis=args[6]
NL=Int(args[7])
target_density=args[8]
temp=args[9]
pairing=Int(args[10])
trytime=args[11]
file_pos=args[12]


eig_set,k_set,k_index,eig_vec_set,Area,single_matrix=get_single_particle(radius,num_points,uD,CNP,NL)

println("finish1")



#initial_density_matrix, BG_density_matrix=get_initial_proj(k_set,eig_vec_set,NL,pairing)

rd, BG_density_matrix=get_initial_proj(k_set,eig_vec_set,NL,2)

scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "double_RMG_full_HF_allowcoherence_constrainedpairing/data_output$(Int(args[12]))/seeds/seed$(mod(Int(trytime),1)+1).jld2")
st=load(savepath)

initial_density_matrix=st["parent_DS"]+rd*10^(-3)




HF_eigenvalues,HF_eigenvectors,energy, DIIS_input_density_matrix,fermi_level,Hartree_matrix,Fock_matrix,eout,renormalized_density=iteration(initial_density_matrix,BG_density_matrix,ϵr,k_set,single_matrix,ildis,Area,NL,target_density,temp,pairing)
println("finish first iterations")

final_density_matrix=DIIS_input_density_matrix[1]
second_input=deepcopy(final_density_matrix)

DIIS_input_density_matrix=nothing
initial_density_matrix=nothing
GC.gc()
println("start second iterations")

parent_HF_eigenvalues,parent_HF_eigenvectors,_, parent_DIIS_input_density_matrix,parent_fermi_level,parent_Hartree_matrix,parent_Fock_matrix,_,_=iteration(second_input,BG_density_matrix,ϵr,k_set,single_matrix,ildis,Area,NL,target_density,temp,-1)

println("finish second iterations")






scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "double_RMG_full_HF_allowcoherence_constrainedpairing/data_output$(Int(args[12]))/$(args[1])radius$(args[2])num_points$(args[3])uD$(args[4])CNP$(args[5])er$(args[6])ildis$(args[7])NL$(args[8])tgden$(args[9])temp$(args[10])pairing$(args[11])trytime.jld2")

  
jldsave(savepath,
             final_density_matrix=final_density_matrix,HF_eigenvalues=HF_eigenvalues,
             k_set=k_set,eig_vec_set=eig_vec_set,HF_eigenvectors=HF_eigenvectors,fermi_level=fermi_level,energy=energy,eig_set=eig_set,
             single_matrix=single_matrix,Hartree_matrix=Hartree_matrix,Fock_matrix=Fock_matrix,k_index=k_index,eout=eout,renormalized_density=renormalized_density,
             parent_HF_eigenvalues=parent_HF_eigenvalues,parent_HF_eigenvectors=parent_HF_eigenvectors,parent_DS=parent_DIIS_input_density_matrix[1],
             parent_Hartree_matrix=parent_Hartree_matrix,parent_Fock_matrix=parent_Fock_matrix)
