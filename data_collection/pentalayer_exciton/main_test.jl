using JLD2
include("../../src/operators_pentalayer_exciton.jl")
using LinearAlgebra
args=parse.(Float64,ARGS)
#args=[1,1,1.8,100.0,30.0,-160.0,5.0,1.0]
vone=Int(args[1])
vtwo=Int(args[2])
radius=args[3]
num_points=Int(args[4])
uD=args[5]
CNP=args[6]
ϵr=args[7]
trytime=args[8]

il_dis=CNP/uD*0.335

eig_set,k_set,k_index,eig_vec_set,Area,single_matrix=get_single_particle(vone,vtwo,radius,num_points,uD,CNP)
println("finish1")
tic=time()
formfactors=get_formfactors(k_set,eig_vec_set,il_dis)
toc=time()
println("finish2",toc-tic)
initial_density_matrix, BG_density_matrix=get_initial_proj(k_set)
HF_eigenvalues,HF_eigenvectors,energy, DIIS_input_density_matrix,fermi_level,Hartree_matrix,Fock_matrix=iteration(formfactors,initial_density_matrix,BG_density_matrix,ϵr,k_set,single_matrix,il_dis)


scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "pentalayer_graphene_exciton/uD=$(args[5])/vone$(vone)vtwo$(vtwo)/$(args[1])vone$(args[2])vtwo$(args[3])radius$(args[4])num_points$(args[5])uD$(args[6])CNP$(args[7])er$(args[8])trytime.jld2")

  
jldsave(savepath,
             DIIS_input_density_matrix=DIIS_input_density_matrix,HF_eigenvalues=HF_eigenvalues,
             k_set=k_set,eig_vec_set=eig_vec_set,HF_eigenvectors=HF_eigenvectors,fermi_level=fermi_level,energy=energy,eig_set=eig_set,
             single_matrix=single_matrix,Hartree_matrix=Hartree_matrix,Fock_matrix=Fock_matrix)
