using JLD2
include("../../src/operators_pentalayer_exciton_screening.jl")
using LinearAlgebra
args=parse.(Float64,ARGS)
#args=[1,1,1.4,20.0,30.0,10.0,1.0,5.0,1.0,-40.0,1.0]
vone=Int(args[1])
vtwo=Int(args[2])
radius=args[3]
num_points=Int(args[4])
uD=args[5]
ildis=args[6]
ϵr=args[7]
stackingtwo=Int(args[8])
CNP=args[9]-4*uD
trytime=args[10]



eig_set,k_set,k_index,eig_vec_set,Area,single_matrix=get_single_particle(vone,vtwo,radius,num_points,uD,stackingtwo,CNP)

println("finish1")

#=
tic=time()
formfactors,Hartree_formfactors=get_formfactors(k_set,eig_vec_set,ildis)
toc=time()
println("finish2",toc-tic)

scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "pentalayer_exciton_phasediagram/FF_$(args[1])vone$(args[2])vtwo$(args[3])radius$(args[4])num_points$(args[5])uD$(args[6])ildis$(args[8])stackingtwo.jld2")
jldsave(savepath,formfactors=formfactors,Hartree_formfactors=Hartree_formfactors)
=#
scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "pentalayer_exciton_phasediagram/FF_$(args[1])vone$(args[2])vtwo$(args[3])radius$(args[4])num_points$(args[5])uD$(args[6])ildis$(args[8])stackingtwo.jld2")
st=load(savepath)
formfactors=st["formfactors"]
Hartree_formfactors=st["Hartree_formfactors"]

initial_density_matrix, BG_density_matrix=get_initial_proj(k_set)
HF_eigenvalues,HF_eigenvectors,energy, DIIS_input_density_matrix,fermi_level,Hartree_matrix,Fock_matrix=iteration(formfactors,initial_density_matrix,BG_density_matrix,ϵr,k_set,single_matrix,Hartree_formfactors)


scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "pentalayer_exciton_phasediagram/data_output1/uD=$(args[5])/vone$(vone)vtwo$(vtwo)/$(args[1])vone$(args[2])vtwo$(args[3])radius$(args[4])num_points$(args[5])uD$(args[6])ildis$(args[7])er$(args[8])stackingtwo$(args[9])CNP$(args[10])trytime.jld2")

  
jldsave(savepath,
             DIIS_input_density_matrix=DIIS_input_density_matrix,HF_eigenvalues=HF_eigenvalues,
             k_set=k_set,eig_vec_set=eig_vec_set,HF_eigenvectors=HF_eigenvectors,fermi_level=fermi_level,energy=energy,eig_set=eig_set,
             single_matrix=single_matrix,Hartree_matrix=Hartree_matrix,Fock_matrix=Fock_matrix,k_index=k_index)
