using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames

# The goal of this is to do HF band projection and also allow me to do HF-projected skyrmion excitation calculation
include("../../src/operators_R5G_contactinteraction_interpolation_withhBN_skyrmionexcitation_step2.jl")

args=parse.(Float64,ARGS)

NL=args[1]

Nq=Int(args[2]);
θ=args[3]/180*pi;
constq=args[4]
ϵr=args[5]
uD=args[6]
nop_filling=Int(args[7])# This needs to be changed
gcutoff=args[8]
λ=args[9]
trytimes=Int(args[10])
enlarge_factor=Int(args[11])
V0_hBN=args[12]
V1_hBN=args[13]
ψ_hBN=args[14]
V2_scalar=args[15]
ϕ=args[16]/180*π
qcutoff=args[17]
Nband=Int(args[18])
filling=Int(args[19])
filepos=Int(args[20])



scratch_dir = ENV["SCRATCH"]
seed_path=joinpath(scratch_dir, "triangle_R5G_contact_interpolation_withhBN_skyrmionexcitation/data_output$(Int(args[20]))/step2/seed/seed_$(args[1])NL$(args[2])Nq$(args[3])theta$(args[4])constq$(args[5])ϵr$(args[6])uD$(args[7])filling$(args[8])cutoff$(args[9])lambda$(args[11])enlarge$(args[12])V0_hBN$(args[13])V1_hBN$(args[14])ψ_hBN$(args[15])V2_scalar$(args[16])ϕ.jld2")

allowedq=seed_file["allowedq"]
spinor_set=seed_file["spinor_set"]
nop_HF_eigenvector=seed_file["HF_eigenvector"]
nop_HF_eigenvalue=seed_file["HFeigenvalue"]
nop_single_Ham=seed_file["single_Ham"]
nop_single_MoirePo=seed_file["single_MoirePo"]
nop_densitymatrix=seed_file["densitymatrix"]
nop_HartreeMatrix=seed_file["HartreeMatrix"]
nop_FockMatrix=seed_file["FockMatrix"]
wave=seed_file["wave"]
T1=seed_file["T1"]
T2=seed_file["T2"]
b1=seed_file["b1"]
b2=seed_file["b2"]
a1m=seed_file["a1m"]
a2m=seed_file["a2m"]
b1T=Int.(round.(inv([T1 T2])*b1))
b2T=Int.(round.(inv([T1 T2])*b2))



allowedq_dic,single_Ham,single_MoirePo,wave_diff=initial_process(allowedq,nop_single_Ham,
                         nop_single_MoirePo,nop_HF_eigenvector,
                         Nband,qcutoff,b1,b2,b1T,b2T)



Area=Nq^2*√3/2*norm(a1m)^2


DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,eout,HartreeMatrix,FockMatrix=iteration_loop(initial_DensityMatrix, allowedq,allowedq_dic,T1,T2,
                                                                                                                                    Nq,wave,wave_diff,single_Ham,
                                                                                                                                    single_MoirePo,constq,ϵr,formfactors,filling,Area)



scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "triangle_R5G_contact_interpolation_withhBN_skyrmionexcitation/data_output$(Int(args[20]))/step2/$(args[1])NL$(args[2])Nq$(args[3])theta$(args[4])constq$(args[5])ϵr$(args[6])uD$(args[7])filling$(args[8])cutoff$(args[9])lambda$(args[10])trytime$(args[11])enlarge$(args[12])V0_hBN$(args[13])V1_hBN$(args[14])ψ_hBN$(args[15])V2_scalar$(args[16])ϕ$(args[17])qcutoff$(args[18])Nband$(args[19])va_filling.jld2")



jldsave(savepath,single_Ham=single_Ham,single_MoirePo=single_MoirePo,
                spinor_set=spinor_set,arguments=args,
                densitymatrix=DIIS_input_DensityMatrix[1],energy=energy,eout=eout,
                HFeigenvalue=HF_eigenvalue,
                HartreeMatrix=HartreeMatrix,FockMatrix=FockMatrix,HF_eigenvector=HF_eigenvector,
                allowedq=allowedq,T1=T1,T2=T2,wave=wave,wave_diff=wave_diff,
                a1m=a1m,a2m=a2m,b1=b1,b2=b2,
                nop_densitymatrix=nop_densitymatrix,nop_HF_eigenvector=nop_HF_eigenvector,
                nop_HF_eigenvalue=nop_HF_eigenvalue,nop_single_Ham=nop_single_Ham,nop_single_MoirePo=nop_single_MoirePo,
                nop_FockMatrix=nop_FockMatrix,nop_HartreeMatrix=nop_HartreeMatrix)


