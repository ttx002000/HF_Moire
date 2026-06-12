using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames

# same  as v2, but for arbitrary geometry
include("../../src/operators_R5G_contactinteraction_interpolation_withhBN_skyrmionexcitation_v5.jl")

args=parse.(Float64,ARGS)

NL=args[1]

geonum=Int(args[2]);
θ=args[3]/180*pi;
constq=args[4]
ϵr=args[5]
uD=args[6]
filling=Int(args[7])
gcutoff=args[8]
λ=args[9]
trytimes=Int(args[10])
V0_hBN=args[11]
V1_hBN=args[12]
ψ_hBN=args[13]
V2_scalar=args[14]
ϕ=args[15]/180*π
gateD=args[16]
filepos=Int(args[17])






overlapmatrix, wave, initial_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,single_eigenvector,allowedq, T1, T2, a1m, a2m, b1,b2,spinor_set,Area,b1T,b2T,_,_,_,_,_,_=triangle_initial_Densitymatrix(Int(NL),θ,geonum,gcutoff,uD,λ,V0_hBN,V1_hBN,ψ_hBN,V2_scalar,ϕ)


scratch_dir = ENV["SCRATCH"]

seed_path=joinpath(scratch_dir, "triangle_R5G_contact_interpolation_withhBN_skyrmionexcitation_v5/data_output$(Int(args[19]))/seed")
if only(rand())>0.0
   seed_file_path=pick_random_jld2_path(seed_path)
   if !(seed_file_path==nothing)
      seed_file=load(seed_file_path)
      initial_DensityMatrix=transform_dm(seed_file["densitymatrix"],wave,allowedq,seed_file["spinor_set"],spinor_set)
      println("Using seed from $seed_file_path")
   end
end





     
DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,eout,HartreeMatrix,FockMatrix=iteration_loop(initial_DensityMatrix,
                                                    allowedq,T1,T2,wave,single_Ham,single_MoirePo,constq,ϵr,overlapmatrix,filling,Area,gateD)


chern,Flink,chern_single,Flink_single,uniform,uniform_single=triangle_chern(wave,allowedq,HF_eigenvector,single_eigenvector,spinor_set,geonum,b1T,b2T,Int(NL),filling)




scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "triangle_R5G_contact_interpolation_withhBN_skyrmionexcitation_v5/data_output$(Int(args[17]))/$(args[1])NL$(args[2])geonum$(args[3])theta$(args[4])constq$(args[5])ϵr$(args[6])uD$(args[7])filling$(args[8])cutoff$(args[9])lambda$(args[10])trytime$(args[11])V0_hBN$(args[12])V1_hBN$(args[13])ψ_hBN$(args[14])V2_scalar$(args[15])ϕ$(args[16])gateD.jld2")

#savepath="test.jld2"


jldsave(savepath,single_Ham=single_Ham,single_MoirePo=single_MoirePo,
               single_eigenvalue=single_eigenvalue,spinor_set=spinor_set,
                chern=chern,Flink=Flink,chern_single=chern_single,
                Flink_single=Flink_single,arguments=args,
                densitymatrix=DIIS_input_DensityMatrix[1],energy=energy,eout=eout,
                HFeigenvalue=HF_eigenvalue,
                uniform=uniform,uniform_single=uniform_single,
                HartreeMatrix=HartreeMatrix,FockMatrix=FockMatrix,HF_eigenvector=HF_eigenvector,
                single_eigenvector=single_eigenvector,allowedq=allowedq,T1=T1,T2=T2,wave=wave,
                a1m=a1m,a2m=a2m,b1=b1,b2=b2)


