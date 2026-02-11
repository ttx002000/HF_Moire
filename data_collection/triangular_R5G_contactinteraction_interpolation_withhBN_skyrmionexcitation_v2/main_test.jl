using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames

# The goal of this script is to do spinless calculation, and also allow me to do spinful HF-projected skyrmion excitation calculation 
include("../../src/operators_R5G_contactinteraction_interpolation_withhBN_skyrmionexcitation_v2.jl")

args=parse.(Float64,ARGS)
#args=[5,3,0.6,30.0*10^4,5,35.0,1,3.51,1.0,1,1.0,28.9,21.0,-0.29,5.0,0.0,1.0]


NL=args[1]

Nq=Int(args[2]);
θ=args[3]/180*pi;
constq=args[4]
ϵr=args[5]
uD=args[6]
filling=Int(args[7])
gcutoff=args[8]
λ=args[9]
trytimes=Int(args[10])
enlarge_factor=Int(args[11])
V0_hBN=args[12]
V1_hBN=args[13]
ψ_hBN=args[14]
V2_scalar=args[15]
ϕ=args[16]/180*π
filepos=Int(args[17])






overlapmatrix, wave, initial_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,single_eigenvector,allowedq, T1, T2, a1m, a2m, b1,b2,spinor_set=triangle_initial_Densitymatrix(Int(NL),θ,Nq,gcutoff,uD,λ,enlarge_factor,V0_hBN,V1_hBN,ψ_hBN,V2_scalar,ϕ)

Area=Nq^2*√3/2*norm(a1m)^2


DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,eout,HartreeMatrix,FockMatrix=iteration_loop(initial_DensityMatrix,
                                                    allowedq,T1,T2,Nq,wave,single_Ham,single_MoirePo,constq,ϵr,overlapmatrix,filling,Area)


chern,Flink,chern_single,Flink_single,trace_condition,trace_condition_single,uniform,uniform_single=triangle_chern(Nq,wave,allowedq,HF_eigenvector,single_eigenvector,spinor_set)




scratch_dir = ENV["SCRATCH"]
savepath=joinpath(@__DIR__, "triangle_R5G_contact_interpolation_withhBN_skyrmionexcitation_v2/data_output$(Int(args[17]))/$(args[1])NL$(args[2])Nq$(args[3])theta$(args[4])constq$(args[5])ϵr$(args[6])uD$(args[7])filling$(args[8])cutoff$(args[9])lambda$(args[10])trytime$(args[11])enlarge$(args[12])V0_hBN$(args[13])V1_hBN$(args[14])ψ_hBN$(args[15])V2_scalar$(args[16])ϕ.jld2")
#savepath="test.jld2"


jldsave(savepath,single_Ham=single_Ham,single_MoirePo=single_MoirePo,
               single_eigenvalue=single_eigenvalue,spinor_set=spinor_set,
                chern=chern,Flink=Flink,chern_single=chern_single,
                Flink_single=Flink_single,arguments=args,
                densitymatrix=DIIS_input_DensityMatrix[1],energy=energy,eout=eout,
                TC=trace_condition,TCS=trace_condition_single,
                HFeigenvalue=HF_eigenvalue,
                uniform=uniform,uniform_single=uniform_single,
                HartreeMatrix=HartreeMatrix,FockMatrix=FockMatrix,HF_eigenvector=HF_eigenvector,
                single_eigenvector=single_eigenvector,allowedq=allowedq,T1=T1,T2=T2,wave=wave,
                a1m=a1m,a2m=a2m,b1=b1,b2=b2)


