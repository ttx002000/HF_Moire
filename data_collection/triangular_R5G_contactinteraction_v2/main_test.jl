using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames

#in this version, we use some arbitrary dispersion, and also I get rid of the Coulomb

include("../../src/operators_R5G_contactinteraction_v2.jl")

args=parse.(Float64,ARGS)



NL=args[1]
V0=args[2]
ϕ=args[3]/180*π
Nq=Int(args[4]);
scale=args[5];
constq=args[6]
filling=Int(args[7])
gcutoff=args[8]
trytimes=Int(args[9])
filepos=Int(args[10])






overlapmatrix, wave, initial_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,single_eigenvector,allowedq, T1, T2, a1m, a2m,spinor_set=triangle_initial_Densitymatrix(Int(NL),V0,ϕ,scale,Nq,gcutoff)

Area=Nq^2*√3/2*norm(a1m)^2


DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,eout,HartreeMatrix,FockMatrix=iteration_loop(initial_DensityMatrix,
                                                    allowedq,T1,T2,Nq,wave,single_Ham,single_MoirePo,constq,overlapmatrix,filling,Area)


chern,Flink,chern_single,Flink_single,trace_condition,trace_condition_single,uniform,uniform_single=triangle_chern(Nq,wave,scale,Int(NL),allowedq,HF_eigenvector,single_eigenvector,spinor_set)




scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "triangle_R5G_contact_v2/data_output$(Int(args[10]))/$(args[1])NL$(args[2])V0$(args[3])phi$(args[4])Nq$(args[5])scale$(args[6])constq$(args[7])filling$(args[8])cutoff$(args[9])trytime.jld2")



jldsave(savepath,single_Ham=single_Ham,single_MoirePo=single_MoirePo,
               single_eigenvalue=single_eigenvalue,
                chern=chern,Flink=Flink,chern_single=chern_single,
                Flink_single=Flink_single,arguments=args,
                densitymatrix=DIIS_input_DensityMatrix[1],energy=energy,eout=eout,
                TC=trace_condition,TCS=trace_condition_single,
                HFeigenvalue=HF_eigenvalue,
                uniform=uniform,uniform_single=uniform_single,
                HartreeMatrix=HartreeMatrix,FockMatrix=FockMatrix,HF_eigenvector=HF_eigenvector,
                single_eigenvector=single_eigenvector,allowedq=allowedq,T1=T1,T2=T2,wave=wave,spinor_set=spinor_set)


