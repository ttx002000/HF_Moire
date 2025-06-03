using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames


include("../../src/operators_R5G_contactinteraction.jl")

args=parse.(Float64,ARGS)

#args=[0.5,0.0,0.0,3,1.0,1.0,1.0,1.0,1.0,1.0]

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






overlapmatrix, wave, initial_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,single_eigenvector,allowedq, T1, T2, a1m, a2m=triangle_initial_Densitymatrix(Int(NL),V0,ϕ,scale,Nq,gcutoff)

Area=Nq^2*√3/2*norm(a1m)^2


DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,eout,HartreeMatrix,FockMatrix=iteration_loop(initial_DensityMatrix,
                                                    allowedq,T1,T2,Nq,wave,single_Ham,single_MoirePo,constq,overlapmatrix,filling,Area)


chern,Flink,chern_single,Flink_single,trace_condition,trace_condition_single,uniform,uniform_single=triangle_chern(Nq,wave,scale,Int(NL),allowedq,HF_eigenvector,single_eigenvector)




scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "triangle_R5G_contact/data_output$(Int(args[10]))/$(args[1])NL$(args[2])V0$(args[3])phi$(args[4])Nq$(args[5])scale$(args[6])constq$(args[7])filling$(args[8])cutoff$(args[9])trytime.jld2")
#savepath="test.jld2"


jldsave(savepath,
                chern=chern,Flink=Flink,chern_single=chern_single,
                Flink_single=Flink_single,arguments=args,
                densitymatrix=DIIS_input_DensityMatrix,energy=energy,eout=eout,
                TC=trace_condition,TCS=trace_condition_single,
                HFeigenvalue=HF_eigenvalue,single_eigenvalue=single_eigenvalue,
                uniform=uniform,uniform_single=uniform_single,
                HartreeMatrix=HartreeMatrix,FockMatrix=FockMatrix,HF_eigenvector=HF_eigenvector,
                single_eigenvector=single_eigenvector,allowedq=allowedq,T1=T1,T2=T2,wave=wave)


