using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames


include("../../src/operators_square_finiteS_contactinteraction.jl")

args=parse.(Float64,ARGS)

#args=[0.5,0.0,0.0,3,1.0,1.0,1.0,1.0,1.0,1.0]

spin=args[1]
V0=args[2]
Nq=Int(args[3]);
scale=args[4];
constq=args[5]/Nq^2
vf=args[6]
filling=Int(args[7])
gcutoff=args[8]
trytimes=Int(args[9])
filepos=Int(args[10])






overlapmatrix, wave, initial_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,single_eigenvector,allowedq, T1, T2, a1m, a2m=triangle_initial_Densitymatrix(spin,vf,V0,ϕ,scale,Nq,filling,gcutoff)

 
for ja in 1:Nq^2
    A=randn(ComplexF64,length(wave),length(wave))
    initial_DensityMatrix[ja]=(A+A')*0.2
end

DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,eout,HartreeMatrix,FockMatrix=iteration_loop(initial_DensityMatrix,
                                                    allowedq,T1,T2,Nq,wave,single_Ham,single_MoirePo,constq,overlapmatrix,filling)


chern,Flink,chern_single,Flink_single,trace_condition,trace_condition_single,uniform,uniform_single=triangle_chern(Nq,wave,scale,spin,vf,allowedq,HF_eigenvector,single_eigenvector)




scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "square_finiteS_contact/data_output$(Int(args[10]))/$(args[1])spin$(args[2])V0$(args[3])Nq$(args[4])scale$(args[5])constq$(args[6])vf$(args[7])filling$(args[8])cutoff$(args[9])trytime.jld2")
#savepath="test.jld2"


jldsave(savepath,
                chern=chern,Flink=Flink,chern_single=chern_single,
                Flink_single=Flink_single,arguments=args,
                densitymatrix=DIIS_input_DensityMatrix,energy=energy,eout=eout,
                TC=trace_condition,TCS=trace_condition_single,
                HFeigenvalue=HF_eigenvalue,single_eigenvalue=single_eigenvalue,
                uniform=uniform,uniform_single=uniform_single,
                HartreeMatrix=HartreeMatrix,FockMatrix=FockMatrix,HF_eigenvector=HF_eigenvector,
                single_eigenvector=single_eigenvector)


