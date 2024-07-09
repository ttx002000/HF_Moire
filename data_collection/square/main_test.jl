using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames


include("../../src/operators.jl")

args=parse.(Float64,ARGS)
#args=[2.0,0.0,0.0,3.0,1.0,1.0,1.0]
flux=args[1]*π
V0=args[2]
ϕ=args[3]/180*π
Nq=Int(args[4]);
scale=args[5];
constq=args[6]/Nq^2
trytimes=Int(args[7])



overlapmatrix, wave, initial_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,allowedq, T1, T2, a1m, a2m=square_initial_Densitymatrix(flux,V0,ϕ,scale,Nq)
NoHFdensity=Densitymap(a1m,a2m,overlapmatrix,wave,initial_DensityMatrix)
 
for ja in 1:Nq^2
    A=randn(ComplexF64,length(wave),length(wave))
    initial_DensityMatrix[ja]=initial_DensityMatrix[ja]+(A+A')*0.01
end

DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,HF_eigenvec=iteration_loop(initial_DensityMatrix,allowedq,T1,T2,Nq,wave,single_Ham,single_MoirePo,constq,overlapmatrix)
HFdensity=Densitymap(a1m,a2m,overlapmatrix,wave,DIIS_input_DensityMatrix[1])
HFdensity_conjugate=Densitymap(a1m,a2m,overlapmatrix,wave,conj.(DIIS_input_DensityMatrix[1]))
chern,Flink,chern_single,Flink_single,trace_condition,trace_condition_single,uniform,uniform_single,chern_conjugate,Flink_conjugate,trace_condition_conjugate,uniform_conjugate=square_chern(Nq,wave,scale,ϕ,flux,DIIS_input_DensityMatrix[1],constq)


jldsave(joinpath(@__DIR__, "data_output/$(args[4])Nq$(args[1])flux$(args[2])V0$(args[3])phi$(args[5])scale$(args[6])constq$(args[7])try.jld2"),chern=chern,Flink=Flink,uniform=uniform,trace_condition=trace_condition,chern_single=chern_single,Flink_single=Flink_single,trace_condition_single=trace_condition_single,uniform_single=uniform_single,HFdensity=HFdensity,HFdensity_conjugate=HFdensity_conjugate,NoHFdensity=NoHFdensity,arguments=args,chern_conjugate=chern_conjugate,Flink_conjugate=Flink_conjugate,uniform_conjugate=uniform_conjugate,trace_condition_conjugate=trace_condition_conjugate,HF_eigenvec=HF_eigenvec)


