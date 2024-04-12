using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames


include("../../src/operators_spinful.jl")

args=parse.(Float64,ARGS)
flux=args[1]*π
V0=args[2]
ϕ=args[3]/180*π
Nq=Int(args[4]);
scale=args[5];
constq=args[6]/Nq^2
ζ=args[7]
trytimes=Int(args[8])


overlapmatrix, wave, initial_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,allowedq, T1, T2, a1m, a2m=triangle_initial_Densitymatrix(flux,V0,ϕ,scale,Nq)
NoHFdensity=Densitymap(a1m,a2m,overlapmatrix,wave,initial_DensityMatrix)
 

for ja in 1:Nq^2, vi in 1:2
   A=randn(ComplexF64,length(wave),length(wave))
    initial_DensityMatrix[ja][vi]=(A+A')*0.001
end


DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue=iteration_loop(initial_DensityMatrix,allowedq,T1,T2,Nq,wave,single_Ham,single_MoirePo,constq,ζ,overlapmatrix)
HFdensity=Densitymap(a1m,a2m,overlapmatrix,wave,DIIS_input_DensityMatrix[1])

chern,Flink,chern_single,Flink_single,trace_condition,trace_condition_single=triangle_chern(Nq,wave,scale,ϕ,flux,DIIS_input_DensityMatrix[1],constq,ζ)
println(trace_condition_single)
energy=calculate_energy(Nq,wave,scale,ϕ,flux,DIIS_input_DensityMatrix[1],constq,ζ,overlapmatrix)
#jldsave(joinpath(@__DIR__, "data_output/$(args[4])Nq$(args[1])flux$(args[2])V0$(args[3])phi$(args[5])scale$(args[6])constq_$(args[7])try.jld2"),chern=chern,Flink=Flink,chern_single=chern_single,Flink_single=Flink_single,HFdensity=HFdensity,NoHFdensity=NoHFdensity,arguments=args,densitymatrix=DIIS_input_DensityMatrix,energy=energy,TC=trace_condition,TCS=trace_condition_single,HFeigenvalue=HF_eigenvalue,single_eigenvalue=single_eigenvalue,uniform=uniform,uniform_single=uniform_single)
jldsave(joinpath(@__DIR__, "data_output/$(args[4])Nq$(args[1])flux$(args[2])V0$(args[3])phi$(args[5])scale$(args[6])constq$(args[7])zeta_$(args[8])try.jld2"),chern=chern,Flink=Flink,chern_single=chern_single,Flink_single=Flink_single,HFdensity=HFdensity,NoHFdensity=NoHFdensity,arguments=args,energy=energy,TC=trace_condition,TCS=trace_condition_single,HFeigenvalue=HF_eigenvalue,single_eigenvalue=single_eigenvalue,uniform=uniform,uniform_single=uniform_single)

