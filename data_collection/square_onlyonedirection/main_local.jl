using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames


include("../../src/operators_onlyonedirection.jl")

args=[2.0,2.0,0.0,3.0,1.0,1.5,1.0]
flux=args[1]*π
Vx=args[2]
Vy=args[3]
Nq=Int(args[4]);
scale=args[5];
constq=args[6]/Nq^2
trytimes=Int(args[7])



overlapmatrix, wave, initial_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,allowedq, T1, T2, a1m, a2m=square_initial_Densitymatrix(flux,Vx,Vy,scale,Nq)
NoHFdensity=Densitymap_kresolved(a1m,a2m,overlapmatrix,wave,initial_DensityMatrix)
 
for ja in 1:Nq^2
    A=randn(ComplexF64,length(wave),length(wave))
    initial_DensityMatrix[ja]=initial_DensityMatrix[ja]+(A+A')*0.001
end

DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue=iteration_loop(initial_DensityMatrix,allowedq,T1,T2,Nq,wave,single_Ham,single_MoirePo,constq,overlapmatrix)
HFdensity=Densitymap_kresolved(a1m,a2m,overlapmatrix,wave,DIIS_input_DensityMatrix[1])

chern,Flink,chern_single,Flink_single,trace_condition,trace_condition_single,uniform,uniform_single=square_chern(Nq,wave,scale,flux,Vx,Vy,DIIS_input_DensityMatrix[1],constq)
energy=calculate_energy_square(Nq,wave,scale,flux,Vx,Vy,DIIS_input_DensityMatrix[1],constq,overlapmatrix)

jldsave(joinpath(@__DIR__, "data_output/$(args[4])Nq$(args[1])flux$(args[2])Vx$(args[3])Vy$(args[5])scale$(args[6])constq$(args[7])try.jld2"),chern=chern,Flink=Flink,chern_single=chern_single,Flink_single=Flink_single,HFdensity=HFdensity,NoHFdensity=NoHFdensity,arguments=args,trace_condition=trace_condition,trace_condition_single=trace_condition_single,energy=energy,single_eigenvalue=single_eigenvalue)

