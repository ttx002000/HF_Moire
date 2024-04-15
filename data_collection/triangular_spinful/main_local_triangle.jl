using JLD2
include("../../src/operators_spinful.jl")
flux=1.0*π
V0=0.0
ϕ=0.0
Nq=3;
scale=1.0;
constq=1.0/Nq^2
ζ=0.5



overlapmatrix, wave, initial_DensityMatrix, single_MoirePo, single_Ham,seed_MoirePo,single_eigenvalue,allowedq, T1, T2, a1m, a2m=triangle_initial_Densitymatrix(flux,V0,ϕ,scale,Nq)
NoHFdensity=Densitymap(a1m,a2m,overlapmatrix,wave,initial_DensityMatrix)
 
#=
for ja in 1:Nq^2, vi in 1:2
   A=randn(ComplexF64,length(wave),length(wave))
    initial_DensityMatrix[ja][vi]=initial_DensityMatrix[ja][vi]+(A+A')*0.1
end
=#

DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue=iteration_loop(initial_DensityMatrix,allowedq,T1,T2,Nq,wave,single_Ham,single_MoirePo,seed_MoirePo,constq,ζ,overlapmatrix)
HFdensity=Densitymap(a1m,a2m,overlapmatrix,wave,DIIS_input_DensityMatrix[1])

chern,Flink,chern_single,Flink_single,trace_condition,trace_condition_single,uniform,uniform_single=triangle_chern(Nq,wave,scale,ϕ,flux,DIIS_input_DensityMatrix[1],constq,ζ)

energy=calculate_energy(Nq,wave,scale,ϕ,flux,DIIS_input_DensityMatrix[1],constq,ζ,overlapmatrix)
println(trace_condition_single)
jldsave(joinpath(@__DIR__, "data_output/trywithpining.jld2"),HFeigenvalue=HF_eigenvalue,chern=chern,Flink=Flink,chern_single=chern_single,Flink_single=Flink_single,HFdensity=HFdensity,NoHFdensity=NoHFdensity,densitymatrix=DIIS_input_DensityMatrix,energy=energy,TC=trace_condition,TCS=trace_condition_single)
