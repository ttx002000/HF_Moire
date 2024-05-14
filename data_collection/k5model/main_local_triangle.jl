using JLD2
include("../../src/operators_k5model.jl")

V0=0.0
ϕ=0.0
Nq=3;
scale=1.0;
constq=0.5/Nq^2
Dfield=1.0



 wave, initial_DensityMatrix, BG_DensityMatrix,single_MoirePo, single_Ham, single_eigenvalue,allowedq, T1, T2, a1m, a2m=triangle_initial_Densitymatrix(Dfield,V0,ϕ,scale,Nq)
NoHFdensity=Densitymap(a1m,a2m,wave,initial_DensityMatrix)
 
for ja in 1:Nq^2
   A=randn(ComplexF64,2*length(wave),2*length(wave))
    initial_DensityMatrix[ja]=initial_DensityMatrix[ja]+(A+A')*0.001
end

DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue=iteration_loop(initial_DensityMatrix,BG_DensityMatrix,allowedq,T1,T2,Nq,wave,single_Ham,single_MoirePo,constq)
HFdensity=Densitymap(a1m,a2m,wave,DIIS_input_DensityMatrix[1])


chern,Flink,chern_single,Flink_single,trace_condition,trace_condition_single,uniform,uniform_single=triangle_chern(Nq,wave,scale,ϕ,Dfield,DIIS_input_DensityMatrix[1],constq)

energy=calculate_energy(Nq,wave,scale,ϕ,flux,DIIS_input_DensityMatrix[1],constq)
print(chern)

jldsave(joinpath(@__DIR__, "data_output/try3.jld2"),chern=chern,Flink=Flink,chern_single=chern_single,Flink_single=Flink_single,HFdensity=HFdensity,NoHFdensity=NoHFdensity,densitymatrix=DIIS_input_DensityMatrix,energy=energy,TC=trace_condition,TCS=trace_condition_single,HFeigenvalue=HF_eigenvalue)
