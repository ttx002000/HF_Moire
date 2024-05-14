using JLD2
include("../../src/operators_k5model.jl")

V0=0.1
ϕ=0.0
Nq=3;
period=10
scale=4*π/(period*√3);
ϵr=7
constq=10447.22667/(ϵr*Nq^2*period^2)
Dfield=10.0



 wave, initial_DensityMatrix, BG_DensityMatrix,single_MoirePo, single_Ham, single_eigenvalue,allowedq, T1, T2, a1m, a2m=triangle_initial_Densitymatrix(Dfield,V0,ϕ,scale,Nq)
NoHFdensity=Densitymap(a1m,a2m,wave,initial_DensityMatrix)
 

for ja in 1:Nq^2
   A=randn(ComplexF64,2*length(wave),2*length(wave))
    initial_DensityMatrix[ja]=initial_DensityMatrix[ja]+(A+A')*0.5
end

DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,bound=iteration_loop(initial_DensityMatrix,BG_DensityMatrix,allowedq,T1,T2,Nq,wave,single_Ham,single_MoirePo,constq)
HFdensity=Densitymap(a1m,a2m,wave,DIIS_input_DensityMatrix[1])


chern,Flink,chern_single,Flink_single,trace_condition,trace_condition_single,uniform,uniform_single=triangle_chern(Nq,wave,scale,ϕ,Dfield,DIIS_input_DensityMatrix[1],constq)

energy=calculate_energy(Nq,wave,scale,ϕ,Dfield,DIIS_input_DensityMatrix[1],constq)
print(chern)

jldsave(joinpath(@__DIR__, "data_output/try3.jld2"),chern=chern,Flink=Flink,chern_single=chern_single,Flink_single=Flink_single,HFdensity=HFdensity,NoHFdensity=NoHFdensity,densitymatrix=DIIS_input_DensityMatrix,energy=energy,TC=trace_condition,TCS=trace_condition_single,HFeigenvalue=HF_eigenvalue)
