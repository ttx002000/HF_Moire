using JLD2
include("../../src/operators_k5model.jl")

args=parse.(Float64,ARGS)
V0=args[1]
ϕ=args[2]/180*π
Nq=Int(args[3]);
period=args[4]
ϵr=args[5]
Dfield=args[6]
trytime=args[7]

scale=4*π/(period*√3);
constq=10447.22667/(ϵr*Nq^2*period^2)

 wave, initial_DensityMatrix, BG_DensityMatrix,single_MoirePo, single_Ham, single_eigenvalue,allowedq, T1, T2, a1m, a2m=triangle_initial_Densitymatrix(Dfield,V0,ϕ,scale,Nq)
NoHFdensity=Densitymap(a1m,a2m,wave,initial_DensityMatrix)
 

for ja in 1:Nq^2
   A=randn(ComplexF64,2*length(wave),2*length(wave))
    initial_DensityMatrix[ja]=initial_DensityMatrix[ja]+(A+A')*0.5
end

DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,bound,eout=iteration_loop(initial_DensityMatrix,BG_DensityMatrix,allowedq,T1,T2,Nq,wave,single_Ham,single_MoirePo,constq)
HFdensity=Densitymap(a1m,a2m,wave,DIIS_input_DensityMatrix[1])


chern,Flink,chern_single,Flink_single,trace_condition,trace_condition_single,uniform,uniform_single=triangle_chern(Nq,wave,scale,ϕ,Dfield,DIIS_input_DensityMatrix[1],constq)

energy=calculate_energy(Nq,wave,scale,ϕ,Dfield,DIIS_input_DensityMatrix[1],constq)
print(chern)

jldsave(joinpath(@__DIR__, "data_output/$(args[1])V0$(args[2])phi$(args[3])Nq$(args[4])period$(args[5])er$(args[6])D$(args[7])try.jld2"),chern=chern,Flink=Flink,chern_single=chern_single,Flink_single=Flink_single,HFdensity=HFdensity,NoHFdensity=NoHFdensity,energy=energy,TC=trace_condition,TCS=trace_condition_single,HFeigenvalue=HF_eigenvalue,single_eigenvalue=single_eigenvalue,eout=eout)
