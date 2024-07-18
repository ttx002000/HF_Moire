using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames


include("../../src/operators_DIISon_physicalunits.jl")

args=parse.(Float64,ARGS)


#args=[1.0,0.0,0.0,3.0,20*π^2/√3,1.0]




flux=args[1]*π
V0=args[2]
ϕ=args[3]/180*π
Nq=Int(args[4]);
Vc=args[5]
trytimes=Int(args[6])


β=4*π/(√3)
scale=sqrt(flux/π)
am=4π/(scale*√3)
Auc=√3/2*am^2
constq=Vc/(Nq^2*Auc)


overlapmatrix, wave, initial_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,allowedq, T1, T2, a1m, a2m=triangle_initial_Densitymatrix(V0,ϕ,scale,Nq)
NoHFdensity=Densitymap(a1m,a2m,overlapmatrix,wave,initial_DensityMatrix)
 
for ja in 1:Nq^2
    A=randn(ComplexF64,length(wave),length(wave))
    initial_DensityMatrix[ja]=initial_DensityMatrix[ja]+(A+A')*0.1
end

DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,HF_eigenvec=iteration_loop(initial_DensityMatrix,allowedq,T1,T2,Nq,wave,single_Ham,single_MoirePo,constq,overlapmatrix)
HFdensity=Densitymap(a1m,a2m,overlapmatrix,wave,DIIS_input_DensityMatrix[1])


chern,Flink,chern_single,Flink_single,trace_condition,trace_condition_single,uniform,uniform_single=triangle_chern(Nq,wave,scale,ϕ,DIIS_input_DensityMatrix[1],constq)
kin_energy,Ha_energy,Fk_energy=calculate_energy_HFresolved(Nq,wave,scale,ϕ,DIIS_input_DensityMatrix[1],constq,overlapmatrix)

jldsave(joinpath(@__DIR__, "data_output/$(args[4])Nq$(args[1])flux$(args[2])V0$(args[3])phi$(args[5])Vc_$(args[6])try.jld2"),chern=chern,Flink=Flink,chern_single=chern_single,Flink_single=Flink_single,HFdensity=HFdensity,NoHFdensity=NoHFdensity,arguments=args,kin_energy=kin_energy,Ha_energy=Ha_energy,Fk_energy=Fk_energy,TC=trace_condition,TCS=trace_condition_single,HFeigenvalue=HF_eigenvalue,single_eigenvalue=single_eigenvalue,uniform=uniform,uniform_single=uniform_single)

