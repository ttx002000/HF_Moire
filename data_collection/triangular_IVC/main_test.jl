using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames


include(joinpath(@__DIR__,"../../src/operators_IVC.jl"))

args=parse.(Float64,ARGS)
flux=args[1]*π
V0=args[2]
ϕ=args[3]/180*π
Nq=Int(args[4]);
scale=args[5];
constq=args[6]/Nq^2
ζ=args[7]
trytimes=Int(args[8])
shiftindex=Int(args[9])

#=
    flux=arg[1]
    V0=arg[2]
    ϕ=arg[3]
    Nq=arg[4];
    scale=arg[5];
    constq=arg[6];
    ζ=arg[7];
    trytimes=arg[8]
    shiftindex=arg[9]
=#


overlapmatrix, wave, initial_DensityMatrix, single_MoirePo, single_Ham,single_eigenvalue,allowedq, T1, T2, a1m, a2m=triangle_initial_Densitymatrix(flux,V0,ϕ,scale,Nq,shift_index)
NoHFdensity=Densitymap(a1m,a2m,overlapmatrix,wave,initial_DensityMatrix)
 

for ja in 1:Nq^2
   A=randn(ComplexF64,2*length(wave),2*length(wave))
    initial_DensityMatrix[ja]=initial_DensityMatrix[ja]+(A+A')*0.1
end


DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue=iteration_loop(initial_DensityMatrix,allowedq,T1,T2,Nq,wave,single_Ham,single_MoirePo,constq,ζ,overlapmatrix)
HFdensity=Densitymap(a1m,a2m,overlapmatrix,wave,DIIS_input_DensityMatrix[1])

chern,Flink,chern_single,Flink_single,trace_condition,trace_condition_single,uniform,uniform_single=triangle_chern(Nq,wave,scale,ϕ,flux,DIIS_input_DensityMatrix[1],constq,ζ,allowedq[shift_index])

energy=calculate_energy(Nq,wave,scale,ϕ,flux,DIIS_input_DensityMatrix[1],constq,ζ,overlapmatrix,allowedq[shift_index])
jldsave(joinpath(@__DIR__, "data_output/$(args[1])flux$(args[2])V0$(args[3])phi$(args[4])Nq$(args[5])scale$(args[6])constq$(args[7])zeta$(args[8])try$(args[9])shift.jld2"),chern=chern,Flink=Flink,chern_single=chern_single,Flink_single=Flink_single,HFdensity=HFdensity,NoHFdensity=NoHFdensity,arguments=args,energy=energy,TC=trace_condition,TCS=trace_condition_single,HFeigenvalue=HF_eigenvalue,single_eigenvalue=single_eigenvalue,uniform=uniform,uniform_single=uniform_single,DM=DIIS_input_DensityMatrix[1])

