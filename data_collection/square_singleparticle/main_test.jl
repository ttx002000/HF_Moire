using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames


include("../../src/operators_singleparticle.jl")

args=parse.(Float64,ARGS)
#args=[2.0,0.0,0.0,3.0,1.0,1.0,1.0]
#args=[2.5,5.0,0.0,12.0,1.0,9.01]
flux=args[1]*π
V0=args[2]
ϕ=args[3]/180*π
Nq=Int(args[4]);
scale=args[5];
maxg=args[6]


β=flux/scale^2
V0=V0*exp(β/4*scale^2)

wave, initial_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,allowedq, T1, T2, a1m, a2m=square_initial_Densitymatrix(flux,V0,ϕ,scale,Nq,maxg)
chern_single,Flink_single,trace_condition_single,uniform_single=square_chern(Nq,wave,scale,ϕ,flux)


jldsave(joinpath(@__DIR__, "data_output/NormalizedV0_$(args[4])Nq$(args[1])flux$(args[2])V0$(args[3])phi$(args[5])scale$(args[6])maxg.jld2"),chern_single=chern_single,Flink_single=Flink_single,trace_condition_single=trace_condition_single,uniform_single=uniform_single,arguments=args, single_eigenvalue=single_eigenvalue)


