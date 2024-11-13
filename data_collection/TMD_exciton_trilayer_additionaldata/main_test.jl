using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2
using Plots
using CSV,DataFrames

include(joinpath(@__DIR__,"../../src/operators_exciton_additionaldata.jl"))

args=parse.(Float64,ARGS)
println("This is the arguments$args")




parameters=args[1:4]
holenum=Int(args[5])
trytime=Int(args[6])
Nq=Int(args[7]);
seednum=Int(args[8])

#=
ϵr=parameters[10]
Eg=parameters[11]
θ=parameters[12]/180*π
w=parameters[13]

holenum=parameters[14]
trytime=parameters[15]
Nq=parameters[16]
=#


wave, initial_DensityMatrix,BG_DensityMatrix, single_Ham, single_eigenvalue,allowedq, T1, T2, a1m, a2m,constq=triangle_initial_Densitymatrix_control(parameters,Nq,seednum)

DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,bound=iteration_loop_control(initial_DensityMatrix,BG_DensityMatrix,allowedq,T1,T2,Nq,wave,single_Ham,constq,holenum,seednum)
xgrid,ygird,HFdensity=Densitymap(a1m,a2m,wave,DIIS_input_DensityMatrix[1])
energy,HF_vectors=calculate_energy(Nq,wave,DIIS_input_DensityMatrix[1],constq,T1,T2,allowedq,single_Ham)


savepath=joinpath(@__DIR__, "data_output/$(args[1])er$(args[2])Eg$(args[3])theta$(args[4])w$(args[5])holenum$(args[6])trytime$(args[7])Nq$(args[8])seed.jld2")

jldsave(savepath,HFdensity=HFdensity,HF_vectors=HF_vectors,energy=energy,parameters=args,HF_eigenvalue=HF_eigenvalue,single_eigenvalue=single_eigenvalue,bound=bound)