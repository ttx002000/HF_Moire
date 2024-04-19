using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2
using Plots
using CSV,DataFrames

include(joinpath(@__DIR__,"../../src/operators_exciton_screened_bilayer.jl"))

#args=parse.(Float64,ARGS)


parameters=[0.3,0.54,10,45.8366,-10,177.6169,5,80.0,2.0,0.0]
Nq=3;
holenum=3
seednum=1

#=
mt=parameters[1]
mm=parameters[2]
mb=parameters[3]
Vt=parameters[4]
ϕt=parameters[5]/180*π
Vm=parameters[6]
ϕm=parameters[7]/180*π
Vb=parameters[8]
ϕb=parameters[9]/180*π
ϵr=parameters[10]
Eg=parameters[11]
θ=parameters[12]/180*π
w=parameters[13]
=#




wave, initial_DensityMatrix,BG_DensityMatrix, single_Ham, single_eigenvalue,allowedq, T1, T2, a1m, a2m,constq=triangle_initial_Densitymatrix_control(parameters,Nq,seednum)

DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,bound=iteration_loop_control(initial_DensityMatrix,BG_DensityMatrix,allowedq,T1,T2,Nq,wave,single_Ham,constq,holenum,seednum)
xgrid,ygird,HFdensity=Densitymap(a1m,a2m,wave,DIIS_input_DensityMatrix[1])
energy=calculate_energy(Nq,wave,DIIS_input_DensityMatrix[1],constq,T1,T2,allowedq,single_Ham)



savepath=joinpath(@__DIR__, "data_output/try.jld2")

jldsave(savepath,HFdensity=HFdensity)



