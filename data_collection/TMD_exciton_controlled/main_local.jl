using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames


include("../../src/operators_exciton.jl")

#args=parse.(Float64,ARGS)


parameters=[0.35,0.4,0.35,8,270.0,14,190.0,10,270.0,8,300.0,2.0,0.0]
Nq=3;
holenum=3
seednum=2

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
DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,bound=iteration_loop(initial_DensityMatrix,BG_DensityMatrix,allowedq,T1,T2,Nq,wave,single_Ham,constq,holenum)
xgrid,ygird,HFdensity=Densitymap(a1m,a2m,wave,DIIS_input_DensityMatrix[1])
energy=calculate_energy(Nq,wave,DIIS_input_DensityMatrix[1],constq,T1,T2,allowedq,single_Ham)



savepath=joinpath(@__DIR__, "data_output/try.jld2")

jldsave(savepath,HFdensity=HFdensity)



