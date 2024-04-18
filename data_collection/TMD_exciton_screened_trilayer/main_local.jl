using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames


include("../../src/operators_exciton.jl")

#args=parse.(Float64,ARGS)


parameters=[0.35,0.4,0.35,8,270.0,14,190.0,10,80.0,5,220.0,2.0,0.0]
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
xgrid,ygird,HFdensity=Densitymap(a1m,a2m,wave,initial_DensityMatrix)




savepath=joinpath(@__DIR__, "data_output/try.jld2")

jldsave(savepath,HFdensity=HFdensity)



