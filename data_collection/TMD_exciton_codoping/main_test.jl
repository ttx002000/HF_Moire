using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2
using Plots
using CSV,DataFrames

include(joinpath(@__DIR__,"../../src/operators_exciton_codoping.jl"))

args=parse.(Float64,ARGS)
println("This is the arguments$args")
#args=[0.36,-0.6,-0.6,0.0,0.0,20.8,107.7,20.8,-107.7,60000,70.0,3.89,-23.8,1.0,1.0,3.0]
parameters=args[1:13]
holenum=Int(args[14])
trytime=Int(args[15])
Nq=Int(args[16]);


#=
mt=parameters[1] actually args not parameters
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

holenum=parameters[14]
trytime=parameters[15]
Nq=parameters[16]
=#


wave, initial_DensityMatrix,BG_DensityMatrix, single_Ham, single_eigenvalue,allowedq, T1, T2, a1m, a2m,constq,single_chern=triangle_initial_Densitymatrix_control(parameters,Nq)

DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,HF_eigenvector,bound,energy,HF_chern=iteration_loop(initial_DensityMatrix,BG_DensityMatrix,allowedq,T1,T2,Nq,wave,single_Ham,constq,holenum)


xgrid,ygird,HFdensity=Densitymap(a1m,a2m,wave,DIIS_input_DensityMatrix[1])
savepath=joinpath(@__DIR__, "data_output/$(args[1])mt$(args[2])mm$(args[3])mb$(args[4])Vt$(args[5])phit$(args[6])Vm$(args[7])phim$(args[8])Vb$(args[9])phib$(args[10])er$(args[11])Eg$(args[12])theta$(args[13])w$(args[14])holenum$(args[15])trytime$(args[16])Nq.jld2")

jldsave(savepath,HFdensity=HFdensity,energy=energy,parameters=args,HF_eigenvalue=HF_eigenvalue,single_eigenvalue=single_eigenvalue,bound=bound,single_chern=single_chern,HF_chern=HF_chern)