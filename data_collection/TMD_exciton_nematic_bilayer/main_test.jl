using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2
using Plots
using CSV,DataFrames

include(joinpath(@__DIR__,"../../src/operators_nematic_bilayer.jl"))

args=parse.(Float64,ARGS)
println("This is the arguments$args")
parameters=args[1:11]
holenum=Int(args[12])
trytime=Int(args[13])
Nq=Int(args[14]);
seednum=Int(args[15])

#=
mt=parameters[1] actually args not parameters
mb=parameters[2]
Vt=parameters[3]
ϕt=parameters[4]/180*π
Vb=parameters[5]
ϕb=parameters[6]/180*π
ϵr=parameters[7]
Eg=parameters[8]
θ=parameters[9]/180*π
w=parameters[10]
omega=parameters[11]

holenum=parameters[12]
trytime=parameters[13]
Nq=parameters[14]
seed=parametesd[15]
=#


wave, initial_DensityMatrix,BG_DensityMatrix, single_Ham, single_eigenvalue,allowedq, T1, T2, a1m, a2m,constq=triangle_initial_Densitymatrix_control(parameters,Nq,seednum)

DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,bound=iteration_loop_control(initial_DensityMatrix,BG_DensityMatrix,allowedq,T1,T2,Nq,wave,single_Ham,constq,holenum,seednum)
xgrid,ygird,HFdensity=Densitymap(a1m,a2m,wave,DIIS_input_DensityMatrix[1])
energy,HF_vectors=calculate_energy(Nq,wave,DIIS_input_DensityMatrix[1],constq,T1,T2,allowedq,single_Ham)


savepath=joinpath(@__DIR__, "data_output/$(args[1])mt$(args[2])mb$(args[3])Vt$(args[4])phit$(args[5])Vb$(args[6])phib$(args[7])er$(args[8])Eg$(args[9])theta$(args[10])w$(args[11])omega$(args[12])holenum$(args[13])trytime$(args[14])Nq$(args[15])seed.jld2")

jldsave(savepath,HF_vectors=HF_vectors,HFdensity=HFdensity,energy=energy,parameters=args,HF_eigenvalue=HF_eigenvalue,single_eigenvalue=single_eigenvalue,bound=bound)