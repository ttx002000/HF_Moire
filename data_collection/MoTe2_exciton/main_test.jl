using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2
using Plots
using CSV,DataFrames

include(joinpath(@__DIR__,"../../src/operators_exciton_MoTe2.jl"))

args=parse.(Float64,ARGS)
println("This is the arguments$args")



#args=[-1.0,0.62,0.62,0.0,0.0,-11.2,91.0,-11.2,-91.0,5.0,100.0,10.0,-13.3,0.0,1.0,1.0,4.01]

parameters=args[1:13]
elecnum=Int(args[14])
trytime=Int(args[15])
geonum=Int(args[16]);
cutoffnum=args[17]




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
period=parameters[12]
w=parameters[13]

elecnum=parameters[14]
trytime=parameters[15]
geonum=parameters[16]
cutoffnum=parameters[17]
=#

wave, initial_DensityMatrix,BG_DensityMatrix, single_Ham, single_eigenvalue,single_eigenvector,allowedq, T1, T2, a1m, a2m,constq=triangle_initial_Densitymatrix(parameters,geonum,cutoffnum)

DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,HF_vector_tosave,bound,energy=iteration_loop(initial_DensityMatrix,BG_DensityMatrix,allowedq,T1,T2,geonum,wave,single_Ham,constq,elecnum)

   

xgrid,ygird,HFdensity=Densitymap(a1m,a2m,wave,DIIS_input_DensityMatrix[1])



savepath=joinpath(@__DIR__, "data_output/$(args[1])mt$(args[2])mm$(args[3])mb$(args[4])Vt$(args[5])phit$(args[6])Vm$(args[7])phim$(args[8])Vb$(args[9])phib$(args[10])er$(args[11])Eg$(args[12])am$(args[13])w$(args[14])elecnum$(args[15])try$(args[16])geo$(args[17])cut.jld2")




jldsave(savepath,HFdensity=HFdensity,HF_vector=HF_vector_tosave,energy=energy,parameters=args,HF_eigenvalue=HF_eigenvalue,single_eigenvalue=single_eigenvalue,bound=bound,a1m=a1m,a2m=a2m,dimension=length(wave),constq=constq)