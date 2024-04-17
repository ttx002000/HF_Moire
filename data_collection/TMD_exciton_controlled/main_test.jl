


include("../../src/operators_exciton.jl")

args=parse.(Float64,ARGS)
println("This is the arguments$args")
parameters=args[1:13]
#parameters=[0.35,0.4,0.35,-10.0,70.0,10.0,80.0,10.0,1.0,10,100.0,2.0,5.0]
holenum=Int(args[14])
trytime=Int(args[15])
Nq=Int(args[16]);
seednum=Int(args[17])

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


wave, initial_DensityMatrix,BG_DensityMatrix, single_Ham, single_eigenvalue,allowedq, T1, T2, a1m, a2m,constq=triangle_initial_Densitymatrix_control(parameters,Nq,seednum)

DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,bound=iteration_loop_control(initial_DensityMatrix,BG_DensityMatrix,allowedq,T1,T2,Nq,wave,single_Ham,constq,holenum,seednum)
xgrid,ygird,HFdensity=Densitymap(a1m,a2m,wave,DIIS_input_DensityMatrix[1])
energy=calculate_energy(Nq,wave,DIIS_input_DensityMatrix[1],constq,T1,T2,allowedq,single_Ham)


savepath=joinpath(@__DIR__, "data_output/$(args[1])mt$(args[2])mm$(args[3])mb$(args[4])Vt$(args[5])phit$(args[6])Vm$(args[7])phim$(args[8])Vb$(args[9])phib$(args[10])er$(args[11])Eg$(args[12])theta$(args[13])w$(args[14])holenum$(args[15])trytime$(args[16])Nq$(args[17])seed.jld2")

jldsave(savepath,HFdensity=HFdensity,energy=energy,parameters=args,HF_eigenvalue=HF_eigenvalue,single_eigenvalue=single_eigenvalue,bound=bound)