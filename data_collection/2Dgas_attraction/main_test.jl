using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2


# The structure is basically the same as the skyrmion excitation v4, but this is for 2D gas
include("../../src/2Dgas_attraction.jl")

args=parse.(Float64,ARGS)
#args=[1,3.51,0.0,10.0,1.0,1.0,1.0,1.0,1.0]

scale=args[1]
gcutoff=args[2]
constq=args[3]
gatedis=args[4]
filling=Int(args[5])
lpo=args[6]
attstr=args[7]
trytime=Int(args[8])
filepos=Int(args[9])




wave, initial_DensityMatrix,  single_Ham, T1, T2, a1m, a2m, b1,b2,Area=triangle_initial_Densitymatrix(scale,gcutoff)



DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,eout,HartreeMatrix,FockMatrix=iteration_loop(initial_DensityMatrix,
                                                                            T1,T2,wave,single_Ham,
                                                                            constq,gatedis,
                                                                            filling,Area,lpo,attstr)





scratch_dir = ENV["SCRATCH"]

savepath=joinpath(scratch_dir, "2Dgas_attraction/data_output$(Int(args[9]))/$(args[1])scale$(args[2])gcut$(args[3])constq$(args[4])gatedis$(args[5])filling$(args[6])lpo$(args[7])attstr$(args[8])trytime.jld2")
#savepath="test.jld2"


jldsave(savepath,single_Ham=single_Ham,
                arguments=args,
                densitymatrix=DIIS_input_DensityMatrix[1],energy=energy,eout=eout,
                HFeigenvalue=HF_eigenvalue,
                HartreeMatrix=HartreeMatrix,FockMatrix=FockMatrix,HF_eigenvector=HF_eigenvector,
                 T1=T1,T2=T2,wave=wave,
                a1m=a1m,a2m=a2m,b1=b1,b2=b2)


