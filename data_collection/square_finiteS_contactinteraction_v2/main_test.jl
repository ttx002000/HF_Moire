using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames


include("../../src/operators_square_finiteS_contactinteraction_v2.jl")

args=parse.(Float64,ARGS)

#args=[0.5,0.0,0.0,3,1.0,1.0,1.0,1.0,1.0,1.0]

spin=args[1]
V0=args[2]
Nx=Int(args[3]);
Ny=Int(args[4]);
scale=args[5];
constq=args[6]/(Nx*Ny)
vf=args[7]
filling=Int(args[8])
gcutoff=args[9]
trytimes=Int(args[10])
filepos=Int(args[11])



 

overlapmatrix, wave, initial_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,single_eigenvector,allowedq, T1, T2, a1m, a2m=triangle_initial_Densitymatrix(spin,vf,V0,scale,Nx,Ny,filling,gcutoff)



DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,eout,HartreeMatrix,FockMatrix=iteration_loop(initial_DensityMatrix,
                                                    allowedq,T1,T2,Nx,Ny,wave,single_Ham,single_MoirePo,constq,overlapmatrix,filling)





scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "square_finiteS_contact_v2/data_output$(Int(args[11]))/$(args[1])spin$(args[2])V0$(args[3])Nx$(args[4])Ny$(args[5])scale$(args[6])constq$(args[7])vf$(args[8])filling$(args[9])cutoff$(args[10])trytime.jld2")
#savepath="test.jld2"


jldsave(savepath,
                arguments=args,
                densitymatrix=DIIS_input_DensityMatrix,energy=energy,eout=eout,
                HFeigenvalue=HF_eigenvalue,single_eigenvalue=single_eigenvalue,
                HartreeMatrix=HartreeMatrix,FockMatrix=FockMatrix,HF_eigenvector=HF_eigenvector,
                single_eigenvector=single_eigenvector,single_MoirePo=single_MoirePo, single_Ham=single_Ham)


