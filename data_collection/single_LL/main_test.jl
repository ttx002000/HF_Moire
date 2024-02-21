using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames


include("../../src/construct_Vmatrix.jl")
include("../../src/single_particle.jl")
args=parse.(Int64,ARGS)

(eigenvalue_single, eigenvector_single, wave, T1, T2, g1T, g3T,allowedq,Nx,Ny,Nparticle,dimension,Lb)=energy_state(1.0);

Vmatrix=BuildVmatrix(Nx,Ny,wave,eigenvector_single,allowedq,T1,T2,g1T,g3T,dimension,Lb);

(MB_state_can, MB_state_integer)=Construct_MBstate(Nx,Ny,Nparticle,allowedq);
(reduced_Vcol,reduced_Vcoor)=reducedV(Vmatrix,Nx,Ny)



values=Construct_Manybodymatrix(reduced_Vcol,reduced_Vcoor,MB_state_can[args[1]],MB_state_integer[args[1]],eigenvalue_single)


jldsave(joinpath(@__DIR__, "data_output/LL_eigenvalue$(args[1])sector.jld2"),variebla1=values)



