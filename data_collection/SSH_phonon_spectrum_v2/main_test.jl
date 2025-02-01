
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2
using Plots
using CSV,DataFrames

include(joinpath(@__DIR__,"../../src/operators_SSH_phonon_v2.jl"))

args=parse.(Float64,ARGS)
#args=[20,20,-0.37,2.0,0.08,-1.0,-1.0,1.0,1.0,1.25,0.01,0.1,20,20,1]
Nx=Int(args[1])
Ny=Int(args[2])
tper=args[3]
tpa=args[4]
tNNN=args[5]
α=args[6]
β=args[7]
K=args[8]
KNNN=args[9]
filling=(args[10])
temp=args[11]
shearstrength=args[12]
Nqx=Int(args[13])
Nqy=Int(args[14])
filepos=Int(args[15])



gmatrix,Keff_set, Kbare_momentum_set,qset, qindex=runrunrun(α,β,Nx,Ny,
                                              Nqx,Nqy,temp,tper,tpa,tNNN,
                                               K,KNNN,shearstrength,filling)



scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "SSH_phonon_v2/data_output$(filepos)/$(Int(args[1]))Nx$(Int(args[2]))Ny$(args[3])tper$(args[4])tpa$(args[5])tNNN$(args[6])alpha$(args[7])beta$(args[8])K$(args[9])KNNN$(args[10])filling$(args[11])temp$(args[12])shear$(args[13])Nqx$(args[14])Nqy.jld2")




jldsave(savepath,gmatrix=gmatrix,Keff_set=Keff_set,Kbare_momentum_set=Kbare_momentum_set,qset=qset, qindex=qindex)
