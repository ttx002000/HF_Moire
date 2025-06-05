
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2



include(joinpath(@__DIR__,"../../src/operators_SSH_mirror.jl"))

args=parse.(Float64,ARGS)
Nx=Int(args[1])
Ny=Int(args[2])
tper=args[3]
tpa=args[4]
tNNN=args[5]

α=args[6]
β=args[7]


K=args[8]
gshear=args[9]

filling=(args[10])
KNNN=args[11]
temp=args[12]
stop_standard=args[13]
filepos=Int(args[14])

Nelec=Nx*Ny*filling

scratch_dir = ENV["SCRATCH"]
seedpath=joinpath(scratch_dir, "SSH_v3/test$(Int(args[14]))/mirror_seed/$(Int(args[1]))Nx$(Int(args[2]))Ny$(args[3])tper$(args[4])tpa$(args[5])tNNN$(args[6])alpha$(args[7])beta$(args[8])K$(args[9])gshear$(args[10])filling$(args[11])KNNN$(args[12])temp$(args[13])stop.jld2")
st1=load(seedpath)
unsym_phonon_coor=st1["phonon_coor"]
unsym_free_energy=st1["free_energy"]

sym_phonon_coor,sym_free_energy=construct_sym(unsym_phonon_coor,Nx,Ny,tper,tpa,tNNN,α,β,K,KNNN,gshear,temp,Nelec)






savepath=joinpath(scratch_dir, "SSH_v3/test$(Int(args[14]))/mirror_result/$(Int(args[1]))Nx$(Int(args[2]))Ny$(args[3])tper$(args[4])tpa$(args[5])tNNN$(args[6])alpha$(args[7])beta$(args[8])K$(args[9])gshear$(args[10])filling$(args[11])KNNN$(args[12])temp$(args[13])stop.jld2")

jldsave(savepath,sym_phonon_coor=sym_phonon_coor,sym_free_energy=sym_free_energy,unsym_phonon_coor=unsym_phonon_coor,unsym_free_energy=unsym_free_energy)
         
