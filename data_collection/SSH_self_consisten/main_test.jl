
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2
using Plots
using CSV,DataFrames

include(joinpath(@__DIR__,"../../src/operators_SSH.jl"))

args=parse.(Float64,ARGS)
#args=[10,10,0.37,2.0,0.16,-1.2,-1.2,1.0,1.0,1.25,1.0]
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
trytime=Int(args[11])
Nelec=Int(round(Nx*Ny*filling))

H0, orbital_id, phonon_id, px_xbond, px_ybond, py_xbond, py_ybond, NNN_sp_d1, NNN_sp_d2=initialize(Nx,Ny,tper,tpa,tNNN)


phonon_coor, Hph, grad_old, E_old, Eelec_new, Egap,ave_npa,FL=iteration(Nx,Ny,Nelec,px_xbond,px_ybond,py_xbond,py_ybond,NNN_sp_d1,NNN_sp_d2,orbital_id,phonon_id,α,β,K,KNNN,H0)
dis_x,dis_y,max_record,average_record=resh_phonon(phonon_coor,phonon_id,Nx,Ny)


FFF=eigen(H0+Hph)
spectrum=FFF.values
savepath=joinpath(@__DIR__, "data_output/$(Int(args[1]))Nx$(Int(args[2]))Ny$(args[3])tper$(args[4])tpa$(args[5])tNNN$(args[6])alpha$(args[7])beta$(args[8])K$(args[9])KNNN$(args[10])filling$(Int(args[11]))try.jld2")
jldsave(savepath,phonon_id=phonon_id,orbital_id=orbital_id,phonon_coor=phonon_coor,E_total=E_old,Eelec=Eelec_new,Echemical=E_old-FL*ave_npa/(Nx*Ny),spectrum=spectrum,Egap=Egap,ave_npa=ave_npa,max_record=max_record,average_record=average_record,FL=FL)