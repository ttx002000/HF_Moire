
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2
using Plots
using CSV,DataFrames

include(joinpath(@__DIR__,"../../src/operators_square_SSH.jl"))

args=parse.(Float64,ARGS)
#args=  [15.0, 15.0, 1.0, -0.9, -0.9, 1.0, 1.0, 0.5, 10.0]
Nx=Int(args[1])
Ny=Int(args[2])
tpa=args[3]
α=args[4]
β=args[5]
K=args[6]
KNNN=args[7]
filling=(args[8])
trytime=Int(args[9])
Nelec=Int(round(Nx*Ny*filling))

H0, orbital_id, phonon_id, px_xbond, px_ybond,  NNN_sp_d1, NNN_sp_d2=initialize(Nx,Ny,tpa)


phonon_coor, Hph, grad_old, E_old, Eelec_new, Egap,ave_npa=iteration(Nx,Ny,Nelec,px_xbond,px_ybond,NNN_sp_d1,NNN_sp_d2,orbital_id,phonon_id,α,β,K,KNNN,H0)
dis_x,dis_y,max_record,average_record=resh_phonon(phonon_coor,phonon_id,Nx,Ny)

FFF=eigen(H0+Hph)
spectrum=FFF.values
savepath=joinpath(@__DIR__, "data_output/$(Int(args[1]))Nx$(Int(args[2]))Ny$(args[3])tpa$(args[4])alpha$(args[5])beta$(args[6])K$(args[7])KNNN$(args[8])filling$(Int(args[9]))try.jld2")
jldsave(savepath,phonon_id=phonon_id,orbital_id=orbital_id,phonon_coor=phonon_coor,E_total=E_old,Eelec=Eelec_new,spectrum=spectrum,Egap=Egap,max_record=max_record,average_record=average_record,ave_npa=ave_npa)