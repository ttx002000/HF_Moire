
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2



include(joinpath(@__DIR__,"../../src/operators_SSH_withshear_v3.jl"))

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
trytime=Int(args[14])
filepos=Int(args[15])

Nelec=Nx*Ny*filling

H0, orbital_id, phonon_id, px_xbond, px_ybond, py_xbond, py_ybond, NNN_sp_d1, NNN_sp_d2=initialize(Nx,Ny,tper,tpa,tNNN)


phonon_coor, Hph, grad_old, E_old, Eelec_new, Egap,ave_npa,FL,free_energy,record_cal=iteration(Nx,Ny,Nelec,px_xbond,px_ybond,py_xbond,py_ybond,NNN_sp_d1,NNN_sp_d2,orbital_id,phonon_id,α,β,K,KNNN,H0,gshear,temp,stop_standard)
dis_x,dis_y,max_record,average_record=resh_phonon(phonon_coor,phonon_id,Nx,Ny)


FFF=eigen(H0+Hph)
spectrum=FFF.values

scratch_dir = ENV["SCRATCH"]

savepath=joinpath(scratch_dir, "SSH_v3/test$(Int(args[15]))/$(Int(args[1]))Nx$(Int(args[2]))Ny$(args[3])tper$(args[4])tpa$(args[5])tNNN$(args[6])alpha$(args[7])beta$(args[8])K$(args[9])gshear$(args[10])filling$(args[11])KNNN$(args[12])temp$(args[13])stop$(Int(args[14]))try.jld2")
#savepath=joinpath(@__DIR__,"test_data_output/test.jld2")
jldsave(savepath,phonon_id=phonon_id,orbital_id=orbital_id,
         phonon_coor=phonon_coor,E_total=E_old,Eelec=Eelec_new,
         free_energy=free_energy,spectrum=spectrum,
         Egap=Egap,ave_npa=ave_npa,max_record=max_record,average_record=average_record,FL=FL,grad=grad_old,record_cal=record_cal)
         
