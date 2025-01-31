
using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2
using Plots
using CSV,DataFrames

include(joinpath(@__DIR__,"../../src/operators_SSH_phonon.jl"))

args=parse.(Float64,ARGS)
#args=[20,20,-0.37,2.0,-0.08,-1.0,-1.0,1.0,1.0,1.25,10^(-6),0.1]
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


#=
Λset=get_Lambdaset(Nx,Ny,α,β)
Λmatrix_bandbasis,electron_spectrum,FL=Lambdaset_to_matrix(Nx,Ny,Λset,tper,tpa,tNNN,temp)
Λset=nothing
Keff_momentum=Lambda_to_Keff(Λmatrix_bandbasis,electron_spectrum,FL,temp)
Λmatrix_bandbasis=nothing
Kbare_momentum=barephonon(Nx,Ny,K,KNNN,shearstrength)
spectrum,unperturbed_spectrum,χspectrum=get_spectrum(Kbare_momentum,Keff_momentum,Nx,Ny)
=#


Λset,gmatrix=get_Lambdaset_reduandant(Nx,Ny,α,β)
Λmatrix_bandbasis,electron_spectrum,FL=Lambdaset_to_matrix_redundant(Nx,Ny,Λset,tper,tpa,tNNN,temp)
Λset=nothing
Keff_momentum=Lambda_to_Keff_redundant(Λmatrix_bandbasis,electron_spectrum,FL,temp)
Λmatrix_bandbasis=nothing
Kbare_momentum=barephonon(Nx,Ny,K,KNNN,shearstrength)
spectrum,unperturbed_spectrum,χspectrum=get_spectrum_redundant(Kbare_momentum,Keff_momentum,Nx,Ny,gmatrix)


scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "SSH_phonon/test3/$(Int(args[1]))Nx$(Int(args[2]))Ny$(args[3])tper$(args[4])tpa$(args[5])tNNN$(args[6])alpha$(args[7])beta$(args[8])K$(args[9])KNNN$(args[10])filling$(args[11])temp$(args[12])shear.jld2")
#savepath=joinpath(@__DIR__, "data_output/$(Int(args[1]))Nx$(Int(args[2]))Ny$(args[3])tper$(args[4])tpa$(args[5])tNNN$(args[6])alpha$(args[7])beta$(args[8])K$(args[9])KNNN$(args[10])filling$(args[11])temp$(args[12])shear.jld2")




jldsave(savepath,spectrum=spectrum,unperturbed_spectrum=unperturbed_spectrum,Keff_momentum_redundant=Keff_momentum,Kbare_momentum=Kbare_momentum,χspectrum_redundant=χspectrum,gmatrix=gmatrix)
#jldsave(savepath,spectrum=spectrum,unperturbed_spectrum=unperturbed_spectrum,Keff_momentum=Keff_momentum,Kbare_momentum=Kbare_momentum,χspectrum=χspectrum)