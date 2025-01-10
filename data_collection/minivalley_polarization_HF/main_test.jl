using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2
using Plots
using CSV,DataFrames

include(joinpath(@__DIR__,"../../src/operators_minivalley_HF.jl"))

args=parse.(Float64,ARGS)


#args=[50.0,5.0,50.0,100,100,0.5,0.5,1.0]
uD=args[1]
ϵr=args[2]
cutoff=args[3]
num_u_grid=Int(args[4])
num_theta_grid=Int(args[5])
density_point=args[6]*0.01
temp=args[7]
trytime=args[8]



eig_set_final,eig_vec_set_final,k_set_final,theta_set_final,u_set_final,umin,umax,Qlength,chemical_potential,ugrid,θgrid=sample_states(uD,num_u_grid,num_theta_grid,density_point,cutoff)





ϵr=1.0
formfactors,measureone=get_formfactors( eig_set_final,eig_vec_set_final,
                                        k_set_final,theta_set_final,
                                        u_set_final,ϵr,ugrid,θgrid)
scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "minivalley/FF_$(args[1])uD1.0er$(args[3])cutoff$(args[4])num_u_grid$(args[5])num_theta_grid.jld2")
jldsave(savepath,formfactors=formfactors,measureone=measureone)                                        


#=
scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "minivalley/FF_$(args[1])uD1.0er$(args[3])cutoff$(args[4])num_u_grid$(args[5])num_theta_grid.jld2")

s1=load(savepath)
formfactors=s1["formfactors"]/ϵr
measureone=s1["measureone"]


fermifactor_final,fermilevel_final,renormalized_density_final,energy_final,quasi_particle_energy_final=do_iterations(measureone,density_point,formfactors,temp,eig_set_final)

scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "minivalley/data_output2/$(args[1])uD$(args[2])er$(args[3])cutoff$(args[4])num_u_grid$(args[5])num_theta_grid$(args[6])density$(args[7])temp$(args[8])try.jld2")


jldsave(savepath,fermifactor=fermifactor_final,energy=energy_final,fermilevel=fermilevel_final,
                renormalized_density=renormalized_density_final,quasi_particle_energy=quasi_particle_energy_final,
                eig_set_final= eig_set_final,k_set_final=k_set_final,theta_set_final=theta_set_final,u_set_final=u_set_final)
=#