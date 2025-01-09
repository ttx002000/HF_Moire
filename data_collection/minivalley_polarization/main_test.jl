using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2
using Plots
using CSV,DataFrames

include(joinpath(@__DIR__,"../../src/operators_minivalley.jl"))

args=parse.(Float64,ARGS)


#args=[50.0,5.0,1.0,150,300,0.03,0.3]
uD=args[1]
ϵr=args[2]
pocket_num=Int(args[3])
num_u_grid=Int(args[4])
num_theta_grid=Int(args[5])
density_start=args[6]
density_stop=args[7]



eig_set_final,eig_vec_set_final,k_set_final,theta_set_final,u_set_final,umin,umax,Qlength,density_point,chemical_potential,ugrid,θgrid=sample_states(uD,num_u_grid,num_theta_grid,density_start,density_stop,pocket_num)




formfactors,Coulommatrix,measureone,renormalized_density=get_formfactors(density_point,chemical_potential,
                                                             eig_set_final,eig_vec_set_final,
                                                              k_set_final,theta_set_final,
                                                              u_set_final,ϵr,ugrid,θgrid)



final_energy=calculate_energy(measureone,density_point,renormalized_density,
                                  formfactors,Coulommatrix)


savepath=joinpath(@__DIR__, "data_output/$(args[1])uD$(args[2])er$(args[3])pocket_num$(args[4])num_u_grid$(args[5])num_theta_grid$(args[6])density_start$(args[7])density_stop.jld2")


jldsave(savepath,renormalized_density=renormalized_density,final_energy=final_energy)