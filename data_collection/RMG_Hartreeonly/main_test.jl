using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames


include("../../src/operators_RMG_Hartreeonly.jl")

args=parse.(Float64,ARGS)


num_kpoints=Int(args[1])
radius=args[2]
uD=args[3]
temp=args[4]
num_layers=Int(args[5])
ϵr=args[6]
target_density=args[7]
tg_dis=args[8]
bg_dis=args[9]
active_flavor=Int(args[10])
trytime=Int(args[11])
num_kpoints_OBM=Int(args[12])
Ham_ver=Int(args[13])
file_pos=Int(args[14])



DIIS_density,eout,potential_profile,kinetic_energy=iteration_loop(num_kpoints,radius,uD,temp,
                  num_layers,ϵr,target_density,tg_dis,
                  bg_dis,active_flavor,Ham_ver)

 OBM=get_OBM(num_kpoints_OBM,radius,num_layers,potential_profile,
               target_density,temp,active_flavor,Ham_ver)



scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "RMG_Hartreeonly/data_output$(Int(args[14]))/$(args[1])nk$(args[2])radius$(args[3])uD$(args[4])T$(args[5])nL$(args[6])er$(args[7])tgden$(args[8])tgdis$(args[9])bgdis$(args[10])active$(args[11])try$(args[12])nkobm$(args[13])hv.jld2")



jldsave(savepath,density_profile=DIIS_density[1],eout=eout,potential_profile=potential_profile,
        kinetic_energy=kinetic_energy, OBM=OBM)


