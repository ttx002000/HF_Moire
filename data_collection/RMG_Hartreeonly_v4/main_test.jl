using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2

# This version allows for partial polarization
include("../../src/operators_RMG_Hartreeonly_v4.jl")

args=parse.(Float64,ARGS)


num_kpoints=Int(args[1])
radius=args[2]
uD=args[3]
temp=args[4]
num_layers=Int(args[5])
ϵr=args[6]
tg1=args[7]
tg2=args[8]
tg3=args[9]
tg4=args[10]
tg_dis=args[11]
bg_dis=args[12]
trytime=Int(args[13])
num_kpoints_OBM=Int(args[14])
Ham_ver=Int(args[15])
file_pos=Int(args[16])

target_density_list=zeros(2,2)
target_density_list=[tg1 tg2;tg3 tg4]


DIIS_density,density_layer_resolved,eout,potential_profile,kinetic_energy,kinetic_energy_resolved,potential_energy,fermi_energy_list=iteration_loop(num_kpoints,radius,uD,temp,
                  num_layers,ϵr,target_density_list,tg_dis,
                  bg_dis,Ham_ver)

OBM,total_M_resolved=get_OBM(num_kpoints_OBM,radius,num_layers,potential_profile,
               target_density_list,temp,Ham_ver)
DOS,DOS_diff=get_DOS(num_kpoints_OBM,radius,num_layers,potential_profile,
               target_density_list,temp,Ham_ver)


scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "RMG_Hartreeonly_v4/data_output$(Int(args[16]))/$(args[1])nk$(args[2])radius$(args[3])uD$(args[4])T$(args[5])nL$(args[6])er$(args[7])tg1$(args[8])tg2$(args[9])tg3$(args[10])tg4$(args[11])tgdis$(args[12])bgdis$(args[13])try$(args[14])nkobm$(args[15])hv.jld2")



jldsave(savepath,density_profile=DIIS_density[1],density_layer_resolved=density_layer_resolved,eout=eout,potential_profile=potential_profile,
        kinetic_energy=kinetic_energy,kinetic_energy_resolved=kinetic_energy_resolved,
         OBM=OBM,total_M_resolved=total_M_resolved,DOS=DOS,DOS_diff=DOS_diff,
        potential_energy=potential_energy,fermi_energy_list=fermi_energy_list)


