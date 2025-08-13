using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames


include("../../src/operators_graphenesusceptability.jl")

args=parse.(Float64,ARGS)

uD=args[1]
whichstack=Int(args[2])
radius=args[3]
Γ=args[4]
ϵrange=Int(args[5])
num_points=Int(args[6])
file_pos=args[7]

loadpath=joinpath(scratch_dir, "suseptability_graphene/data_output$(Int(args[7]))/epsilon_range/epsilon_range$(ϵrange).jld2")
st1=load(loadpath)
ϵspace=st1["ϵrange"]

Fz=main(radius,Γ,uD,whichstack,ϵspace,num_points)

scratch_dir = ENV["SCRATCH"]

savepath=joinpath(scratch_dir, "suseptability_graphene/data_output$(Int(args[7]))/$(args[1])uD$(args[2])whichstack$(args[3])radius$(args[4])Gamma$(args[5])ϵrange$(args[6])num_points.jld2")




jldsave(savepath,Fz=Fz)


