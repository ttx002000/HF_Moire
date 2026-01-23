using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames


include("../../src/RMG_stoner_magnetization.jl")

args=parse.(Float64,ARGS)

radius=args[1]
num_kpoints=Int(args[2])
uD=args[3]
num_layers=Int(args[4])
target_density=args[5]
Us=args[6]
ϵr=args[7]
tol=args[8]
trytime=args[9]
file_pos=Int(args[10])




ppp=construct_parameters(radius,num_kpoints,uD,num_layers,target_density,Us,ϵr)


val_record,Gr_record,Δ_current=GR_descent(tol,ppp)











scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "RMG_stoner_magnetization/data_output$(Int(args[10]))/$(args[1])radius$(args[2])nk$(args[3])uD$(args[4])nL$(args[5])tgden$(args[6])Us$(args[7])ϵr$(args[8])tol$(args[9])trytime.jld2")



jldsave(savepath,val_record=val_record,Gr_record=Gr_record,Δ_current=Δ_current)


