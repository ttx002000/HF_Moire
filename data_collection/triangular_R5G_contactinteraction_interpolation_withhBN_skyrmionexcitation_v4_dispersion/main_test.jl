using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames

# Basically the same as v1, but this one specifically deals with Nq=1 case, and also spinless
include("../../src/operators_R5G_contactinteraction_interpolation_withhBN_skyrmionexcitation_v4.jl")

args=parse.(Float64,ARGS)



NL=args[1]
θ=args[2]/180*pi;
constq=args[3]
ϵr=args[4]
uD=args[5]
filling=Int(args[6])
gcutoff=args[7]
λ=args[8]
trytimes=Int(args[9])
enlarge_factor=Int(args[10])
V0_hBN=args[11]
V1_hBN=args[12]
ψ_hBN=args[13]
V2_scalar=args[14]
ϕ=args[15]/180*π
filepos=Int(args[16])





manybody_overlap,H_matrixelement,shift_set=get_manybodyoverlap(args)




                                                      


scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "triangle_R5G_contact_interpolation_withhBN_skyrmionexcitation_v4/data_output$(Int(args[16]))/dispersion/result/$(args[1])NL$(args[2])theta$(args[3])constq$(args[4])ϵr$(args[5])uD$(args[6])filling$(args[7])cutoff$(args[8])lambda$(args[9])trytime$(args[10])enlarge$(args[11])V0_hBN$(args[12])V1_hBN$(args[13])ψ_hBN$(args[14])V2_scalar$(args[15])ϕ.jld2")



jldsave(savepath,manybody_overlap=manybody_overlap,H_matrixelement=H_matrixelement,shift_set=shift_set)


