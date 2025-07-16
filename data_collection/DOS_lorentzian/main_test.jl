using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames


include("../../src/operators_DOS_Lorentzian.jl")

#args=parse.(Float64,ARGS)
args=[50.0,50.0,1.0,1]

uD=args[1]
Nq=Int(args[2])
Γ=args[3]
type=Int(args[4])
file_pos=Int(args[5])

ac=0.246
R1=ac*[1,0]
R2=ac*[1/2,√3/2]
G1=2π/ac*[1,-1/√3]
G2=2π/ac*[0,2/√3]
KGr=4π/(3*ac)*[1.0,0.0]
Area=Nq^2*ac^2*√3/2


Elist,DOS=record_values(Nq,type,Γ,uD,Area)



scratch_dir = ENV["SCRATCH"]

savepath=joinpath(scratch_dir, "DOS_Lorenztian/data_output$(Int(args[5]))/$(args[1])uD$(args[2])Nq$(args[3])Gamma$(args[4])type.jld2")




jldsave(savepath,Elist,DOS)


