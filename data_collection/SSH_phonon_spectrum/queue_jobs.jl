using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2
include(joinpath(@__DIR__,"submit_job.jl"))
filepath = joinpath(@__DIR__, "main_test.jl")

job_prefix = "SSH"

Nx=30
Ny=30
tper=0.37
tpa=2
tNNN=0.16
#α=-0.05
#β=-0.05
K=1.0;
KNNN=1.0;
#filling=1.35;
temp=10^(-6)


#=
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
=#



for filling in [1.22,1.24,1.25], α in [-0.9,-1.1]
  arguments=Float64.([30,30,tper,tpa,tNNN,α,α,K,KNNN,filling,temp])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="4:00:00",ntasks=1,mem=16)
end
