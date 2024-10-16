using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2
include(joinpath(@__DIR__,"submit_job.jl"))
filepath = joinpath(@__DIR__, "main_test.jl")

job_prefix = "s_SSH"

Nx=18
Ny=18
tpa=2
#α=-0.05
#β=-0.05
K=0.5;
KNNN=0.2;
#filling=1.35;



#=
Nx=Int(args[1])
Ny=Int(args[2])
tpa=args[3]
α=args[4]
β=args[5]
K=args[6]
KNNN=args[7]
filling=(args[8])
trytime=Int(args[9])
Nelec=Int(round(Nx*Ny*filling))
=#


for trytime in 1:6, filling in [0.2,0.3,0.4], α in [-0.001,-0.01,-0.05,-0.1,-0.2,-0.4,-0.8]
  arguments=Float64.([Nx,Ny,tpa,α,α,K,KNNN,filling,trytime])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="20:00",ntasks=1,mem=8)
end

for trytime in 1:6, filling in [0.2,0.3,0.4], α in [-0.001,-0.01,-0.05,-0.1,-0.2,-0.4,-0.8]
  arguments=Float64.([Nx,Ny,tpa,α,α,K,KNNN,filling,trytime])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="20:00",ntasks=1,mem=8)
end
