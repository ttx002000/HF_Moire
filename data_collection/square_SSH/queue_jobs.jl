using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2
include(joinpath(@__DIR__,"submit_job.jl"))
filepath = joinpath(@__DIR__, "main_test.jl")

job_prefix = "s_SSH"

Nx=20
Ny=20
tpa=1.0
#α=-0.05
#β=-0.05
K=1.0;
KNNN=1.0;
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


for trytime in 1:10, filling in [0.1,0.2], α in collect(-0.0:-0.05:-1.0)
  arguments=Float64.([Nx,Ny,tpa,α,α,K,KNNN,filling,trytime])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="30:00",ntasks=1,mem=8)
end


#=
st=load(joinpath(@__DIR__, "missedjobs.jld2"))
index=st["index"]
for ja in eachindex(index)
  arguments=index[ja]
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="40:00",ntasks=1,mem=8)
end
=#