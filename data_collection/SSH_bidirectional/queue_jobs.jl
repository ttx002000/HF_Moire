using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2
include(joinpath(@__DIR__,"submit_job.jl"))
filepath = joinpath(@__DIR__, "main_test.jl")

job_prefix = "SSH_bidirection"

Nx=20
Ny=20
tper=0.37
tpa=2
tNNN=0.16
#α=-0.05
#β=-0.05
K=1.0;
KNNN=1.0;
#filling=1.35;
seed=1.0



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
trytime=Int(args[11])
seednum=Int(args[12])
γ=args[13]
Nelec=Int(round(Nx*Ny*filling))
=#




for trytime in 1:20, filling in [1.25], α in [-1.2,-1.1,-1.3], γ in collect(0:0.5:5.0), temp in [0.1]
  arguments=Float64.([Nx,Ny,tper,tpa,tNNN,α,α,K,KNNN,filling,trytime,seed,γ,temp])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="30:00",ntasks=1,mem=8)
end


#=
st=load(joinpath(@__DIR__, "missedjobs.jld2"))
index=st["index"]
for ja in eachindex(index)
  arguments=index[ja]
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="3:00:00",ntasks=1,mem=8)
end
=#


