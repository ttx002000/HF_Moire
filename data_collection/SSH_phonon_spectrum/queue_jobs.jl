using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2
include(joinpath(@__DIR__,"submit_job.jl"))
filepath = joinpath(@__DIR__, "main_test.jl")

job_prefix = "SSH"

#Nx=30
#Ny=30
tper=-0.37
tpa=2
tNNN=-0.08
#α=-0.05
#β=-0.05
#K=1.0;
#KNNN=1.0;
#filling=1.35;
temp=10^(-6)
Nx=80
couplingset=[[-1.0,1.0],[-2.0,2.0],[-1.0,2.0],[-2.0,1.0],[-0.5,2.0],[-2.0,0.5],[-0.5,1.5],[-1.5,0.5]]

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



for filling in [1.25], temp in [0.1], shear in [2.0,1.0,0.5,0.1,0.0], KNNN in [2.0,1.0,0.5,0.1,0.0], K in [2.0,1.0,0.5,0.1]
  α=couplingset[1][1]
  β=couplingset[1][2]
  arguments=Float64.([Nx,Nx,tper,tpa,tNNN,α,β,K,KNNN,filling,temp,shear])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="2:30:00",ntasks=1,mem=32)
end


for filling in [1.25], temp in [0.1], shear in [2.0,1.0,0.5,0.1,0.0], KNNN in [2.0,1.0,0.5,0.1,0.0], K in [2.0,1.0,0.5,0.1]
  α=couplingset[1][1]
  β=couplingset[1][2]
  arguments=Float64.([Nx,Nx,tper,tpa,tNNN,α,β,K,KNNN,filling,temp,shear])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="2:30:00",ntasks=1,mem=32)
end


#=
st=load(joinpath(@__DIR__, "missedjobs.jld2"))
index=st["index"]
for ja in eachindex(index)
  arguments=index[ja]
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="3:00:00",ntasks=1,mem=32)
end


st=load(joinpath(@__DIR__, "missedjobs.jld2"))
index=st["index"]
for ja in eachindex(index)
  arguments=index[ja]
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="3:00:00",ntasks=1,mem=32)
end
=#