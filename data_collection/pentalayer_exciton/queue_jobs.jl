using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "sing"


vone=1
vtwo=1
radius=1.8
num_points=100
uD=30.0
#CNP=-100.0
ϵr=5.0

#=
vone=Int(args[1])
vtwo=Int(args[2])
radius=args[3]
num_points=Int(args[4])
uD=args[5]
CNP=args[6]
ϵr=args[7]
trytime=args[8]
=#

#=
for jtry in 1:20, CNP in collect(-60.0:-10.0:-200.0)
  arguments=Float64.([vone,-1,radius,num_points,uD,CNP,ϵr,jtry])
    submit_job(filepath, @__DIR__, job_prefix,arguments; time="1:00:00",ntasks=8,mem=32)
end
=#
#=
for jtry in 1:20, CNP in collect(-60.0:-10.0:-200.0)
  arguments=Float64.([vone,vtwo,radius,num_points,uD,CNP,ϵr,jtry])
    submit_job(filepath, @__DIR__, job_prefix,arguments; time="30:00",ntasks=8,mem=32)
end
=#
#=
for jtry in 1:20, CNP in collect(-100.0:-10.0:-250.0)
  arguments=Float64.([vone,-1,radius,num_points,50.0,CNP,ϵr,jtry])
    submit_job(filepath, @__DIR__, job_prefix,arguments; time="30:00",ntasks=8,mem=32)
end
=#



st=load("missedjobs.jld2")
index=st["index"]
for ja in eachindex(index)
  arguments=index[ja]
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="30:00",ntasks=8,mem=32)
end
