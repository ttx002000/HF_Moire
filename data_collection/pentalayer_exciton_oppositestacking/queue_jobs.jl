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

ϵr=5.0

#=
vone=Int(args[1])
vtwo=Int(args[2])
radius=args[3]
num_points=Int(args[4])
uD=args[5]
ildis=args[6]
ϵr=args[7]
trytime=args[8]
=#





st=load("missedjobs.jld2")
index=st["index"]
for ja in eachindex(index)
  arguments=index[ja]
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="50:00",ntasks=4,mem=16)
end
