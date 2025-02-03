using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "double_pentalayer"



radius=1.6
num_points=41
#uD=20.0
ϵr=5.0
#=
radius=args[1]
num_points=Int(args[2])
uD=args[3]
CNP=args[4]
ϵr=args[5]
ildis=args[6]
trytime=args[7]
file_pos=args[8]
=#



st=load("missedjobs.jld2")
index=st["index"]
for ja in 1:500
  arguments=index[ja]
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="2:30:00",ntasks=8,mem=64)
end

