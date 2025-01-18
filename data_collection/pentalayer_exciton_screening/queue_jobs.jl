using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "exciton"


vone=1
vtwo=1
radius=1.8
num_points=100
#uD=30.0
ildis=1.0
ϵr=5.0
stackingtwo=1
CNP=0.0


#=
vone=Int(args[1])
vtwo=Int(args[2])
radius=args[3]
num_points=Int(args[4])
uD=args[5]
ildis=args[6]
ϵr=args[7]
stackingtwo=Int(args[8])
CNP=args[9]-4*uD
trytime=args[10]
=#

for uD=[10.0,20.0,30.0,40.0,50.0], trytime in [1.0]
  arguments=[vone,-1,radius,num_points,uD,ildis,ϵr,stackingtwo,CNP,trytime]
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="1:00:00",ntasks=16,mem=32)
end


#=
st=load("missedjobs.jld2")
index=st["index"]
for ja in eachindex(index)
  arguments=index[ja]
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="50:00",ntasks=4,mem=16)
end
=#