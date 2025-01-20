using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "exciton"



radius=1.6
num_points=41
#uD=20.0
ϵr=5.0
temp=0.1
target_density=-0.02


#=
radius=args[1]
num_points=Int(args[2])
uD=args[3]
target_density=args[4]
ϵr=args[5]
temp=args[6]
trytime=args[7]

=#

for uD in [10.0],ϵr in [5.0,10.0], trytime in 21:40
  arguments=Float64.([radius,num_points,uD,target_density,ϵr,temp,trytime])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="120:00",ntasks=16,mem=24)
end

for uD in [30.0],ϵr in [5.0,10.0], trytime in 1:40
  arguments=Float64.([radius,num_points,uD,target_density,ϵr,temp,trytime])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="120:00",ntasks=16,mem=24)
end
#=
st=load("missedjobs.jld2")
index=st["index"]
for ja in eachindex(index)
  arguments=index[ja]
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="50:00",ntasks=4,mem=24)
end
=#
