using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "finiteS"


ϕ=0.0 #convert this to radian 
Nq=6.0; 
scale=1.0;
constq=[2.0] #divide this by Nq^2
V0=0.0
#vfspace=collect(3.5:0.05:4.95)
spinspace=collect(1.0:0.5:10.0)
for jb in eachindex(vfspace), jc in eachindex(constq)
   arguments=[spinspace[jb],V0,ϕ,Nq,scale,constq[jc]]
for ja in 1:3
     submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="5:00:00",cpus_per_task=36)
end
end




