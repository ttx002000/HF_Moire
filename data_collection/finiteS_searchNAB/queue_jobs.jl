using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "search"


ϕ=0.0 #convert this to radian 
Nq=9.0; 
scale=1.0;
constq=[1.5] #divide this by Nq^2
V0=0.0
vfspace=collect(3.5:0.05:4.95)
for jb in eachindex(vfspace), jc in eachindex(constq)
   arguments=[vf,spin,V0,ϕ]

   submit_job(filepath, @__DIR__, job_prefix,arguments; time="18:00:00",cpus_per_task=36)

end




