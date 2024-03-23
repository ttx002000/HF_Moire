using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "HF"

fluxspace=collect(0.0:0.1:1.0) #multiply this by pi
ϕ=0.0 #convert this to radian 
Nq=9.0; 
scale=1.0;
constq=2.0 #divide this by Nq^2
V0=0.0
for jb in eachindex(fluxspace)
  arguments=[fluxspace[jb],V0,ϕ,Nq,scale,constq]
for ja in 1:4
    submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="20:00:00",cpus_per_task=36)
end
end
