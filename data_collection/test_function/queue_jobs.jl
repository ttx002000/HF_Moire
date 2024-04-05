using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "TDHF"

fluxspace=collect(0.0:0.1:1.0) #multiply this by pi
ϕ=0.0 #convert this to radian 
Nq=3.0; 
scale=1.0;
constq=1.0 #divide this by Nq^2
V0=0.0

for jb in 1:1
  arguments=[fluxspace[jb],V0,ϕ,Nq,scale,constq]
for ja in 1:1, bigQindex in 1:1
    submit_job(filepath, @__DIR__, job_prefix,arguments,ja,bigQindex; time="2:00",cpus_per_task=3)
end
end


