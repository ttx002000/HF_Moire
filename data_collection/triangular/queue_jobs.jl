using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "HF"

fluxspace=[2.0] #multiply this by pi
ϕ=0.0 #convert this to radian 
Nqspace=[6.0]; 
scale=1.0*√2;
constq=2.0 #divide this by Nq^2
V0=0.0

for ja in eachindex(fluxspace), jNq in eachindex(Nqspace)
  arguments=[fluxspace[ja],V0,ϕ,Nqspace[jNq],scale,constq]
for ja in 1:5
    submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="$(Int(Nqspace[jNq])+2):00:00",cpus_per_task=5)
end
end

constq=2.4
for ja in eachindex(fluxspace), jNq in eachindex(Nqspace)
  arguments=[fluxspace[ja],V0,ϕ,Nqspace[jNq],scale,constq]
for ja in 1:5
    submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="$(Int(Nqspace[jNq])+2):00:00",cpus_per_task=5)
end
end