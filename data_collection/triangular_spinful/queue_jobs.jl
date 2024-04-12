using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "HF"

fluxspace=[1.0] #multiply this by pi
ϕ=0.0 #convert this to radian 
Nq=6.0; 
scale=1.0;
constq=1.0 #divide this by Nq^2
V0=0.0
ζspace=collect(0.0:0.2:1.0)


for jb in eachindex(fluxspace), jz in eachindex(ζspace)
  arguments=[fluxspace[jb],V0,ϕ,Nq,scale,constq,ζspace[jz]]
 for ja in 1:8
    submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="10:00:00",cpus_per_task=5)
 end
end

#=
st=load(joinpath(@__DIR__, "missedjobs.jld2"))
index=st["index"]
for ja in eachindex(index)
    arguments=[index[ja][1],V0,ϕ,Nq,scale,constq,ζ]
    trytime=Int(index[ja][2])
    submit_job(filepath, @__DIR__, job_prefix,arguments,trytime; time="38:00:00",cpus_per_task=36)
end
=#