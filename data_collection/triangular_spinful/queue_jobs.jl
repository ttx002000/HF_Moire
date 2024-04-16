using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "spinful"

fluxspace=[1.0] #multiply this by pi
ϕ=0.0 #convert this to radian 
Nqspace=collect(3.0:1.0:10.0); 
scale=1.0;
constq=1.2 #divide this by Nq^2
V0=0.0
ζspace=[0.5]


for jb in eachindex(fluxspace), jz in eachindex(ζspace), jNq in eachindex(Nqspace)
  arguments=[fluxspace[jb],V0,ϕ,Nqspace[jNq],scale,constq,ζspace[jz]]
 for ja in 1:2
    submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="$(Int(Nqspace[jNq])+4):00:00",cpus_per_task=12,mem=128)
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