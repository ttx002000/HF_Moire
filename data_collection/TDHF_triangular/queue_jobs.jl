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

for jb in [0.0]
  arguments=[jb,V0,ϕ,Nq,scale,constq]
for ja in 1:1, bigQindex in 1:Nq^2
    submit_job(filepath, @__DIR__, job_prefix,arguments,ja,bigQindex; time="15:00",cpus_per_task=4)
end
end


#=
miss=load(joinpath(@__DIR__, "missedindex.jld2"))
index=miss["index"]
for ja in eachindex(index)
  arguments=[index[ja][1],V0,ϕ,Nq,scale,constq]
  submit_job(filepath, @__DIR__, job_prefix,arguments,Int(index[ja][3]),Int(index[ja][2]); time="25:00",cpus_per_task=1)
end
=#