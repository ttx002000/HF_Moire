using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "piflux"

fluxspace=collect(0.9:0.02:1.1) #multiply this by pi
ϕ=0.0 #convert this to radian 
Nqspace=collect(3.0:1.0:9.0); 
scale=1.0;
constq=5.0 #divide this by Nq^2
V0=0.0



for ja in eachindex(fluxspace), jNq in 1:4
  arguments=[fluxspace[ja],V0,ϕ,Nqspace[jNq],scale,constq]
for jb in 1:3
    submit_job(filepath, @__DIR__, job_prefix,arguments,jb; time="3:00:00",ntasks=16,mem=32)
end
end

for ja in eachindex(fluxspace), jNq in 5:7
  arguments=[fluxspace[ja],V0,ϕ,Nqspace[jNq],scale,constq]
for jb in 1:3
    submit_job(filepath, @__DIR__, job_prefix,arguments,jb; time="6:00:00",ntasks=16,mem=64)
end
end



#=
st=load(joinpath(@__DIR__, "missedjobs1.jld2"))
index=st["index"]
for ja in eachindex(index)
  arguments=[0.7,index[ja][1],20.0,9.0,1.0,index[ja][3]]
    submit_job(filepath, @__DIR__, job_prefix,arguments,Int(index[ja][4]); time="12:00:00",ntasks=16,mem=32)
end
=#