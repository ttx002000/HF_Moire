using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "piflux"

fluxspace=[0.7] #multiply this by pi
ϕ=20.0 #convert this to radian 
Nqspace=[9.0]; 
scale=1.0;
constqspace=collect(0.1:0.1:0.8) #divide this by Nq^2
V0space=collect(0.1:0.1:1.5)

for ja in eachindex(fluxspace), jNq in eachindex(Nqspace), constq in constqspace,V0 in V0space
  arguments=[fluxspace[ja],V0,ϕ,Nqspace[jNq],scale,constq]
for jb in 1:3
    submit_job(filepath, @__DIR__, job_prefix,arguments,jb; time="$(Int(Nqspace[jNq])):00:00",ntasks=16,mem=32)
end
end

