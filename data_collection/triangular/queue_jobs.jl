using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "piflux"

fluxspace=[2.0] #multiply this by pi
ϕ=0.0 #convert this to radian 
Nqspace=collect(3.0:1.0:9.0); 
scale=1.0*√2;
constqspace=collect(2.0:0.5:6.0) #divide this by Nq^2
V0=0.0

for ja in eachindex(fluxspace), jNq in eachindex(Nqspace), constq in constqspace
  arguments=[fluxspace[ja],V0,ϕ,Nqspace[jNq],scale,constq]
for jb in 1:2
    submit_job(filepath, @__DIR__, job_prefix,arguments,jb; time="$(Int(Nqspace[jNq])):00:00",ntasks=Int(Nqspace[jNq])^2,mem=128)
end
end

