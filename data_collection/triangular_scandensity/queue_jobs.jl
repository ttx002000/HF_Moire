using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "HF"

#fluxspace=[2.0] #multiply this by pi
ϕ=0.0 #convert this to radian 
Nq=9.0; 
scalespace=collect(0.3:0.1:1.8);
#constq=2.4 #divide this by Nq^2
V0=0.0

for jb in eachindex(scalespace)
  flux=1.0*scale^2
  constq=1.0*scale^2
  arguments=[flux,V0,ϕ,Nq,scalespace[jb],constq]
 for ja in 1:4
    submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="10:00:00",cpus_per_task=3)
 end
end

for jb in eachindex(scalespace)
  flux=1.0*scale^2
  constq=1.2*scale^2
  arguments=[flux,V0,ϕ,Nq,scalespace[jb],constq]
 for ja in 1:4
    submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="10:00:00",cpus_per_task=3)
 end
end

