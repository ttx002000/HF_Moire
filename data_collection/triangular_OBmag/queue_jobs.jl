using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "piflux"

flux=1.0 #multiply this by pi
ϕ=0.0 #convert this to radian 
Nqspace=[6.0,9.0]; 
scale=1.0;
constqspace=[1.0,1.5,2.0] #divide this by Nq^2
V0=0.0

for  jNq in eachindex(Nqspace), constq in constqspace
  arguments=[flux,V0,ϕ,Nqspace[jNq],scale,constq]
for ja in 1:2
    submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="$(Int(Nqspace[jNq])):00:00",ntasks=Int(Nqspace[jNq])^2,mem=256)
end
end


flux=1.0 #multiply this by pi
ϕ=0.0 #convert this to radian 
Nqspace=[6.0,9.0]; 
scale=1.0;
constqspace=[0.0,0.2] #divide this by Nq^2
V0space=[0.5,1.0,1.5]

for  jNq in eachindex(Nqspace), constq in constqspace, jV in eachindex(V0space)
  arguments=[flux,V0space[jV],ϕ,Nqspace[jNq],scale,constq]
for ja in 1:2
    submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="$(Int(Nqspace[jNq])):00:00",ntasks=Int(Nqspace[jNq])^2,mem=256)
end
end
