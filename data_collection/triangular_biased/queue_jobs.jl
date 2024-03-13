using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "HF"

flux=0.7 #multiply this by pi
V0=0.6
ϕ=0.0 #convert this to radian 
Nq=9.0; 
scale=1.0;
#constq=1.3 #divide this by Nq^2

constqspace=[0.6,0.7,0.8,0.9,1.0]
#=
for jb in eachindex(V0space)
  arguments=[flux,V0space[jb],ϕ,Nq,scale,constq]
 for ja in 1:3
     submit_job(filepath, @__DIR__, job_prefix,arguments,ja,1; time="400:00",cpus_per_task=36)
 end
end
=#

for jb in eachindex(constqspace)
     arguments=[flux,V0,ϕ,Nq,scale,constqspace[jb]]
    for ja in 1:3
        submit_job(filepath, @__DIR__, job_prefix,arguments,ja,1; time="30:00:00",cpus_per_task=36)
    end
end


