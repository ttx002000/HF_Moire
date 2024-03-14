using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test2.jl")
job_prefix = "HF"

flux=0.7 #multiply this by pi
V0=0.6
ϕ=0.0 #convert this to radian 
Nq=9.0; 
scale=1.0;
#constq=1.3 #divide this by Nq^2

constqspace=[1.32,1.34,1.36,1.38,1.4,1.42,1.44,1.46,1.48,1.5,1.52,1.54]
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


