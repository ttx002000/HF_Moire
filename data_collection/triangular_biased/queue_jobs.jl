using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "flux_threading"

flux=1.0 #multiply this by pi
#V0=30
ϕ=0.0 #convert this to radian 
Nq=6.0; 
scale=1.0;
constq=0.5 #divide this by Nq^2

V0space=[0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5,6,6.5,7]

for jb in eachindex(V0space)
  arguments=[flux,V0space[jb],ϕ,Nq,scale,constq]
 for ja in 1:3
     submit_job(filepath, @__DIR__, job_prefix,arguments,ja,1; time="200:00",cpus_per_task=8)
 end
end

for jb in eachindex(V0space)
     arguments=[flux,V0space[jb],ϕ,Nq,scale,constq]
    for ja in 1:3
        submit_job(filepath, @__DIR__, job_prefix,arguments,ja,0; time="200:00",cpus_per_task=8)
    end
end


