using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "HF"

flux=1.0 #multiply this by pi
ϕ=0.0 #convert this to radian 
Nq=15.0; 
scale=1.0;
constq=[0.0,0.5,1.0,0.7,0.2] #divide this by Nq^2
V0space=[1.5]
for jb in eachindex(V0space)
   arguments=[flux,V0space[jb],ϕ,Nq,scale,constq]
for ja in 1:3
     submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="1000:00",cpus_per_task=48)
end
end




