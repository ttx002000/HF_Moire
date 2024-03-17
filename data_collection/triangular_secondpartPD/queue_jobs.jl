using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "HF"

flux=0.7 #multiply this by pi
ϕ=0.0 #convert this to radian 
Nq=9.0; 
scale=1.0;
constq=[0.1,0.2] #divide this by Nq^2
V0space=[0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5]
for jb in eachindex(V0space), jc in eachindex(constq)
   arguments=[flux,V0space[jb],ϕ,Nq,scale,constq[jc]]
for ja in 1:3
    submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="12:00:00",cpus_per_task=36)
end
end

