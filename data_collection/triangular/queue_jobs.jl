using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "flux_threading"

flux=3*π
V0=0.5*exp(π/2)
ϕ=0.0
Nq=3;
scale=1.0;
constq=0.1/Nq^2




arguments=[flux,V0,ϕ,Nq,scale,constq]
for ja in 1:1
     submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="120:00",cpus_per_task=5)
end