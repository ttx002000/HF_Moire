using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra


flux=2 #I multiply it by pi when doing the calculation
V0=0.5
ϕ=0.0 #I convert this degree to randian
Nq=3.0;
scale=1.0;
constq=0.1 #I divide it by Nq^2 in the actual calculation
arguments=[flux,V0,ϕ,Nq,scale,constq]

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "HF"


for ja in 1:5
     submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="120:00",cpus_per_task=5)
end