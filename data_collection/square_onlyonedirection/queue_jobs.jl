using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra


flux=2 #I multiply it by pi when doing the calculation
Vx=1.0
Vy=0.0
Nq=6.0;
scale=1.0;
constqspace=collect(0.1:0.1:1.5) #I divide it by Nq^2 in the actual calculation


include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "square_HF"


for ja in 1:3
     arguments=Float64.([flux,Vx,Vy,Nq,scale,constq])
     submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="3:00:00",ntasks=6,mem=32)
end