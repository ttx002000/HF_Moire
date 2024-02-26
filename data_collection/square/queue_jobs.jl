using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "flux_threading"


Nx=5;
Ny=6;

for ja in 1:Nx*Ny
     submit_job(filepath, @__DIR__, job_prefix,ja; time="120:00")
end