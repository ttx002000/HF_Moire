using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra





include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "HF_TB"


for trytime in collect(1.0:1.0:30.0)
     arguments=[Int(trytime)]
     submit_job(filepath, @__DIR__, job_prefix,arguments; time="3:00:00",ntasks=16,mem=16)
end