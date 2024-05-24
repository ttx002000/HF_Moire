using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "thinfilm"



  submit_job(filepath, @__DIR__, job_prefix; time="1:00:00",ntasks=16,mem=64)
