using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "thinfilm"
#=
for Nb in [4,8,12]
  submit_job(filepath, @__DIR__, job_prefix,Nb; time="5:00:00",ntasks=64,mem=128)
end
=#


submit_job(filepath, @__DIR__, job_prefix,[4,0]; time="30:00",ntasks=16,mem=64)
#submit_job(filepath, @__DIR__, job_prefix,8; time="4:00:00",ntasks=32,mem=64)
#submit_job(filepath, @__DIR__, job_prefix,12; time="6:00:00",ntasks=32,mem=128)