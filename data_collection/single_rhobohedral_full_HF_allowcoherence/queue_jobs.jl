using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "single_RMG"




st=load("missedjobs.jld2")
index=st["index"]
count=0
for ja in eachindex(index)
  global count+=1
  arguments=index[ja]
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="1:30:00",ntasks=4,mem=24)
  if mod(count,200)==0
    sleep(1800)
  end
  println(count)
end

