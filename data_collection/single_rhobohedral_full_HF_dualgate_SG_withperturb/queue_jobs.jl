using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "SG"




st=load("missedjobs.jld2")
index=st["index"]
count=0

for ja in 5100:length(index)
  global count+=1
  arguments=index[ja]
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="1:00:00",ntasks=4,mem=16)
  if mod(count,100)==0
    sleep(300)
  end
  println(ja)
end

