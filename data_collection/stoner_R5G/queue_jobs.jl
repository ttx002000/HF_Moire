using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job_v2.jl")

filepath = joinpath(@__DIR__, "main_test.wl")
job_prefix = "stoner"





st=load("missedjobs.jld2")
index=st["index"]
start=1
ee=length(index)
ba_size=2
count=1
while (count-1)*ba_size+1<=ee
  submit_job(filepath, @__DIR__, job_prefix,index[(count-1)*ba_size+start:min(count*ba_size+start-1,ee)]; time="30:00",cpus_per_task=1,mem=8)
  println((count-1)*ba_size+1,min(count*ba_size+start,ee))
  sleep(5)
  global count+=1
end

