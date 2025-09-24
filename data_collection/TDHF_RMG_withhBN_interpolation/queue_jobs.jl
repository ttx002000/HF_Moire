using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2


include("submit_job_v2.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "bandproject"



st=load("missedjobs.jld2")
index=st["index"]
start=1
ee=length(index)
ba_size=300
count=1
while (count-1)*ba_size+1<=ee
  submit_job(filepath, @__DIR__, job_prefix,index[(count-1)*ba_size+start:min(count*ba_size+start,ee)]; time="6:00:00",ntasks=16,mem=128)
  println((count-1)*ba_size+start,min(count*ba_size+start,ee))
  sleep(5)

  global count+=1
  if mod(count,5)==0
    sleep(0)
  end


end


