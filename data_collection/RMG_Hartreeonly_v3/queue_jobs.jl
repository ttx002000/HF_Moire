using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job_v2.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "stoner"

aa= parse.(Int, ARGS)

scratch_dir = ENV["SCRATCH"]
misspath=joinpath(scratch_dir, "RMG_Hartreeonly_v3/data_output$(Int(aa[1]))/missedjobs.jld2")

st=load(misspath)

index=st["index"]
ee=length(index)
ba_size=300
count=1
while (count-1)*ba_size+1<=ee
  lo = (count-1)*ba_size + 1
  hi = min(count*ba_size, ee)
 index[lo:hi]
  submit_job(filepath, @__DIR__, job_prefix,index[lo:hi]; time="2:00:00",cpus_per_task=4,mem=16)
  println(lo,hi)
  sleep(5)

  global count+=1
  if mod(count,5)==0
    sleep(3600)
  end


end

