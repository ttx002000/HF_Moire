using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job_v2.jl")


aa= parse.(Int, ARGS)
time_required=aa[2]
cups_required=aa[3]
memory_required=aa[4]

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "MWCv2_$(Int(aa[1]))"

scratch_dir = ENV["SCRATCH"]
misspath=joinpath(scratch_dir, "MWC_v4/data_output$(Int(aa[1]))/missedjobs.jld2")

st=load(misspath)

index=st["index"]
start=1
ee=length(index)
ba_size=300
count=1
while (count-1)*ba_size+1<=ee
  submit_job(filepath, @__DIR__, job_prefix,index[(count-1)*ba_size+start:min(count*ba_size+start,ee)]; time="$(time_required):00:00",cpus_per_task=cups_required,mem=memory_required)
  println((count-1)*ba_size+start,min(count*ba_size+start,ee))
  sleep(5)

  global count+=1
  if mod(count,4)==0
    sleep(0)
  end


end

