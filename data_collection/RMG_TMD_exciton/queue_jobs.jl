using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job_v2.jl")


aa= parse.(Int, ARGS)


filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "RMG_TMD_exciton_$(Int(aa[1]))"

scratch_dir = ENV["SCRATCH"]
misspath=joinpath(scratch_dir, "RMG_TMD_exciton/data_output$(Int(aa[1]))/missedjobs.jld2")

st=load(misspath)

index=st["index"]
start=1
ee=length(index)
ba_size=300
count=1
while (count-1)*ba_size+1<=ee
  submit_job(filepath, @__DIR__, job_prefix,index[(count-1)*ba_size+start:min(count*ba_size+start,ee)]; time="2:00:00",ntasks=8,mem=36)
  println((count-1)*ba_size+start,min(count*ba_size+start,ee))
  sleep(5)

  global count+=1
  if mod(count,3)==0
    sleep(0)
  end


end



