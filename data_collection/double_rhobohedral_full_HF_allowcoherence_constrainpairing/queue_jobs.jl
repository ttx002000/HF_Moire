using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job_v2.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "constrained"



radius=1.6
num_points=41
#uD=20.0
ϵr=5.0
#=
radius=args[1]
num_points=Int(args[2])
uD=args[3]
CNP=args[4]
ϵr=args[5]
ildis=args[6]
trytime=args[7]
file_pos=args[8]
=#

#=
st=load("missedjobs.jld2")
index=st["index"]
count=0
for ja in 1:length(index)
  global count+=1
  arguments=index[ja]
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="3:00:00",ntasks=8,mem=36)
  if mod(count,200)==0
    sleep(500)
  end
  println(count)
end
=#

st=load("missedjobs.jld2")
index=st["index"]
start=1
ee=length(index)
ba_size=300
count=1
while (count-1)*ba_size+1<=ee
  submit_job(filepath, @__DIR__, job_prefix,index[(count-1)*ba_size+start:min(count*ba_size+start,ee)]; time="3:00:00",ntasks=8,mem=36)
  println((count-1)*ba_size+start,min(count*ba_size+start,ee))
  sleep(5)

  global count+=1

  if count==5
    sleep(7200)
  end
end


  #submit_job(filepath, @__DIR__, job_prefix,index[1001:2000]; time="3:00:00",ntasks=8,mem=36)
  #sleep(5)
  #submit_job(filepath, @__DIR__, job_prefix,index[2001:length(index)]; time="3:00:00",ntasks=8,mem=36)