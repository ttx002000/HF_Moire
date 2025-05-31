using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2
include(joinpath(@__DIR__,"submit_job.jl"))
filepath = joinpath(@__DIR__, "main_test.jl")

job_prefix = "SSH_withshear"




#=
Nx=Int(args[1])
Ny=Int(args[2])
tper=args[3]
tpa=args[4]
tNNN=args[5]
α=args[6]
β=args[7]

K=args[8]
gshear=args[9]


filling=(args[10])
KNNN=args[11]
temp=args[12]
trytime=Int(args[13])
Nelec=Int(round(Nx*Ny*filling))
=#







st=load("missedjobs.jld2")
index=st["index"]
count=0
for ja in 1:length(index)
  global count+=1
  arguments=index[ja]
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="1:00:00",ntasks=1,mem=4)
  if mod(count,200)==0
    sleep(600)
  end
  println(count)
end


#=
st=load("missedjobs.jld2")
index=st["index"]
start=1
ee=length(index)
ba_size=200
count=1
while (count-1)*ba_size+1<=ee
  submit_job(filepath, @__DIR__, job_prefix,index[(count-1)*ba_size+start:min(count*ba_size+start,ee)]; time="20:00:00",ntasks=4,mem=16)
  println((count-1)*ba_size+1,min(count*ba_size+start,ee))
  sleep(5)
  global count+=1
end
=#





