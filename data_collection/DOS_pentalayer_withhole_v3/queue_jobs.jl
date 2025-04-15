using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2
include(joinpath(@__DIR__,"submit_job_v2.jl"))
filepath = joinpath(@__DIR__, "main_test.jl")

job_prefix = "DOS_v3"

st=load("missedjobs.jld2")
index=st["index"]
count=0
#=
for ja in 1:length(index)
  global count+=1
  arguments=index[ja]
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="2:30:00",ntasks=16,mem=48)
  if mod(count,300)==0
    sleep(500)
  end
  println(ja)
end
=#

submit_job(filepath, @__DIR__, job_prefix,index; time="4:30:00",ntasks=16,mem=128)



