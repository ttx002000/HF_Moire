using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "search"


ϕspace=collect(0.0:5.0:60.0) #convert this to radian 
vfspace=collect(0.5:0.2:4.5)
V0=[1.0,1.5,2.0,2.5,3.0,3.5]
spinspace=[1.5]
bandspace=[2.0]

for jb in eachindex(vfspace), jc in eachindex(ϕspace), jd in eachindex(V0), spin in spinspace, bandindex in bandspace
   arguments=[vfspace[jb],spin,V0[jd],ϕspace[jc],bandindex]

   submit_job(filepath, @__DIR__, job_prefix,arguments; time="25:00",cpus_per_task=1)

end



#st=load(joinpath(@__DIR__, "missedindex.jld2"))
#missedindex=st["missedindex"]
#for tt=eachindex(missedindex)
 #  arguments=[missedindex[tt][3],missedindex[tt][1],missedindex[tt][2],missedindex[tt][4],missedindex[tt][5]]
  # submit_job(filepath, @__DIR__, job_prefix,arguments; time="25:00",cpus_per_task=1)
#end




