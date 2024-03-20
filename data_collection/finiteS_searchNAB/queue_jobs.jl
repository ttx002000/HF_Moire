using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "search"


ϕspace=collect(0.0:5.0:60.0) #convert this to radian 
vfspace=collect(0.5:0.2:4.5)
V0=[1.0,2.0,3.0]
spin=3.0


for jb in eachindex(vfspace), jc in eachindex(ϕspace), jd in eachindex(V0)
   arguments=[vfspace[jb],spin,V0[jd],ϕspace[jc]]

   submit_job(filepath, @__DIR__, job_prefix,arguments; time="10:00",cpus_per_task=1)

end




