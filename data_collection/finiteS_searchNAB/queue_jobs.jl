using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "search"


ϕspace=collect(0:2.5:60) #convert this to radian 
vfspace=collect(0.5:0.2:5.1)
V0=2.0
spin=2.0

vfspace=collect(3.5:0.05:4.95)
for jb in eachindex(vfspace), jc in eachindex(ϕspace)
   arguments=[vfspace[jb],spin,V0,ϕspace[jc]]

   submit_job(filepath, @__DIR__, job_prefix,arguments; time="20:00",cpus_per_task=1)

end




