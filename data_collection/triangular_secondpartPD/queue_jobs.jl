using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2
include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "HF"

flux=0.7 #multiply this by pi
ϕ=0.0 #convert this to radian 
Nq=9.0; 
scale=1.0;
#constq=collect(0.3:0.1:0.8) #divide this by Nq^2
#V0space=collect(0.1:0.1:1.5)
#for jb in eachindex(V0space), jc in eachindex(constq)
 #  arguments=[flux,V0space[jb],ϕ,Nq,scale,constq[jc]]
#for ja in 1:3
 #   submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="18:00:00",cpus_per_task=36)
#end
#end
st=load(joinpath(@__DIR__, "missedjobs.jld2"))
index=st["index"]
for ja in eachindex(index)
    arguments=[flux,index[ja][1],ϕ,Nq,scale,index[ja][2]]
    trytime=Int(index[ja][3])
    submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="28:00:00",cpus_per_task=36)
end

constq=collect(0.9:0.1:1.0) #divide this by Nq^2
V0space=collect(0.1:0.1:1.5)
for jb in eachindex(V0space), jc in eachindex(constq)
   arguments=[flux,V0space[jb],ϕ,Nq,scale,constq[jc]]
for ja in 1:3
   submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="28:00:00",cpus_per_task=36)
end
end