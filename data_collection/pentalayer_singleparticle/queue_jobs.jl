using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "sing"



#θspace=collect(0.95:0.05:1.0)
uDspace=collect(20.0:2.5:60.0)
Nq=15
#Vperiodspace=collect(2.5:2.5:20.0)
#periodspace=collect(8.0:0.35:15.0)
phaseanglespace=collect(0.0:2.5:60.0)
Vperiod=10.0
period=10.0

for  uD in uDspace, jtry in 1:2, phaseangle in phaseanglespace
  arguments=Float64.([0.0,Vperiod,uD,period,phaseangle,jtry])
    submit_job(filepath, @__DIR__, job_prefix,arguments; time="15:00",ntasks=1,mem=16)
end
#=
st=load("missedjobs.jld2")
index=st["index"]
for ja in eachindex(index)
  arguments=index[ja]
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="10:00",ntasks=1,mem=8)
end
=#