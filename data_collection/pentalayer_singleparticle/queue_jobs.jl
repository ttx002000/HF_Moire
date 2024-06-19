using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "sing"



θspace=collect(0.95:0.05:1.0)
uDspace=collect(20.0:5.0:60.0)
Nq=15
Vperiodspace=collect(2.5:2.5:20.0)


for θ in θspace, uD in uDspace, jtry in 1:2, Vperiod in Vperiodspace
  arguments=Float64.([θ,Vperiod,uD,jtry])
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