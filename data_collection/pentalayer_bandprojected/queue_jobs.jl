using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "pentalayer"


erspace=collect(3.0:0.5:10.0)
θspace=[0.0]
Nband=7
uDspace=[50.0]
Nq=15



for  jtry in 1:15, er in erspace, θ in θspace, uD in uDspace
  arguments=Float64.([er,θ,Nq,uD,Nband,jtry])
    submit_job(filepath, @__DIR__, job_prefix,arguments; time="3:00:00",ntasks=16,mem=32)
end



#=
st=load(joinpath(@__DIR__, "missedjobs.jld2"))
index=st["index"]
for ja in eachindex(index)
  arguments=Float64.(index[ja])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="10:00:00",ntasks=16,mem=64)
end
=#