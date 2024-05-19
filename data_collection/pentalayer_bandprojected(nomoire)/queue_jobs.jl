using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "pentalayer"


erspace=collect(3.0:0.5:10.0)
θ=0.77
Nband=7
uD=50.0
Nq=15

#=
ϵr=args[1]
θ=args[2]/180*π
Nq=Int(args[3]);
uD=args[4]
Nband=Int(args[5])
trytime=Int(args[6])
=#


for er in erspace, jtry in 1:10
  arguments=Float64.([er,θ,Nq,uD,Nband,jtry])
    submit_job(filepath, @__DIR__, job_prefix,arguments; time="7:00:00",ntasks=16,mem=64)
end

