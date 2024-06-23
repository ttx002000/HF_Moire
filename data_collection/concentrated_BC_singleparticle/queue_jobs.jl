using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "sing"


#ϕspace=collect(0.0:2.5:60.0)
ϕ=0.0
spinspace=collect(0.0:0.5:10.0)
Nq=180
V0=2.0
mass=1.0
M=0.5
scale=9.0


#=
Nq=Int(args[1])
spin=args[2]
scale=args[3]
ϕ=args[4]/180*π
V0=args[5]
mass=args[6]
M=args[7]

=#

for  jtry in 1:2, spin in spinspace
  arguments=Float64.([Nq,spin,scale,ϕ,V0,mass,M,jtry])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="20:00",ntasks=1,mem=64)
end
for  jtry in 1:2, spin in spinspace
  arguments=Float64.([Nq,spin,scale,ϕ,3.0,mass,M,jtry])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="20:00",ntasks=1,mem=64)
end
