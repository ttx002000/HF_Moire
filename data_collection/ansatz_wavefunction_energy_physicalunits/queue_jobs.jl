using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2


fluxspace=collect(0.9:0.02:1.1)#I multiply it by pi when doing the calculation
Nqspace=[3.0,4.0,5.0,6.0,7.0,8.0,9.0]
V0=0.0
ϕ=0.0
chistart=0.1
chiend=0.7
Vc=40*π^2/√3


#=
flux=args[1]*π
Nq=Int(args[2])
V0=args[3]
ϕ=args[4]
chistart=args[5]
chiend=args[6]
Vc=args[7]
=#

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "physical_ansatz_energy"
#=
for Nq in 3:6,flux in fluxspace
     arguments=Float64.([flux,Nq,V0,ϕ,chistart,chiend,Vc])
     submit_job(filepath, @__DIR__, job_prefix,arguments; time="1:00:00",ntasks=4,mem=32)
end

for Nq in 7:9,flux in fluxspace
     arguments=Float64.([flux,Nq,V0,ϕ,chistart,chiend,Vc])
     submit_job(filepath, @__DIR__, job_prefix,arguments; time="4:00:00",ntasks=16,mem=256)
end

=#

st=load("missedjobs.jld2")
index=st["index"]
for ja in eachindex(index)
  arguments=index[ja]
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="4:00:00",ntasks=16,mem=128)
end
