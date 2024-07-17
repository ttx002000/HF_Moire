using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

scale=1.0
fluxspace=collect(0.9:0.02:1.1)#I multiply it by pi when doing the calculation
Nqspace=[3.0,4.0,5.0,6.0,7.0,8.0,9.0]
V0=0.0
ϕ=0.0
chistart=0.1
chiend=1.0
constq=5.0 #I divide it by Nq^2 in the actual calculation


#=
scale=args[1]
flux=args[2]*π
Nq=Int(args[3])
V0=args[4]
ϕ=args[5]
chistart=args[6]
chiend=args[7]
constq=args[8]/Nq^2
=#

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "ansatz_energy"

#=
for Nq in 3:6,flux in fluxspace
     arguments=Float64.([scale,flux,Nq,V0,ϕ,chistart,chiend,constq])
     submit_job(filepath, @__DIR__, job_prefix,arguments; time="1:00:00",ntasks=4,mem=32)
end

for Nq in 7:9,flux in fluxspace
     arguments=Float64.([scale,flux,Nq,V0,ϕ,chistart,chiend,constq])
     submit_job(filepath, @__DIR__, job_prefix,arguments; time="4:00:00",ntasks=16,mem=256)
end
=#
st=load("missedjobs.jld2")
index=st["index"]
for ja in eachindex(index)
  arguments=index[ja]
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="4:00:00",ntasks=32,mem=256)
end