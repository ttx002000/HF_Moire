using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2
include(joinpath(@__DIR__,"submit_job.jl"))
filepath = joinpath(@__DIR__, "main_test.jl")

job_prefix = "DOS"

numsample=2*10^6
θ=1.5
rad=1.0
DOS_n_binnum=50
DOS_E_binnum=5*10^4
Density_start=0.2/4
Density_end=5/4

#=
uD=args[1]
numsample=Int(args[2])
θ=args[3]/180*π
rad=args[4]
DOS_n_binnum=Int(args[5])
DOS_E_binnum=Int(args[6])
Density_start=args[7]
Density_end=args[8]
=#




for trytime in 1:10, filling in [0.1,0.2,0.5], α in collect(-0.0:-0.05:-1.0)
  arguments=Float64.([uD,numsample,θ,rad,DOS_n_binnum,DOS_E_binnum,Density_start,Density_end])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="30:00",ntasks=8,mem=32)
end


#=
st=load(joinpath(@__DIR__, "missedjobs.jld2"))
index=st["index"]
for ja in eachindex(index)
  arguments=index[ja]
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="2:00:00",ntasks=1,mem=8)
end
=#