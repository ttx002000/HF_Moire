using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2
include(joinpath(@__DIR__,"submit_job.jl"))
filepath = joinpath(@__DIR__, "main_test.jl")

job_prefix = "DOS"
uDspace=collect(-60.0:1.0:60.0)
numsample=2*10^6
θ=1.5
DOS_n_binnum=120
DOS_E_binnum=4*10^4
Density_start=-1.0
Density_end=1.0
Ecutoff=200

#=
uD=args[1]
numsample=Int(args[2])
θ=args[3]/180*π
DOS_n_binnum=Int(args[4])
DOS_E_binnum=Int(args[5])
Density_start=args[6]
Density_end=args[7]
Ecutoff=args[8]
=#




for uD in uDspace
  arguments=Float64.([uD,numsample,θ,DOS_n_binnum,DOS_E_binnum,Density_start,Density_end,Ecutoff])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="5:00:00",ntasks=32,mem=64)
end


#=
st=load(joinpath(@__DIR__, "missedjobs.jld2"))
index=st["index"]
for ja in eachindex(index)
  arguments=index[ja]
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="30:00",ntasks=8,mem=32)
end

=#

