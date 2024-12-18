using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2
include(joinpath(@__DIR__,"submit_job.jl"))
filepath = joinpath(@__DIR__, "main_test.jl")

job_prefix = "DOS"
uDspace=collect(-50.0:1.0:50.0)
numsample=8*10^7
θ=1.5
rad=1.5
DOS_n_binnum=120
DOS_E_binnum=6*10^4
Density_start=-1.0
Density_end=1.0
Ecutoff=200

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




for uD in collect(-60.0:1.0:-50.0)
  arguments=Float64.([uD,numsample,θ,rad,DOS_n_binnum,DOS_E_binnum,Density_start,Density_end,Ecutoff])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="1:30:00",ntasks=16,mem=32)
end

for uD in collect(50.0:1.0:60.0)
  arguments=Float64.([uD,numsample,θ,rad,DOS_n_binnum,DOS_E_binnum,Density_start,Density_end,Ecutoff])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="1:30:00",ntasks=16,mem=32)
end


#=
st=load(joinpath(@__DIR__, "missedjobs.jld2"))
index=st["index"]
for ja in eachindex(index)
  arguments=index[ja]
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="30:00",ntasks=8,mem=32)
end

=#

