using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2
include(joinpath(@__DIR__,"submit_job.jl"))
filepath = joinpath(@__DIR__, "main_test.jl")

job_prefix = "DOS"
uDspace=collect(10.0:0.5:40.0)
numsample=10*10^6
θ=1.5
rad=4.0
DOS_n_binnum=50
DOS_E_binnum=10^5
Density_start=0.2/4
Density_end=4.0/4
Ecutoff=100

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




for uD in uDspace
  arguments=Float64.([uD,numsample,θ,rad,DOS_n_binnum,DOS_E_binnum,Density_start,Density_end,Ecutoff])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="25:00",ntasks=8,mem=8)
end


#=
st=load(joinpath(@__DIR__, "missedjobs.jld2"))
index=st["index"]
for ja in eachindex(index)
  arguments=index[ja]
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="2:00:00",ntasks=1,mem=8)
end
=#