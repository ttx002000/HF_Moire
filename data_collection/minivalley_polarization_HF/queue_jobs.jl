using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2
include(joinpath(@__DIR__,"submit_job.jl"))
filepath = joinpath(@__DIR__, "main_test.jl")

job_prefix = "minivalley"

#=
uD=args[1]
ϵr=args[2]
pocket_num=Int(args[3])
num_u_grid=Int(args[4])
num_theta_grid=Int(args[5])
density_start=args[6]
density_stop=args[7]
=#

uD=50.0
erspace=[5.0,10.0,15.0]
num_u_grid=250
num_theta_grid=250
density_range=collect(0.05:0.05:1.5)
cutoff=108.0
temp=1.0

#er=5.0








for  er in [1.0], trytime in collect(1:1:10),density_point in eachindex(density_range)
  arguments=Float64.([uD,er,cutoff,num_u_grid,num_theta_grid,density_point,temp,trytime])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="10:00",ntasks=16,mem=64)
end




