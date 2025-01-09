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
num_theta_grid=300
density_start=0.05
density_stop=0.5

#pocket_num=1
#er=5.0








for  er in [5.0,10.0,15.0], pocket_num in [1,2,3]
  arguments=Float64.([uD,er,pocket_num,num_u_grid,num_theta_grid,density_start,density_stop])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="2:00:00",ntasks=16,mem=64)
end




