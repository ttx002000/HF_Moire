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
density_range=collect(0.05:0.025:1.5)
cutoff=108.0 #use 75.0 if uD=20.0, use 120.0 if uD=20.0
temp=1.0

#er=5.0



#=
for  er in [5.0,10.0,15.0], trytime in collect(1:1:10), jb in eachindex(density_range)
  arguments=Float64.([uD,er,cutoff,num_u_grid,num_theta_grid,density_range[jb],temp,trytime])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="90:00",ntasks=4,mem=64)
end
=#

#=
for  er in [5.0,10.0,15.0], trytime in collect(1:1:1), jb in [1]
  arguments=Float64.([20.0,er,70.0,200,num_theta_grid,density_range[jb],temp,trytime])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="20:00",ntasks=16,mem=64)
end
=#
#=
for  er in [5.0,10.0,15.0], trytime in collect(1:1:1), jb in [1]
  arguments=Float64.([80.0,er,120.0,num_u_grid,num_theta_grid,density_range[jb],temp,trytime])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="90:00",ntasks=4,mem=64)
end
=#


st=load(joinpath(@__DIR__, "missedjobs.jld2"))
index=st["index"]
for ja in eachindex(index)
  arguments=Float64.(index[ja])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="120:00",ntasks=4,mem=72)
end



