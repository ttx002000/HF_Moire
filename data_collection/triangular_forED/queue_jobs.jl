using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "ED"

flux=2.0 #multiply this by pi
ϕ=0.0 #convert this to radian 
Nx=8.0
Ny=4.0
scale=1.0;
constq=6.0 #divide this by Nq^2
V0=0.0


#=
flux=args[1]*π
V0=args[2]
ϕ=args[3]/180*π
Nx=Int(args[4]);
Ny=Int(args[5]);
scale=args[6];
constq=args[7]/Nq^2
trytimes=Int(args[8])
=#



arguments=Float64.([flux,V0,ϕ,Nx,Ny,scale,constq])
 for jb in 1:3
   submit_job(filepath, @__DIR__, job_prefix,arguments,jb; time="3:00:00",ntasks=16,mem=128)
end

arguments=Float64.([flux,V0,ϕ,Nx,Ny,scale,5.0])
 for jb in 1:3
   submit_job(filepath, @__DIR__, job_prefix,arguments,jb; time="3:00:00",ntasks=16,mem=128)
end