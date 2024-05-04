using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "ED"

fluxspace=[0.0] #multiply this by pi
ϕ=0.0 #convert this to radian 
Nqspace=[3.0,6.0,9.0]; 
scale=1.0*√2;
constq=2.0 #divide this by Nq^2
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

for ja in eachindex(fluxspace), jNq in eachindex(Nqspace)
  arguments=[fluxspace[ja],V0,ϕ,Nx,Ny,scale,constq]
for jb in 1:2
    submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="$(Int(Nqspace[jNq])+2):00:00",ntasks=Int(Nqspace[jNq])^2,mem=256)
end
end

