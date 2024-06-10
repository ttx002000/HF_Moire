using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "finiteS"




#=
spin=args[1]
V0=args[2]
ϕ=args[3]/180*π
Geonum=Int(args[4]);
Nx,Ny,l1,l2=Geometry(Geonum)
scale=args[5];
constq=args[6]/(Nx*Ny)
trytimes=Int(args[7])
=#

spinspace=collect(1.0:0.5:10.0)
constq=2.0
scale=1.0
Geonum=8
for spin in spinspace,jtry in 1:3
   arguments=[spin,0.0,0.0,Geomnum,scale,constq,jtry]
   submit_job(filepath, @__DIR__, job_prefix,arguments; time="5:00:00",cpus_per_task=48)

end




