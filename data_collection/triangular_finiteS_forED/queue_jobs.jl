using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "ED_finiteS"




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

spinspace=collect(10.5:0.5:13.0)
constq=2.0
scale=1.0
Geonum=10
for spin in spinspace,jtry in 1:3
   arguments=[spin,0.0,0.0,Geonum,scale,constq,jtry]
   submit_job(filepath, @__DIR__, job_prefix,arguments; time="8:00:00",ntasks=32,mem=128)

end

Geonum=8
for spin in spinspace,jtry in 1:3
   arguments=[spin,0.0,0.0,Geonum,scale,constq,jtry]
   submit_job(filepath, @__DIR__, job_prefix,arguments; time="8:00:00",ntasks=48,mem=128)

end

Geonum=9
for spin in spinspace,jtry in 1:3
   arguments=[spin,0.0,0.0,Geonum,scale,constq,jtry]
   submit_job(filepath, @__DIR__, job_prefix,arguments; time="8:00:00",ntasks=48,mem=128)

end




