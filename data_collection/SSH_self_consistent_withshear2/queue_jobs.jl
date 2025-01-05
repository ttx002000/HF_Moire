using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2
include(joinpath(@__DIR__,"submit_job.jl"))
filepath = joinpath(@__DIR__, "main_test.jl")

job_prefix = "SSH_withshear"

Nx=20
Ny=20
tper=-0.37
tpa=2
tNNN=-0.08
#α=-0.05
#β=-0.05
K=1.0;
KNNN=0.0;
filling=1.25;
#gshear=0.0
temp=0.01



#=
Nx=Int(args[1])
Ny=Int(args[2])
tper=args[3]
tpa=args[4]
tNNN=args[5]
α=args[6]
β=args[7]

K=args[8]
gshear=args[9]


filling=(args[10])
KNNN=args[11]
temp=args[12]
trytime=Int(args[13])
Nelec=Int(round(Nx*Ny*filling))
=#




#=
for trytime in 1:10, gshear in [0.1,0.5,1.0,2.0], α in collect(-1.0:-0.1:-1.5)
 
  β=-α

  arguments=Float64.([Nx,Ny,tper,tpa,tNNN,α,β,K,gshear,filling,KNNN,temp,trytime])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="40:00",ntasks=1,mem=8)
end
=#



st=load(joinpath(@__DIR__, "missedjobs.jld2"))
index=st["index"]
for ja in eachindex(index)
  arguments=index[ja]
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="3:00:00",ntasks=1,mem=8)
end




