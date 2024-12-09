using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2
include(joinpath(@__DIR__,"submit_job.jl"))
filepath = joinpath(@__DIR__, "main_test.jl")

job_prefix = "SSH_withshear"

Nx=20
Ny=20
tper=0.37
tpa=2
tNNN=0.16
#α=-0.05
#β=-0.05
#K=1.0;
KNNN=0.0;
#filling=1.35;
#gshear=0.0
temp=10^(-4)



#=
Nx=Int(args[1])
Ny=Int(args[2])
tper=args[3]
tpa=args[4]
tNNN=args[5]
sum1=args[6]
pol1=args[7]
α=sum1*pol1/(1+pol1)
β=sum1/(1+pol1)

sum2=args[8]
pol2=args[9]
K=sum2/(1+pol2)
gshear=sum2*pol2/(1+pol2)


filling=(args[10])
KNNN=args[11]
temp=args[12]
trytime=Int(args[13])
Nelec=Int(round(Nx*Ny*filling))
=#


for trytime in 7:8, filling in [1.25], pol1 in collect(0.1:0.1:2.0), pol2 in collect(0.1:0.1:2.0)
 
   sum1=-2.5
   sum2=2.0

  arguments=Float64.([Nx,Ny,tper,tpa,tNNN,sum1,pol1,sum2,pol2,filling,KNNN,temp,trytime])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="40:00",ntasks=1,mem=8)
end



#=
st=load(joinpath(@__DIR__, "missedjobs.jld2"))
index=st["index"]
for ja in eachindex(index)
  arguments=index[ja]
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="3:00:00",ntasks=1,mem=16)
end
=#


