using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "piflux"

V0=0.0
ϕ=0.0 #convert this to radian 
Nq=9.0; 
periodspace=[10.0:1.0:20.0]
erspace=collect(3.0:1.0:8.0) #divide this by Nq^2

#=
V0=args[1]
ϕ=args[2]/180*π
Nq=Int(args[3]);
period=args[4]
ϵr=args[5]
Dfield=args[6]
trytime=args[7]
=#


for  period in periodspace, ϵr in erspace,jtry in 1:3
  arguments=Float64.([V0,ϕ,Nq,period,ϵr,jtry])
    submit_job(filepath, @__DIR__, job_prefix,arguments; time="$(Int(Nqspace[jNq])):00:00",ntasks=16,mem=32)
end


