using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "pentalayer"


erspace=[4.0,6.0,8.0]
θspace=[0.88,0.9]
Nband=13
uDspace=[30.0,40.0,50.0,60.0,70.0]
Nq=15

#=
ϵr=args[1]
θ=args[2]/180*π
Nq=Int(args[3]);
uD=args[4]
Nband=Int(args[5])
trytime=Int(args[6])
=#
#=
for  jtry in 1:15, er in erspace
  arguments=Float64.([er,0.0,Nq,20.0,Nband,jtry])
    submit_job(filepath, @__DIR__, job_prefix,arguments; time="5:00:00",ntasks=16,mem=64)
end
=#

for  jtry in 1:15, er in erspace, θ in θspace, uD in uDspace
  arguments=Float64.([er,θ,Nq,uD,Nband,jtry])
    submit_job(filepath, @__DIR__, job_prefix,arguments; time="6:00:00",ntasks=32,mem=128)
end



#=
st=load(joinpath(@__DIR__, "missedjobs.jld2"))
index=st["index"]
for ja in eachindex(index)
  arguments=Float64.(index[ja])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="10:00:00",ntasks=16,mem=64)
end
=#