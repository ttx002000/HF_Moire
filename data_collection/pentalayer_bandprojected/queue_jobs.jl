using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "pentalayer"


erspace=[4.0,5.0,6.0]
θspace=[0.77]
Nband=12
uDspace=[50.0]
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

for  jtry in 1:5, er in erspace, θ in θspace, uD in uDspace
  arguments=Float64.([er,θ,Nq,uD,Nband,jtry])
    submit_job(filepath, @__DIR__, job_prefix,arguments; time="6:00:00",ntasks=16,mem=64)
end

for  jtry in 1:5, er in erspace, θ in θspace, uD in uDspace
  arguments=Float64.([er,θ,Nq,uD,7,jtry])
    submit_job(filepath, @__DIR__, job_prefix,arguments; time="6:00:00",ntasks=16,mem=64)
end

for  jtry in 1:5, er in erspace, θ in θspace, uD in uDspace
  arguments=Float64.([er,θ,Nq,uD,6,jtry])
    submit_job(filepath, @__DIR__, job_prefix,arguments; time="6:00:00",ntasks=16,mem=64)
end

for  jtry in 1:5, er in erspace, θ in θspace, uD in uDspace
  arguments=Float64.([er,θ,Nq,uD,8,jtry])
    submit_job(filepath, @__DIR__, job_prefix,arguments; time="6:00:00",ntasks=16,mem=64)
end


#=
st=load(joinpath(@__DIR__, "missedjobs.jld2"))
index=st["index"]
for ja in eachindex(index)
  arguments=Float64.(index[ja])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="10:00:00",ntasks=16,mem=64)
end
=#