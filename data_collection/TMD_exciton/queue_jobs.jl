using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "exciton"

#parameters=[0.35,0.4,0.35,-10.0,70.0,10.0,80.0,10.0,1.0,10,100.0,2.0,5.0]
mt=0.35
mm=0.4
mb=0.35
Vt=-10
ϕt=70
Vm=10
ϕm=80
Vb=10
ϕb=1
ϵr=7
Egspace=[110,120,130]
θ=2
w=0
holenum=2
Nq=3



#=
mt=parameters[1] actually args not parameters
mm=parameters[2]
mb=parameters[3]
Vt=parameters[4]
ϕt=parameters[5]/180*π
Vm=parameters[6]
ϕm=parameters[7]/180*π
Vb=parameters[8]
ϕb=parameters[9]/180*π
ϵr=parameters[10]
Eg=parameters[11]
θ=parameters[12]/180*π
w=parameters[13]

holenum=parameters[14]
trytime=parameters[15]
Nq=parameters[16]
=#

for trytime in 1:15, Eg in Egspace
  arguments=Float64.([mt,mm,mb,Vt,ϕt,Vm,ϕm,Vb,ϕb,ϵr,Eg,θ,w,holenum,trytime,Nq])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="5:00:00",cpus_per_task=1)
end





#=
st=load(joinpath(@__DIR__, "missedjobs.jld2"))
index=st["index"]
for ja in eachindex(index)
    arguments=[index[ja][1],V0,ϕ,Nq,scale,constq]
    trytime=Int(index[ja][2])
    submit_job(filepath, @__DIR__, job_prefix,arguments,trytime; time="38:00:00",cpus_per_task=36)
end
=#