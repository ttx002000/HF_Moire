using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2
include(joinpath(@__DIR__,"submit_job.jl"))
filepath = joinpath(@__DIR__, "main_test.jl")

job_prefix = "exciton"

mt=0.36
mm=-0.6
mb=-0.6
Vt=0.0
ϕt=0.0
Vm=20.8
ϕm=107.7
Vb=20.8
ϕb=-107.7

#ϵr=6
#Eg=70.0
θ=3.89
w=-23.8
holenum=1.0
geonum=1.0

#args=[0.36,-0.6,-0.6,0.0,0.0,20.8,107.7,20.8,-107.7,6,70.0,3.89,-23.8,1.0,1.0,3.0]



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
geonum=parameters[16]
=#


for trytime in 1:10,ϵr in [15], Eg in [140,160,180]
  arguments=Float64.([mt,mm,mb,Vt,ϕt,Vm,ϕm,Vb,ϕb,ϵr,Eg,θ,w,holenum,trytime,geonum])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="4:00:00",ntasks=4,mem=8)
end

for trytime in 1:10,ϵr in [18], Eg in [160,180,200]
  arguments=Float64.([mt,mm,mb,Vt,ϕt,Vm,ϕm,Vb,ϕb,ϵr,Eg,θ,w,holenum,trytime,geonum])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="4:00:00",ntasks=4,mem=8)
end



#=
st=load(joinpath(@__DIR__, "missedjobs.jld2"))
index=st["vector"]
for ja in eachindex(index)
  arguments=Float64.([mt,mm,mb,Vt,ϕt,Vm,ϕm,Vb,ϕb,ϵr,index[ja][1],θ,w,holenum,1.0,Nq,index[ja][2]])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="12:00:00",ntasks=9)
end
=#