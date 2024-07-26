using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2
include(joinpath(@__DIR__,"submit_job.jl"))
filepath = joinpath(@__DIR__, "main_test.jl")

job_prefix = "exciton"


mtspace=collect(-0.6:0.1:-0.4)
mm=0.62
mb=0.62
Vt=0.0
ϕt=0.0
Vm=-11.2
ϕm=91.0
Vb=-11.2
ϕb=-91.0
ϵrspace=[10.0]
Egspace=collect(100.0:5.0:160.0)
period=10.84
w=-13.3
elecnum=1.0

geonum=3
cutoffnum=4.01







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
period=parameters[12]
w=parameters[13]

elecnum=parameters[14]
trytime=parameters[15]
geonum=parameters[16]
cutoffnum=parameters[17]
=#



for  Eg in Egspace,trytime in 1:10,ϵr in ϵrspace, mt in mtspace
  arguments=Float64.([mt,mm,mb,Vt,ϕt,Vm,ϕm,Vb,ϕb,ϵr,Eg,period,w,elecnum,trytime,geonum,cutoffnum])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="8:00:00",ntasks=8,mem=16)
end




