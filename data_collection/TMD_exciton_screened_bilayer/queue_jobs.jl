using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2
include(joinpath(@__DIR__,"submit_job.jl"))
filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "bialyerexciton"


mt=0.2552
mb=0.4585
Vt=-10
ϕt=177.6169
Vb=10
ϕb=45.8366
ϵrspace=[5]
Egspace=[75,110]
ϵr=5
θ=2
w=0
holenum=1
Nqspace=[3,4,5,6,7]




#=
mt=parameters[1]
mb=parameters[2]
Vt=parameters[3]
ϕt=parameters[4]/180*π
Vb=parameters[5]
ϕb=parameters[6]/180*π
ϵr=parameters[7]
Eg=parameters[8]
θ=parameters[9]/180*π
w=parameters[10]


holenum=parameters[11]
trytime=parameters[12]
Nq=parameters[13]
=#

for trytime in 1:2, Eg in Egspace,seednum in 1:3, Nq in Nqspace
  arguments=Float64.([mt,mb,Vt,ϕt,Vb,ϕb,ϵr,Eg,θ,w,holenum,trytime,Nq,seednum])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="$(Nq):00:00",ntasks=Nq^2)
end





#=
st=load(joinpath(@__DIR__, "missedjobs.jld2"))
index=st["vector"]
for ja in eachindex(index)
  arguments=Float64.([mt,mb,Vt,ϕt,Vb,ϕb,ϵr,index[ja][1],θ,w,holenum,index[ja][2],Nq,index[ja][3]])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="4:00:00",ntasks=9)
end
=#