using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "onetwo"

wAA=75.0
wAB=110.0
vF=579.2265 #unit meV*nm

lambda_MDT=-0.2
Nb_up=3
Nb_down=3
θ=1.23
er=8.0
Nband=Nb_up+Nb_down
geonum=2
filling=3*4+7


#=
wAA=75.0
wAB=110.0
vF=579.2265 #unit meV*nm

lambda_MDT=-0.2
Nb_up=3
Nb_down=3
θ=1.23/180*π
ϵr=8.0
Nband=Nb_up+Nb_down
geonum=Int(1)
filling=3*4+8
=#


for  jtry in 61:90, er in [8.0], Dfield in [0.0,10.0,20.0], filling in [19,20,21]
  arguments=Float64.([wAA,wAB,vF,lambda_MDT,Nb_up,Nb_down,θ,er,geonum,filling,Dfield,jtry])
    submit_job(filepath, @__DIR__, job_prefix,arguments; time="5:00:00",ntasks=16,mem=64)
end



#=
st=load(joinpath(@__DIR__, "missedjobs.jld2"))
index=st["index"]
for ja in eachindex(index)
  arguments=Float64.(index[ja])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="10:00:00",ntasks=16,mem=64)
end
=#