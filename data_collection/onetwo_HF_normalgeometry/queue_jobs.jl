using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "onetwo_normal"

wAA=75.0
wAB=110.0
vF=579.2265 #unit meV*nm

lambda_MDT=-0.2
Nb_up=1
Nb_down=1
θ=1.23
er=8.0
Nband=Nb_up+Nb_down
geonum=3
filling=6
#shift=1


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


for  jtry in 1:100, er in [8.0], Dfield in collect(0.0:10.0:20.0),shift in [2.0]
  arguments=Float64.([wAA,wAB,vF,lambda_MDT,Nb_up,Nb_down,θ,er,geonum,filling,Dfield,shift,jtry])
    submit_job(filepath, @__DIR__, job_prefix,arguments; time="1:30:00",ntasks=16,mem=32)
end









#=
st=load(joinpath(@__DIR__, "missedjobs.jld2"))
index=st["index"]
for ja in eachindex(index)
  arguments=Float64.(index[ja])
  submit_job(filepath, @__DIR__, job_prefix,arguments; time="3:00:00",ntasks=16,mem=32)
  #=
  if arguments[5]==4.0
    submit_job(filepath, @__DIR__, job_prefix,arguments; time="3:30:00",ntasks=16,mem=32)
  else
    submit_job(filepath, @__DIR__, job_prefix,arguments; time="2:30:00",ntasks=16,mem=32)
  end
  =#

end
=#

