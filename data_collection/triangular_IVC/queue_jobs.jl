using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "shift"

flux=0.0 #multiply this by pi
ϕ=60.0 #convert this to radian 
Nq=3.0; 
scale=1.0;
constqspace=[0.1,0.2,0.3,0.4] #divide this by Nq^2
V0space=[0.1,0.2,0.3,0.4]
ζ=1.0
shiftspace=[[1.0,1.0],[6.0,8.0],[8.0,6.0]]
#=
    flux=arg[1]
    V0=arg[2]
    ϕ=arg[3]
    Nq=arg[4];
    scale=arg[5];
    constq=arg[6];
    ζ=arg[7];
    trytimes=arg[8]
    shiftindex=arg[9]
    =#

for jtry in 1:3, constq in constqspace, shift in shiftspace, V0 in V0space
  arguments=[flux,V0,ϕ,Nq,scale,constq,ζ,jtry,shift[1],shift[2]]

    submit_job(filepath, @__DIR__, job_prefix,arguments; time="8:00:00",ntasks=Int(Nq)^2,mem=256)

end



#=
st=load(joinpath(@__DIR__, "missedjobs.jld2"))
index=st["index"]
for ja in eachindex(index)
    arguments=[index[ja][1],V0,ϕ,Nq,scale,constq,ζ]
    trytime=Int(index[ja][2])
    submit_job(filepath, @__DIR__, job_prefix,arguments,trytime; time="38:00:00",cpus_per_task=36)
end
=#