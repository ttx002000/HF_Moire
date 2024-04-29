using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "shift"

flux=1.0 #multiply this by pi
ϕ=60.0 #convert this to radian 
Nq=6.0; 
scale=1.0;
constqspace=[0.4,0.6,0.8,1.0] #divide this by Nq^2
V0space=[0.0]
ζ=1.0
shiftspace=[1,2,3]
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
  arguments=[flux,V0,ϕ,Nq,scale,constq,ζ,jtry,shift]

    submit_job(filepath, @__DIR__, job_prefix,arguments; time="12:00:00",ntasks=Int(Nq)^2,mem=256)

end

for jtry in 1:3, constq in constqspace, shift in shiftspace, V0 in V0space
  arguments=[flux,V0,ϕ,3.0,scale,constq,ζ,jtry,shift]

    submit_job(filepath, @__DIR__, job_prefix,arguments; time="12:00:00",ntasks=Int(Nq)^2,mem=256)

end



#=
st=load(joinpath(@__DIR__, "missedindex.jld2"))
index=st["index"]
for ja in eachindex(index)
    arguments=Float64.([flux,index[ja][1],ϕ,index[ja][6],scale,index[ja][2],ζ,index[ja][3],index[ja][4],index[ja][5]])
    submit_job(filepath, @__DIR__, job_prefix,arguments; time="15:00:00",ntasks=Int(Nq)^2,mem=256)
end
=#


