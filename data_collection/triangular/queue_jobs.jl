using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "piflux"

fluxspace=[1.0] #multiply this by pi
ϕ=0.0 #convert this to radian 
Nqspace=[9.0]; 
scale=1.0;
constqspace=collect(0.2:0.2:2.0) #divide this by Nq^2
V0space=[0.0]


#=
for ja in eachindex(fluxspace), jNq in eachindex(Nqspace), constq in constqspace,V0 in V0space
  arguments=[fluxspace[ja],V0,ϕ,Nqspace[jNq],scale,constq]
for jb in 1:3
    submit_job(filepath, @__DIR__, job_prefix,arguments,jb; time="$(Int(Nqspace[jNq])):00:00",ntasks=16,mem=32)
end
end
=#
V0space=collect(0.2:0.2:3.0)
constqspace=[0.0]

for ja in eachindex(fluxspace), jNq in eachindex(Nqspace), constq in constqspace,V0 in V0space
  arguments=[fluxspace[ja],V0,ϕ,Nqspace[jNq],scale,constq]
for jb in 1:2
    submit_job(filepath, @__DIR__, job_prefix,arguments,jb; time="40:00",ntasks=10,mem=32)
end
end


#=
index=load(joinpath(@__DIR__, "missedjobs1.jld2"))["index"]
for ja in eachindex(index), jb in eachindex(fluxspace), jNq in eachindex(Nqspace)

  arguments=[0.7,index[ja][1],index[ja][2],Nqspace[jNq],scale,index[ja][3]]
  submit_job(filepath, @__DIR__, job_prefix,arguments,Int(index[ja][4]); time="10:00:00",ntasks=32,mem=32)

end

index=load(joinpath(@__DIR__, "missedjobs2.jld2"))["index"]
for ja in eachindex(index), jb in eachindex(fluxspace), jNq in eachindex(Nqspace)

  arguments=[1.0,index[ja][1],index[ja][2],Nqspace[jNq],scale,index[ja][3]]
  submit_job(filepath, @__DIR__, job_prefix,arguments,Int(index[ja][4]); time="10:00:00",ntasks=32,mem=32)

end
=#