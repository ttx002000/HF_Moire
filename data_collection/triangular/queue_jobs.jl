using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "piflux"

fluxspace=[0.0] #multiply this by pi
ϕ=0.0 #convert this to radian 
Nqspace=[3.0,4.0,5.0,6.0,7.0,8.0,9.0]; 
scale=1.0*√2;
constq=2.0 #divide this by Nq^2
V0=0.0

for ja in eachindex(fluxspace), jNq in eachindex(Nqspace)
  arguments=[fluxspace[ja],V0,ϕ,Nqspace[jNq],scale,constq]
for ja in 1:2
    submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="$(Int(Nqspace[jNq])+2):00:00",ntasks=Int(Nqspace[jNq])^2,mem=256)
end
end

#=
fluxspace=[0.7] #multiply this by pi
ϕ=0.0 #convert this to radian 
Nqspace=[3.0,4.0,5.0,6.0,7.0,8.0,9.0]; 
scale=1.0*√0.7;
constq=1.5*0.7 #divide this by Nq^2
V0=0.0

for ja in eachindex(fluxspace), jNq in eachindex(Nqspace)
  arguments=[fluxspace[ja],V0,ϕ,Nqspace[jNq],scale,constq]
for ja in 1:2
    submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="$(Int(Nqspace[jNq])+2):00:00",ntasks=Int(Nqspace[jNq])^2,mem=256)
end
end

fluxspace=[1.3] #multiply this by pi
ϕ=0.0 #convert this to radian 
Nqspace=[3.0,4.0,5.0,6.0,7.0,8.0,9.0]; 
scale=1.0*√1.3;
constq=1.5*1.3 #divide this by Nq^2
V0=0.0

for ja in eachindex(fluxspace), jNq in eachindex(Nqspace)
  arguments=[fluxspace[ja],V0,ϕ,Nqspace[jNq],scale,constq]
for ja in 1:2
    submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="$(Int(Nqspace[jNq])+2):00:00",ntasks=Int(Nqspace[jNq])^2,mem=256)
end
end

fluxspace=[1.2] #multiply this by pi
ϕ=0.0 #convert this to radian 
Nqspace=[3.0,4.0,5.0,6.0,7.0,8.0,9.0]; 
scale=1.0*√1.2;
constq=1.5*1.2 #divide this by Nq^2
V0=0.0

for ja in eachindex(fluxspace), jNq in eachindex(Nqspace)
  arguments=[fluxspace[ja],V0,ϕ,Nqspace[jNq],scale,constq]
for ja in 1:2
    submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="$(Int(Nqspace[jNq])+2):00:00",ntasks=Int(Nqspace[jNq])^2,mem=256)
end
end
=#