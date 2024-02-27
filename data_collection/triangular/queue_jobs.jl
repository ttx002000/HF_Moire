using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra

include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "HF"

flux=1.0 #multiply this by pi
V0=3.0
ϕ=0.0 #convert this to radian 
Nq=10.0; 
scale=1.0;
constq=0.1 #divide this by Nq^2


arguments=[flux,V0,ϕ,Nq,scale,constq]
for ja in 1:3
     submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="500:00",cpus_per_task=50)
end


flux=1.0 #multiply this by pi
V0=3.0
ϕ=0.0 #convert this to radian 
Nq=10.0; 
scale=1.0;
constq=1.5 #divide this by Nq^2

arguments=[flux,V0,ϕ,Nq,scale,constq]
for ja in 1:3
     submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="500:00",cpus_per_task=50)
end

flux=1.0 #multiply this by pi
V0=3.0
ϕ=0.0 #convert this to radian 
Nq=10.0; 
scale=1.0;
constq=2.5 #divide this by Nq^2

arguments=[flux,V0,ϕ,Nq,scale,constq]
for ja in 1:3
     submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="500:00",cpus_per_task=50)
end




flux=1.0 #multiply this by pi
V0=0.0
ϕ=0.0 #convert this to radian 
Nq=10.0; 
scale=1.0;
constq=1.0 #divide this by Nq^2

arguments=[flux,V0,ϕ,Nq,scale,constq]
for ja in 1:3
     submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="500:00",cpus_per_task=50)
end

flux=1.0 #multiply this by pi
V0=4.0
ϕ=0.0 #convert this to radian 
Nq=10.0; 
scale=1.0;
constq=0.1 #divide this by Nq^2


arguments=[flux,V0,ϕ,Nq,scale,constq]
for ja in 1:3
     submit_job(filepath, @__DIR__, job_prefix,arguments,ja; time="500:00",cpus_per_task=50)
end