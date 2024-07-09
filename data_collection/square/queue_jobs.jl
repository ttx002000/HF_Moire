using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra


fluxspace=[0.3,0.5] #I multiply it by pi when doing the calculation
V0=0.0
ϕ=0.0 #I convert this degree to randian
Nq=12.0;
scale=1.0;
constqspace=[0.5,1.0,1.5,2.0] #I divide it by Nq^2 in the actual calculation


include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "square_HF"


for jtry in 1:3,constq in constqspace,flux in fluxspace
     arguments=[flux,V0,ϕ,Nq,scale,constq,jtry]
     submit_job(filepath, @__DIR__, job_prefix,arguments; time="4:00:00",ntasks=32,mem=64)
end