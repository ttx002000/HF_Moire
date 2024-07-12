using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra


fluxspace=collect(0.0:0.05:2.0) #I multiply it by pi when doing the calculation
Vspace=[2.0,10.0]
ϕ=0.0 #I convert this degree to randian
Nq=30.0;
scale=1.0;
maxg=9.01


include("submit_job.jl")

filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "square_singleparticle"


for flux in fluxspace, V0 in Vspace
     arguments=[flux,V0,ϕ,Nq,scale,maxg,jtry]
     submit_job(filepath, @__DIR__, job_prefix,arguments; time="30:00",ntasks=16,mem=16)
end