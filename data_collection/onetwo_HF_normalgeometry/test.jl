using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2,Plots
s=randn(5,5)
jldsave("test.jld2",
       s=s)
