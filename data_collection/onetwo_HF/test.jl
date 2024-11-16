using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2
s=2
scratch_dir = ENV["SCRATCH"]
output_path = joinpath(scratch_dir, "my_output_file.txt")
jldsave(output_path,s=s)