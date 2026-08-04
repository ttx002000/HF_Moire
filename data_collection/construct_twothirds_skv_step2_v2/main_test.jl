using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

using JLD2

include("../../src/construct_twothirds_skv_step2_v2.jl")

args = parse.(Float64, ARGS)

results = main_func(args)

scratch_dir = ENV["SCRATCH"]
file_pos = Int(args[14])

filename =
    "Gq_" *
    "$(args[1])f1$(args[2])f2$(args[3])NL$(args[4])am" *
    "$(args[5])N1$(args[6])N2$(args[7])N1f$(args[8])N2f" *
    "$(args[9])xir$(args[10])xii$(args[11])grid" *
    "$(args[12])type$(args[13])gaud.jld2"

savepath = joinpath(
    scratch_dir,
    "constrcut_twothirds_skv_v2/data_output$(file_pos)/structure_factor",
    filename,
)

mkpath(dirname(savepath))

jldsave(savepath; results...)

println("Saved to:")
println(savepath)

