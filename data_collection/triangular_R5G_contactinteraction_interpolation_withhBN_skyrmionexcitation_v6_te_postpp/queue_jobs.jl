using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using LinearAlgebra
using JLD2

include("submit_job_v2.jl")


aa= parse.(Int, ARGS)


filepath = joinpath(@__DIR__, "main_test.jl")
job_prefix = "time_$(Int(aa[1]))"

scratch_dir = ENV["SCRATCH"]
misspath=joinpath(scratch_dir, "triangle_R5G_contact_interpolation_withhBN_skyrmionexcitation_v6_time_evolution/data_output$(Int(aa[1]))/Wilsondata/missedjobs.jld2")

st=load(misspath)

index=st["index"]
start=1
ee=length(index)
ba_size = 300
count = 1

while (count - 1) * ba_size + 1 <= ee
    batch_start = (count - 1) * ba_size + 1
    batch_end = min(count * ba_size, ee)

    submit_job(
        filepath,
        @__DIR__,
        job_prefix,
        index[batch_start:batch_end];
        time = "12:00:00",
        cpus_per_task = 16,
        mem = 32
    )

    println(batch_start, " ", batch_end)

    sleep(5)
    global count += 1
end
