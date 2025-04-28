

function submit_job(filepath, dirpath, job_prefix,args; nodes=1, ntasks=1, time="00:120:00", cpus_per_task=1, mem=64, partition="owners,simes")
    #outpath = joinpath(dirpath, "out")
    #slurmpath = joinpath(dirpath, "slurmfiles")# Why is there a job_prefix semicolon there?
    #mkpath(outpath)
    #mkpath(slurmpath)

    scratch_dir = ENV["SCRATCH"]
    outpath=joinpath(scratch_dir, "bandprojected_RMG/data_output$(Int(args[12]))/out")
    slurmpath=joinpath(scratch_dir, "bandprojected_RMG/data_output$(Int(args[12]))/slurmfiles")
    mkpath(outpath)
    mkpath(slurmpath)
 

    name = "$(args[1])radius$(args[2])num_points$(args[3])uD$(args[4])er$(args[5])NL$(args[6])tgden$(args[7])temp$(args[8])JH$(args[9])SOC$(args[10])whichside$(args[11])trytime$(args[12])filepos"
    filestr = """#!/bin/bash
    #SBATCH --job-name=$(job_prefix*"_"*name)
    #SBATCH --partition=$partition
    #SBATCH --time=$time
    #SBATCH --nodes=$nodes
    #SBATCH --extra-node-info 2-2:*:*
    #SBATCH --ntasks=$ntasks
    #SBATCH --cpus-per-task=$cpus_per_task
    #SBATCH --mem=$(mem)G
    #SBATCH --mail-type=BEGIN,FAIL,END
    #SBATCH --mail-user=ttx2000@stanford.edu
    #SBATCH --output=$outpath/$(job_prefix*"_"*name)_output.txt
    #SBATCH --error=$outpath/$(job_prefix*"_"*name)_error.txt
    #SBATCH --open-mode=append
    #SBATCH --sockets-per-node=2

    # load Julia module
    ml julia/1.10.0

    # multithreading
    export JULIA_NUM_THREADS=\$SLURM_CPUS_ON_NODE
    export JULIA_DEPOT_PATH=\$SCRATCH/julia_depot

    # run the script
    julia  $filepath $(args[1]) $(args[2]) $(args[3]) $(args[4]) $(args[5]) $(args[6]) $(args[7]) $(args[8]) $(args[9]) $(args[10]) $(args[11]) $(args[12])"""

    open("$slurmpath/$(name).slurm", "w") do io
        write(io, filestr)
    end
    run(`sbatch $(slurmpath)/$(name).slurm`)
end

