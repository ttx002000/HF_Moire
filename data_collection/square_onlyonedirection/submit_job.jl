

function submit_job(filepath, dirpath, job_prefix,args,trytimes; nodes=1, ntasks=1, time="00:120:00", cpus_per_task=1, mem=64, partition="owners,simes")
    outpath = joinpath(dirpath, "out")
    slurmpath = joinpath(dirpath, "slurmfiles")# Why is there a job_prefix semicolon there?
    mkpath(outpath)
    mkpath(slurmpath)

    flux=args[1]
    Vx=args[2]
    Vy=args[3]
    Nq=args[4];
    scale=args[5];
    constq=args[6];
 
    

    name = "$(args[1])flux$(args[2])Vx$(args[3])Vy$(args[4])Nq$(args[5])scale$(args[6])constq$(trytimes)try"

    filestr = """#!/bin/bash
    #SBATCH --job-name=$(job_prefix*"_"*name)
    #SBATCH --partition=$partition
    #SBATCH --time=$time
    #SBATCH --nodes=$nodes
    #SBATCH --ntasks=$ntasks
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
    export JULIA_NUM_THREADS=$(ntasks)

    # run the script
    julia $filepath $(args[1]) $(args[2]) $(args[3]) $(args[4]) $(args[5]) $(args[6]) $(Float64(trytimes))"""

    open("$slurmpath/$(name).slurm", "w") do io
        write(io, filestr)
    end
    run(`sbatch $(slurmpath)/$(name).slurm`)
end

