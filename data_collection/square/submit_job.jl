

function submit_job(filepath, dirpath, job_prefix,arguments,trytimes; nodes=1, ntasks=1, time="00:120:00", cpus_per_task=1, mem=64, partition="owners,simes")
    outpath = joinpath(dirpath, "out")
    slurmpath = joinpath(dirpath, "slurmfiles")# Why is there a job_prefix semicolon there?
    mkpath(outpath)
    mkpath(slurmpath)

    flux=arguments[1]
    V0=arguments[2]
    ϕ=arguments[3]
    Nq=arguments[4];
    scale=arguments[5];
    constq=arguments[6];
 
    

    name = "$(Nq)Nq$(flux)flux$(V0)V0$(ϕ)ϕ$(scale)scale$(constq)constq$(trytimes)try"

    filestr = """#!/bin/bash
    #SBATCH --job-name=$(job_prefix*"_"*name)
    #SBATCH --partition=$partition
    #SBATCH --time=$time
    #SBATCH --nodes=$nodes
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

    # run the script
    julia $filepath $flux $V0 $ϕ $Nq $scale $constq $(Float64(trytimes))"""

    open("$slurmpath/$(name).slurm", "w") do io
        write(io, filestr)
    end
    run(`sbatch $(slurmpath)/$(name).slurm`)
end

