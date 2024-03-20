

function submit_job(filepath, dirpath, job_prefix,arguments; nodes=1, ntasks=1, time="00:120:00", cpus_per_task=1, mem=64, partition="owners,simes")
    outpath = joinpath(dirpath, "out")
    slurmpath = joinpath(dirpath, "slurmfiles")# Why is there a job_prefix semicolon there?
    mkpath(outpath)
    mkpath(slurmpath)

    arguments=[vf,spin,V0,ϕ]
    vf=arguments[1]
    V0=arguments[3]
    ϕ=arguments[4]
    spin=arguments[2];
    
 

    name = "$(vf)vf$(spin)spin$(V0)V0$(ϕ)phi"
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
    julia $filepath $vf $spin $V0 $ϕ"""

    open("$slurmpath/$(name).slurm", "w") do io
        write(io, filestr)
    end
    run(`sbatch $(slurmpath)/$(name).slurm`)
end

