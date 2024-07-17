

function submit_job(filepath, dirpath, job_prefix,args; nodes=1, ntasks=1, time="00:120:00", mem=64, partition="owners,simes")
    outpath = joinpath(dirpath, "out")
    slurmpath = joinpath(dirpath, "slurmfiles")# Why is there a job_prefix semicolon there?
    mkpath(outpath)
    mkpath(slurmpath)


                
    
    

    name = "$(args[1])scale$(args[2])flux$(args[3])Nq$(args[4])V0$(args[5])phi$(args[6])chist$(args[7])chien$(args[8])Cq"

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
    julia $filepath $(args[1]) $(args[2]) $(args[3]) $(args[4]) $(args[5]) $(args[6]) $(args[7]) $(args[8])"""

    open("$slurmpath/$(name).slurm", "w") do io
        write(io, filestr)
    end
    run(`sbatch $(slurmpath)/$(name).slurm`)
end

