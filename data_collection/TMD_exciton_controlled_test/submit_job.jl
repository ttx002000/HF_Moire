

function submit_job(filepath, dirpath, job_prefix,args; nodes=1, ntasks=1, time="00:120:00", cpus_per_task=1, mem=64, partition="owners,simes")
    outpath = joinpath(dirpath, "out")
    slurmpath = joinpath(dirpath, "slurmfiles")# Why is there a job_prefix semicolon there?
    mkpath(outpath)
    mkpath(slurmpath)

 
    
    name = "$(args[1])mt$(args[2])mm$(args[3])mb$(args[4])Vt$(args[5])phit$(args[6])Vm$(args[7])phim$(args[8])Vb$(args[9])phib$(args[10])er$(args[11])Eg$(args[12])theta$(args[13])w$(args[14])holenum$(args[15])trytime$(args[16])Nq$(args[17])seed"
    filestr = """#!/bin/bash
    #SBATCH --job-name=$(job_prefix*"_"*name)
    #SBATCH --partition=$partition
    #SBATCH --time=$time
    #SBATCH --nodes=$nodes
    #SBATCH --ntasks=5
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
    export JULIA_NUM_THREADS=5

    # run the script
    julia $filepath $(args[1]) $(args[2]) $(args[3]) $(args[4]) $(args[5]) $(args[6]) $(args[7]) $(args[8]) $(args[9]) $(args[10]) $(args[11]) $(args[12]) $(args[13]) $(args[14]) $(args[15]) $(args[16]) $(args[17])"""

    open("$slurmpath/$(name).slurm", "w") do io
        write(io, filestr)
    end
    run(`sbatch $(slurmpath)/$(name).slurm`)
end

