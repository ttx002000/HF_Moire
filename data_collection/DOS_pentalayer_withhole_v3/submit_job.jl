

function submit_job(filepath, dirpath, job_prefix,args; nodes=1, ntasks=1, time="00:120:00", cpus_per_task=1, mem=64, partition="owners,simes")

    
    scratch_dir = ENV["SCRATCH"]
    outpath=joinpath(scratch_dir, "pentalayer_DOS_v3/data_output$(Int(args[11]))/out")
    slurmpath=joinpath(scratch_dir, "pentalayer_DOS_v3/data_output$(Int(args[11]))/slurmfiles")
    mkpath(outpath)
    mkpath(slurmpath)
 
 
    
    name = "$(args[1])uD$(args[2])sample$(args[3])nstart$(args[4])nend$(args[5])perturb$(args[6])DosEbin$(args[7])Eint$(args[8])kradius$(args[9])Elow$(args[10])Ehigh$(args[11])filepos"
    filestr = """#!/bin/bash
    #SBATCH --job-name=$(job_prefix*"_"*name)
    #SBATCH --partition=$partition
    #SBATCH --time=$time
    #SBATCH --nodes=$nodes
    #SBATCH --requeue
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
    export JULIA_NUM_THREADS=$ntasks

    # run the script
    julia $filepath $(args[1]) $(args[2]) $(args[3]) $(args[4]) $(args[5]) $(args[6]) $(args[7]) $(args[8]) """

    open("$slurmpath/$(name).slurm", "w") do io
        write(io, filestr)
    end
    run(`sbatch $(slurmpath)/$(name).slurm`)
end

