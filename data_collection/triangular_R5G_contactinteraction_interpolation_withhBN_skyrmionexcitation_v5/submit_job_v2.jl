
using Dates
function submit_job(filepath, dirpath, job_prefix,args_list; nodes=1, ntasks=1, time="00:120:00", cpus_per_task=1, mem=64, partition="owners,simes")
   

  

    scratch_dir = ENV["SCRATCH"]
    filepos=Int(args_list[1][16]) #ensure that all the jobs are to the same file positions
    outpath=joinpath(scratch_dir, "triangle_R5G_contact_interpolation_withhBN_skyrmionexcitation_v5/data_output$(filepos)/out")
    slurmpath=joinpath(scratch_dir, "triangle_R5G_contact_interpolation_withhBN_skyrmionexcitation_v5/data_output$(filepos)/slurmfiles")
    mkpath(outpath)
    mkpath(slurmpath)
    num_jobs = length(args_list)
    unique_id = Dates.format(now(), "yyyy-mm-dd_THH:MM:SS") 
    param_file = joinpath(slurmpath, "job_parameters_$(unique_id).txt")

    open(param_file, "w") do io
        for args in args_list
            println(io, join(args, " "))  # Store each job's arguments as a line
        end
    end
 
    filestr = """#!/bin/bash
    #SBATCH --job-name=$(job_prefix*"_"*"array")
    #SBATCH --partition=$partition
    #SBATCH --time=$time
    #SBATCH --array=1-$num_jobs%300
    #SBATCH --nodes=$nodes
    #SBATCH --requeue
    #SBATCH --ntasks=$ntasks
    #SBATCH --mem=$(mem)G
    #SBATCH --mail-type=BEGIN,FAIL,END
    #SBATCH --mail-user=ttx2000@stanford.edu
    #SBATCH --output=$outpath/$(job_prefix)_%A_%a_output.txt
    #SBATCH --error=$outpath/$(job_prefix)_%A_%a_error.txt
    #SBATCH --open-mode=append
    #SBATCH --sockets-per-node=2

    # load Julia module
    ml julia/1.10.0

    # multithreading
    export JULIA_NUM_THREADS=$ntasks
    export JULIA_DEPOT_PATH=\$OAK/ttx/julia_depot
     
    PARAMS_FILE=$(param_file)
    PARAMS=\$(sed -n "\${SLURM_ARRAY_TASK_ID}p" \$PARAMS_FILE)
    # run the script
    julia $filepath \$PARAMS """
    
    slurmfile = joinpath(slurmpath, "$(job_prefix)_array_$(unique_id)")
    open(slurmfile, "w") do io
        write(io, filestr)
    end
    run(`sbatch $slurmfile`)
end

