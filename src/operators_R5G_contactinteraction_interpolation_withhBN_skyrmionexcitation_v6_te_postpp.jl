using LinearAlgebra
using Arpack
using Combinatorics
using Random
using JLD2
using StaticArrays
using Plots
BLAS.set_num_threads(1)


function main(args)
    # -------------------------
    # Fixed model parameters
    # -------------------------
    NL = Int(args[1])
    theta=args[2]
    constq = args[3]
    ϵr = args[4]
    uD = args[5]
    filling = Int(args[6])
    λ = args[8]
    which_try=Int(args[9])
    enlarge_factor = Int(args[10])

    V0_hBN = args[11]
    V1_hBN = args[12]
    ψ_hBN = args[13]
    V2_scalar = args[14]
    ϕ = args[15] / 180 * π

    pin_coeff = args[16]
    dedis = args[17]
    defec_pos = Int(args[18])
    filepos=Int(args[19])

    # -------------------------
    # TDHF parameters
    # -------------------------
    total_steps = Int(args[20])
    save_every = Int(args[21])
    refresh_every = Int(args[22])
    plot_every = Int(args[23])

    dt = args[24]

    electric_field_magnitude = args[25]
    electric_field_angle_degree = args[26]

    gamma = args[27]
    temp = args[28]
    gcutoff_work = args[29]

    start_step = Int(args[30])
    direction=[Int(args[31]),Int(args[32])]


    ckp_list=collect(start_step:save_every:total_steps)

    scratch_dir = ENV["SCRATCH"]

    base_folder=joinpath(
        scratch_dir,
        "triangle_R5G_contact_interpolation_withhBN_skyrmionexcitation_v6_time_evolution",
        "data_output$(filepos)"
    )


    seed_folder=joinpath(
        scratch_dir,
        "triangle_R5G_contact_interpolation_withhBN_skyrmionexcitation_v6_time_evolution",
        "data_output$(filepos)",
        "seed"
    )

       
    seed_path=joinpath(seed_folder,"$(args[1])NL$(args[2])theta$(args[3])constq$(args[4])ϵr$(args[5])uD$(args[6])filling$(args[7])cutoff$(args[8])lambda$(args[9])trytime$(args[10])enlarge$(args[11])V0_hBN$(args[12])V1_hBN$(args[13])ψ_hBN$(args[14])V2_scalar$(args[15])ϕ$(args[16])pincof$(args[17])dedis$(args[18])depos.jld2")
     
   if isfile(seed_path)
    st_seed=load(joinpath("seed",seed_path))
   else
    error("no seed file")
   end

    T1=st_seed["T1"]
    T2=st_seed["T2"]
    a1m=st_seed["a1m"]
    a2m=st_seed["a2m"]
    Area=abs(a1m[2]*a2m[1]-a1m[1]*a2m[2])

    Z_record = fill(NaN + im * NaN, length(ckp_list))
    missing_record = falses(length(ckp_list))
    error_record = fill("", length(ckp_list))



     Threads.@threads :greedy for kkk in eachindex(ckp_list)
                try
                  ckpoint = ckp_list[kkk]

             

                checkpoint_path =
                    joinpath(basefolder,"$(args[1])NL$(args[2])theta$(args[3])constq$(args[4])ϵr$(args[5])uD$(args[6])filling$(args[7])cutoff$(args[8])lambda$(args[9])trytime$(args[10])enlarge$(args[11])V0_hBN$(args[12])V1_hBN$(args[13])ψ_hBN$(args[14])V2_scalar$(args[15])ϕ$(args[16])pincof$(args[17])dedis$(args[18])depos" *
                    "$(args[22])refreshevery" *
                    "$(args[24])dt$(args[25])Emag$(args[26])Eag$(args[27])gamma$(args[28])temp$(args[29])workcutoff" *
                    "$(ckpoint)stepnum.jld2")

                    if !isfile(checkpoint_path)
                        missing_record[kkk] = true
                        continue
                    end

                wave_work = nothing
                Ashift = nothing
                Prj = nothing
                spinor_set = nothing

                jldopen(checkpoint_path, "r") do checkpoint_file
                    wave_work = checkpoint_file["wave_work"]
                    Ashift = checkpoint_file["Ashift"]
                    Prj = checkpoint_file["Prj"]

                    if haskey(checkpoint_file, "spinor_set")
                        spinor_set = checkpoint_file["spinor_set"]
                    else
                        error("we don't have spinor")
                    end
                end

         
                wave_work_dic = Dict{Tuple{Vararg{Int64}}, Int}()

                for wave_work_index in eachindex(wave_work)
                    wave_work_dic[Tuple(wave_work[wave_work_index])] = wave_work_index
                end

                Uoperator = zeros(
                    ComplexF64,
                    length(spinor_set),
                    length(spinor_set)
                )

                for source_index in eachindex(wave_work)
                    shifted_momentum = wave_work[source_index] .+ direction
                    target_index = get(wave_work_dic, Tuple(shifted_momentum), 0)

                    if target_index != 0
                        Uoperator[target_index, source_index] =
                            spinor_set[target_index]' * spinor_set[source_index]
                    end
                end

                matrix_size = size(Prj, 1)

                manybody_operator = Matrix{ComplexF64}(I, matrix_size, matrix_size)
                manybody_operator .-= Prj

                mul!(
                    manybody_operator,
                    Prj,
                    Uoperator,
                    1.0 + 0.0im,
                    1.0 + 0.0im
                )

                Z = det(manybody_operator)


                Z_record[kkk] = Z
            catch err
             error_record[kkk] = sprint(showerror, err)
            end
            
            end


            missing_indices = findall(missing_record)
            error_indices = findall(error_record .!= "")

            println("number of missing checkpoints = ", length(missing_indices))
            println("number of errored checkpoints = ", length(error_indices))

               for error_index in error_indices
                println("errored checkpoint step = ", ckp_list[error_index])
                println(error_record[error_index])
                println()
            end


             # Get the phase unwinded

            phase_record=angle.(Z_record)/(2π)
            phase_unwrapped = copy(phase_record)
            valid_indices = findall(!isnan, phase_record)

            for valid_position in 2:length(valid_indices)
                previous_index = valid_indices[valid_position - 1]
                current_index = valid_indices[valid_position]

                phase_jump = phase_unwrapped[current_index] - phase_unwrapped[previous_index]

                if phase_jump > 0.5
                    for later_valid_position in valid_position:length(valid_indices)
                        later_index = valid_indices[later_valid_position]
                        phase_unwrapped[later_index] -= 1.0
                    end
                end

                if phase_jump < -0.5
                    for later_valid_position in valid_position:length(valid_indices)
                        later_index = valid_indices[later_valid_position]
                        phase_unwrapped[later_index] += 1.0
                    end
                end
            end

            if !isdir(joinpath(base_folder,"Wilsondata"))
                mkdir(joinpath(base_folder,"Wilsondata"))
            end

            save_path=joinpath(base_folder,"Wilsondata","Wilson_dir$(direction[1])dir$(direction[2])_$(args[1])NL$(args[2])theta$(args[3])constq$(args[4])ϵr$(args[5])uD$(args[6])filling$(args[7])cutoff$(args[8])lambda$(args[9])trytime$(args[10])enlarge$(args[11])V0_hBN$(args[12])V1_hBN$(args[13])ψ_hBN$(args[14])V2_scalar$(args[15])ϕ$(args[16])pincof$(args[17])dedis$(args[18])depos" *
                                        "$(args[22])refreshevery" *
                                        "$(args[24])dt$(args[25])Emag$(args[26])Eag$(args[27])gamma$(args[28])temp$(args[29])workcutoff" *
                                        ".jld2")

            last_checkpoint_path=joinpath(base_folder,"$(args[1])NL$(args[2])theta$(args[3])constq$(args[4])ϵr$(args[5])uD$(args[6])filling$(args[7])cutoff$(args[8])lambda$(args[9])trytime$(args[10])enlarge$(args[11])V0_hBN$(args[12])V1_hBN$(args[13])ψ_hBN$(args[14])V2_scalar$(args[15])ϕ$(args[16])pincof$(args[17])dedis$(args[18])depos" *
                    "$(args[22])refreshevery" *
                    "$(args[24])dt$(args[25])Emag$(args[26])Eag$(args[27])gamma$(args[28])temp$(args[29])workcutoff" *
                    "$(ckp_list[valid_indices[end]])stepnum.jld2")

            last_checkpoint=load(last_checkpoint_path)


            diagnostics_record = last_checkpoint["diagnostics_record"]

            step_record = [diagnostics.step_index for diagnostics in diagnostics_record]
            trace_record = [real(diagnostics.trace) for diagnostics in diagnostics_record]




            jldsave(save_path,ckp_list=ckp_list,Z_record=Z_record,valid_indices=valid_indices,phase_unwrapped=phase_unwrapped,step_record=step_record,trace_record=trace_record)

         

            fig1=plot(
                step_record,
                trace_record,
                xlabel = "step",
                ylabel = "Tr(P)",
                label = "Tr(P)",
                linewidth = 2,
                marker = :circle,
                title="nu$filling, E$E_mag,γ$gamma,uD$(uD),VBg$V2_scalar,dt$dt"
            )

            fig2=plot(ckp_list,phase_record,xlabel="step",ylabel="wrapped phase/2pi",title="nu=$filling, E=$E_mag,γ=$gamma,uD=$(uD),VBg=$V2_scalar,dt=$dt",legend=false)


            fig3=plot(dt*E_mag*ckp_list*norm(a1m)/(2π),abs.(Z_record),xlabel="flux/(2π)",ylabel="abs(Z)",title="nu=$filling, E=$E_mag,γ=$gamma,uD=$(uD),VBg=$V2_scalar,dt=$dt",legend=false)


            fig4 = plot(
                dt * E_mag * ckp_list[valid_indices] * norm(a1m) / (2π),
                phase_unwrapped[valid_indices],
                xlabel = "flux/(2π)",
                ylabel = "unwrapped phase/2pi",
                title = "nu=$filling, E=$E_mag,γ=$gamma,uD=$(uD),VBg=$V2_scalar,dt=$dt",
                legend = false
            )


            fig5=plot(dt*E_mag*ckp_list*norm(a1m)/(2π),phase_record,xlabel="flux/(2π)",ylabel="wrapped phase/2pi",title="nu=$filling, E=$E_mag,γ=$gamma,uD=$(uD),VBg=$V2_scalar,dt=$dt",legend=false)


            
            base_plot_name="$(args[1])NL$(args[2])theta$(args[3])constq$(args[4])ϵr$(args[5])uD$(args[6])filling$(args[7])cutoff$(args[8])lambda$(args[9])trytime$(args[10])enlarge$(args[11])V0_hBN$(args[12])V1_hBN$(args[13])ψ_hBN$(args[14])V2_scalar$(args[15])ϕ$(args[16])pincof$(args[17])dedis$(args[18])depos" *
                    "$(args[22])refreshevery" *
                    "$(args[24])dt$(args[25])Emag$(args[26])Eag$(args[27])gamma$(args[28])temp$(args[29])workcutoff" 

        savefig(fig1, joinpath(base_folder,"Wilsondata","Tr_di$(direction[1])_$(direction[2])" * base_plot_name * ".png"))
        savefig(fig2, joinpath(base_folder,"Wilsondata","wp_step_di$(direction[1])_$(direction[2])" * base_plot_name * ".png"))
        savefig(fig3, joinpath(base_folder,"Wilsondata","absZ_flux_di$(direction[1])_$(direction[2])" * base_plot_name * ".png"))
        savefig(fig4, joinpath(base_folder,"Wilsondata","uwp_flux_di$(direction[1])_$(direction[2])" * base_plot_name * ".png"))
        savefig(fig5, joinpath(base_folder,"Wilsondata","wp_flux_di$(direction[1])_$(direction[2])" * base_plot_name * ".png"))
        

        return nothing


end