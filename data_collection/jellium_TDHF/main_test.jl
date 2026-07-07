using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

using LinearAlgebra
using JLD2
using Plots

# args[1]  = gcutoff_seed        # only for bookkeeping / matching old convention
# args[2]  = λ                   # actual value read from seed["arguments"]
# args[3]  = keep_band_1_tdhf    # 1.0 or 0.0
# args[4]  = keep_band_2_tdhf    # 1.0 or 0.0
# args[5]  = enlarge_factor      # actual value read from seed["arguments"]
# args[6]  = V2_scalar           # actual value read from seed["arguments"]
# args[7]  = ϕ                   # actual value read from seed["arguments"]
# args[8]  = pin_coeff           # actual value read from seed["arguments"]
# args[9]  = rs                  # actual value read from seed["arguments"]
# args[10] = dedis               # actual value read from seed["arguments"]
# args[11] = defec_pos           # actual value read from seed["arguments"]
# args[12] = constq              # actual value read from seed["arguments"]
# args[13] = filling             # actual value read from seed["arguments"]
# args[14] = Δ                   # actual value read from seed["arguments"]

# args[15] = output_index        # data_output$(output_index)
# args[16] = start_step          # 0 means seed/checkpoint_step_0; >0 means checkpoint_step_start_step
# args[17] = total_fluxes
# args[18] = steps_per_flux
# args[19] = save_every
# args[20] = refresh_every
# args[21] = gamma
# args[22] = temp
# args[23] = dt
# args[24] = direction_angle_degree
# args[25] = gcutoff_work        # optional; <=0 means use seed cutoff



BLAS.set_num_threads(1)

include("../../src/lambda_jellium_TDHF.jl")

args = parse.(Float64, ARGS)

run_tdhf_from_args!(args)


