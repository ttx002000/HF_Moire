using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2


include("../../src/construct_twothirds_Lauhglin_ED_step3_v2.jl")

args=parse.(Float64,ARGS)


#=
    flux1=args[1]
     flux2=args[2]
     q1=args[3]
     q2=args[4]
     N1=Int(args[5])
   N2=Int(args[6])
   Npa_LL=Int(args[7])
   Nvec=Int(args[8])
   file_pos=Int(args[9])
   grid_cutoff=args[10]
   NL=Int(args[11])
   moiream=args[12]
   gausdamp=args[13]
=#

st = do_Gq(args)



     



scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "constrcut_twothirds_Lauhglin_ED_v2/data_output$(Int(args[9]))/second_minor/SM_$(args[1])f1$(args[2])f2$(args[3])q1$(args[4])q2$(args[5])N1$(args[6])N2$(args[7])Npa$(args[8])Nvec$(args[10])gcut$(args[11])NL$(args[12])am$(args[13])gda.jld2")



mkpath(dirname(savepath))

jldsave(
    savepath;

    params = st.params,

    blocks = st.blocks,
    block_config_ranges = st.block_config_ranges,
    possible_config = st.possible_config,

    orbital_basis_orthogonal_list =
        st.orbital_basis_orthogonal_list,

    overall_mag_matrix =
        st.overall_mag_matrix,

    reconstruct_matrix =
        st.reconstruct_matrix,

    first_minor_matrix =
        st.first_minor_matrix,

    deno_raw =
        st.deno_raw,

    Gq_raw_array =
        st.Gq_raw_array,

    rho_raw =
        st.rho_raw,

    deno_orth =
        st.deno_orth,

    Gq_orth_array =
        st.Gq_orth_array,

    rho_orth =
        st.rho_orth,

    qkeys =
        st.qkeys,

    qindex =
        st.qindex,

    minus_q_index =
        st.minus_q_index,

    spinor_set =
        st.spinor_set,

    second_minor_matrix=st.second_minor_matrix,
)


