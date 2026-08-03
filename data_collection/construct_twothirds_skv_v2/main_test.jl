using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2


include("../../src/construct_twothirds_skv_v2.jl")

args=parse.(Float64,ARGS)

#args =
 #   parse.(Float64, ARGS)


flux1 =args[1]

flux2 =args[2]

NL=Int(args[3])

moiream=args[4]

N1 =Int(args[5])

N2 =Int(args[6])

N1f =Int(args[7])

N2f =Int(args[8])

xi00_re =args[9]

xi00_im =args[10]

grid_cutoff =args[11]

type =Int(args[12])

gaus_damp=args[13]

file_pos =Int(args[14])


(params,
    topo_sectors,
    overall_mag_matrix,
    spinor_set,
    possible_config,
    vac_fac_matrix,
    orbital_Rmatrix_list,
    orbital_basis_norm_list,
    orbital_basis_orthogonal_list,
    M_eta_list,
    Mkkmatrix_list,
    N_coeff_list) =big_func(
        flux1,
        flux2,
        NL,
        moiream,
        N1,
        N2,
        N1f,
        N2f,
        xi00_re + im*xi00_im,
        grid_cutoff,
        type,gaus_damp
    )


scratch_dir =
    ENV["SCRATCH"]


savepath =
    joinpath(
        scratch_dir,
        "constrcut_twothirds_skv_v2/data_output$(file_pos)/" *
        "$(args[1])f1" *
        "$(args[2])f2" *
        "$(args[3])NL" *
        "$(args[4])am" *
        "$(args[5])N1" *
        "$(args[6])N2" *
        "$(args[7])N1f" *
        "$(args[8])N2f" *
        "$(args[9])xir" *
        "$(args[10])xii" *
        "$(args[11])grid" *
        "$(args[12])type$(args[13])gaud.jld2",
    )



jldsave(
    savepath;

    params = params,

    topo_sectors = topo_sectors,

    # --------------------------------------------------------
    # Three many-body states:
    #
    # |Psi_s> = sum_cc overall_mag_matrix[cc,s] |D_cc>
    # --------------------------------------------------------

    overall_mag_matrix =
        overall_mag_matrix,

    # Common single-particle information
    spinor_set =
        spinor_set,

    possible_config =
        possible_config,

    # Sector-dependent scalar factors
    vac_fac_matrix =
        vac_fac_matrix,

    # --------------------------------------------------------
    # Common fixed-vacancy Slater determinant data.
    # These are NOT duplicated for the three sectors.
    # --------------------------------------------------------

    orbital_Rmatrix_list =
        orbital_Rmatrix_list,

    orbital_basis_norm_list =
        orbital_basis_norm_list,

    orbital_basis_orthogonal_list =
        orbital_basis_orthogonal_list,

    M_eta_list =
        M_eta_list,

    Mkkmatrix_list =
        Mkkmatrix_list,

    N_coeff_list =
        N_coeff_list,
)
