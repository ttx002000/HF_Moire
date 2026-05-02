using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

using JLD2,LinearAlgebra
include("../../src/operators_RMG_TMD_exciton_bandprojected.jl")

args=parse.(Float64,ARGS)
radius=args[1]
num_points=Int(args[2])
uD=args[3]
deltaE=args[4]
ϵr=args[5]
z_TMD=args[6]
m_TMD=args[7]
NL=Int(args[8])
target_density=args[9]
temp=args[10]
sel=Int(args[11])
trytime=args[12]
file_pos=args[13]


selection = get_projected_selection(
    sel,
    NL,
    deltaE,
    m_TMD,
    z_TMD;
) 




model = build_projected_model(
    radius,
    num_points,
    uD,
    NL,
    selection,
    ϵr;
    stacking = 1,
    
)

initial_density_matrix, BG_density_matrix =
    get_initial_proj(
        model,
        NL;
        init_amp = 0.01,
    )

HF_eigenvalues,
HF_eigenvectors,
energy,
DIIS_input_density_matrix,
fermi_level,
Hartree_matrix,
Fock_matrix,
eout,
renormalized_density = iteration(
    initial_density_matrix,
    BG_density_matrix,
    model,
    target_density,
    temp;
    DIIS_size = 5,
    mixing = 0.5,
)

scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "RMG_TMD_exciton_bandprojected/data_output$(Int(args[13]))/$(args[1])radius$(args[2])num_points$(args[3])uD$(args[4])deltaE$(args[5])er$(args[6])ztmd$(args[7])mTMD$(args[8])NL$(args[9])tgden$(args[10])temp$(args[11])sel$(args[12])trytime.jld2")

  
jldsave(savepath;
    # lightweight model data
    active = model.active,
    k_set = model.k_set,
    k_index = model.k_index,
    Area = model.Area,
    single_matrix = model.single_matrix,
    U_by_flavor = model.U_by_flavor,

    active_flavor = model.active_flavor,
    local_band_of_active = model.local_band_of_active,
    active_in_flavor = model.active_in_flavor,
    flavor_ranges = model.flavor_ranges,
    z_flavor = model.z_flavor,
    Vzero = model.Vzero,
    flavor_key_to_id = model.flavor_key_to_id,

    # initial / background / converged density
    BG_density_matrix = BG_density_matrix,
    density_matrix = DIIS_input_density_matrix[1],

    # HF output
    HF_eigenvalues = HF_eigenvalues,
    HF_eigenvectors = HF_eigenvectors,
    Hartree_matrix = Hartree_matrix,
    Fock_matrix = Fock_matrix,
    fermi_level = fermi_level,
    energy = energy,
    eout = eout,
    renormalized_density = renormalized_density,

    # parameters
    arguments = args,
)