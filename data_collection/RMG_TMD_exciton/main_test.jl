using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
# I think in this calculation I should still be keeping two spin and two valley of RMG. Since we want to see in the case of C=5 phase,
# where will the pairing happen.  We should also keep two spins for the TMD. Since it's spin valley locked.

# This code is modified from double_RMG_allow_coherence

#update:I allow arbitrary flavor to be kept
using JLD2,LinearAlgebra
include("../../src/RMG_TMD_exciton.jl")

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


selection = get_selection(sel) 



basis, k_set, k_index, Area, single_matrix, single_particle_wf =get_single_particle(radius,
        num_points,
        uD,
        deltaE,
        NL,
        m_TMD,
        z_TMD;
        selection = selection,
    )

println("finish1")

initial_density_matrix, BG_density_matrix =
    get_initial_proj(
        k_set,
        basis,
        single_matrix,
        NL,
        m_TMD,
    )


HF_eigenvalues,HF_eigenvectors,energy,
DIIS_input_density_matrix,fermi_level,
Hartree_matrix,Fock_matrix,
eout,renormalized_density = iteration(
    initial_density_matrix,
    BG_density_matrix,
    ϵr,
    k_set,
    single_matrix,
    basis,
    Area,
    NL,
    z_TMD,
    target_density,
    temp,
)

scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "RMG_TMD_exciton/data_output$(Int(args[13]))/$(args[1])radius$(args[2])num_points$(args[3])uD$(args[4])deltaE$(args[5])er$(args[6])ztmd$(args[7])mTMD$(args[8])NL$(args[9])tgden$(args[10])temp$(args[11])sel$(args[12])trytime.jld2")

  
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

    # HF output
    BG_density_matrix = BG_density_matrix,
    density_matrix = DIIS_input_density_matrix[1],
    HF_eigenvalues = HF_eigenvalues,
    HF_eigenvectors = HF_eigenvectors,
    Hartree_matrix = Hartree_matrix,
    Fock_matrix = Fock_matrix,
    fermi_level = fermi_level,
    energy = energy,
    eout = eout,
    renormalized_density = renormalized_density,

    # reproducibility
    arguments = args,
)