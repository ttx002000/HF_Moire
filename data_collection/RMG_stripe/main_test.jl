using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

using JLD2
using LinearAlgebra
using Random

include("../../src/RMG_stripe_phase.jl")
args = parse.(Float64, ARGS)

radius = args[1]
num_points = Int(args[2])
uD = args[3]
ϵr = args[4]
NL = Int(args[5])
target_density = args[6]
temp = args[7]

grid_angle = args[8]
Qvec_index = Int[args[9], args[10]]
band_index = Int(args[11])

trytime = Int(args[12])
file_pos = Int(args[13])


println("Building stripe single-particle bookkeeping...")
flush(stdout)

book = get_single_particle(
    radius,
    num_points,
    uD,
    NL,
    grid_angle,
    Qvec_index;
    band_index = band_index,
)


println("Finished single-particle bookkeeping.")
println("Nk full = ", length(book.kvecs_cart))
println("Nk reduced = ", length(book.k_equiv_id))
println("Nn max = ", length(book.nvals))
println("Qint = ", book.Qint)
println("Qvec_cart = ", book.Qvec_cart)
println("Area = ", book.Area)
flush(stdout)



initial_density_matrix = random_hermitian_density_matrix(
    book
)


println("Finished random density matrix initialization.")
flush(stdout)


HF_eigenvalues,
HF_eigenvectors,
energy,
final_density_matrix,
fermi_level,
Hartree_matrix,
Fock_matrix,
eout,
renormalized_density = iteration_stripe(
    initial_density_matrix,
    book,
    ϵr,
    target_density,
    temp,
    band_index,
    NL;
)


scratch_dir = ENV["SCRATCH"]

output_dir = joinpath(
    scratch_dir,
    "RMG_stripe_phase/data_output$(Int(args[13]))",
)

mkpath(output_dir)

savepath = joinpath(
    output_dir,
    "$(args[1])radius" *
    "$(args[2])num_points" *
    "$(args[3])uD" *
    "$(args[4])er" *
    "$(args[5])NL" *
    "$(args[6])tgden" *
    "$(args[7])temp" *
    "$(args[8])angle" *
    "$(args[9])Qx" *
    "$(args[10])Qy" *
    "$(args[11])band" *
    "$(args[12])trytime.jld2"
)

println("Saving to: ", savepath)
flush(stdout)

jldsave(
    savepath;

    # final HF output
    final_density_matrix = final_density_matrix,
    initial_density_matrix = initial_density_matrix,
    HF_eigenvalues = HF_eigenvalues,
    HF_eigenvectors = HF_eigenvectors,
    Hartree_matrix = Hartree_matrix,
    Fock_matrix = Fock_matrix,
    energy = energy,
    fermi_level = fermi_level,
    eout = eout,
    renormalized_density = renormalized_density,

    # just save raw args
    args = args,

    # unpacked StripeSPBookkeeping fields, excluding huge form-factor objects
    band_energies = book.band_energies,
    band_vecs = book.band_vecs,
    single_particle_matrix = book.single_particle_matrix,

    kvecs_cart = book.kvecs_cart,
    k_int = book.k_int,
    kint_to_gid = book.kint_to_gid,

    Qint = book.Qint,
    Qvec_cart = book.Qvec_cart,

    k_equiv_id = book.k_equiv_id,
    k_equiv_n = book.k_equiv_n,
    k_repre_id = book.k_repre_id,

    nvals = book.nvals,
    gid_padded = book.gid_padded,
    valid_ns = book.valid_ns,

    Area = book.Area,
)

println("Done.")
flush(stdout)