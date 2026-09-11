

 #V2 implements flavor-dependent displacement field, and valence/conduction band projection
using Pkg
Pkg.activate(joinpath(@__DIR__,"../.."))
using JLD2

include("../../src/MWC_v3.jl")

args=parse.(Float64,ARGS)
n=args[1]                              # requested total signed carrier density, nm^-2
Area=args[2]                           # crystal unit-cell area, nm^2
ϵr=args[3]
uD=args[4]
flavor_code=Int(args[5])
carrier_population_code=Int(args[6])
lattice_code=Int(args[7])
θ=args[8]/180*π
gate_distance=args[9]                  # nm
NL=Int(args[10])
Nq=Int(args[11])
Nb=Int(args[12])
gcutoff=args[13]                       # nm^-1
trytime=Int(args[14])
displacement_sign_code=Int(args[15])
band_projection_code=Int(args[16])     # 1: conduction, 2: valence
filepos=Int(args[17])

stacking=1
maximum_iterations=500000000
density_tolerance=10^(-14)
energy_tolerance=1e-10





#Set the parameters

valley_by_flavor=get_valley_by_flavor(flavor_code)
Nf=length(valley_by_flavor)


displacement_sign_by_flavor =
    get_displacement_sign_by_flavor(
        displacement_sign_code,
        Nf
    )

if band_projection_code == 1
    project_to_valence = false
    reference_occupation = 0
elseif band_projection_code == 2
    project_to_valence = true
    reference_occupation = 1
else
    error("band_projection_code must be 1 for conduction or 2 for valence")
end




Area>0 || error("Area must be positive")

a1,a2,b1,b2=get_lattice_vectors(lattice_code,Area,θ)

momentum_mesh=MomentumMesh(b1,b2,Nq,Nq)
Nk=momentum_mesh.momentum_count

carrier_fraction_by_flavor =
    get_carrier_fraction_by_flavor(carrier_population_code, valley_by_flavor)

requested_carrier_density_by_flavor =
    n .* carrier_fraction_by_flavor

signed_carrier_count_by_flavor =
    round.(Int, requested_carrier_density_by_flavor .* (Area * Nk))

reference_occupied_state_count_by_flavor =
    fill(reference_occupation * Nb * Nk, Nf)

occupied_state_count_by_flavor =
    reference_occupied_state_count_by_flavor .+ signed_carrier_count_by_flavor


occupied_state_count = sum(occupied_state_count_by_flavor)
reference_occupied_state_count = sum(reference_occupied_state_count_by_flavor)

signed_carrier_count =
    sum(signed_carrier_count_by_flavor)

carrier_count =
    abs(signed_carrier_count)

filling =
    n * Area

filling_real_by_flavor =
    signed_carrier_count_by_flavor ./ Nk

filling_real =
    sum(filling_real_by_flavor)

nreal_by_flavor =
    signed_carrier_count_by_flavor ./ (Nk * Area)

nreal =
    sum(nreal_by_flavor)




level_count_per_flavor = Nb * Nk

all(occupied_state_count -> 0 <= occupied_state_count <= level_count_per_flavor,
    occupied_state_count_by_flavor) ||
    error("At least one flavor occupation is outside the projected Hilbert space")

carrier_count > 0 ||
    error("carrier_count must be positive")

folding_coordinates=get_folding_coordinates(momentum_mesh,gcutoff)

projected_basis = get_projected_basis(
    momentum_mesh,
    folding_coordinates,
    Nb,
    NL,
    uD,
    valley_by_flavor,
    displacement_sign_by_flavor,
    stacking,
    project_to_valence
)

single_Ham=get_single_particle_hamiltonian(projected_basis)
form_factors=precompute_form_factors(projected_basis)



result = run_hartree_fock_with_oda(
    projected_basis, momentum_mesh, single_Ham, form_factors,
    ϵr, gate_distance, occupied_state_count_by_flavor;
    reference_occupation=reference_occupation,
    maximum_iterations=maximum_iterations,
    density_tolerance=density_tolerance,
    energy_tolerance=energy_tolerance,
    verbose=true
)

println(result.solution)






densitymatrix=result.density_matrix
energy=result.solution.energy_per_carrier
eout=result.solution

raw_HF_eigenvalue = result.filling_work.eigenvalues
raw_HF_eigenvector = result.filling_work.eigenvectors

HFeigenvalue = similar(raw_HF_eigenvalue)
HF_eigenvector = similar(raw_HF_eigenvector)
HF_flavor_index = zeros(Int, size(raw_HF_eigenvalue))

for k_index in axes(raw_HF_eigenvalue, 2)
    energy_order = sortperm(@view raw_HF_eigenvalue[:, k_index])

    HFeigenvalue[:, k_index] .= raw_HF_eigenvalue[energy_order, k_index]
    HF_eigenvector[:, :, k_index] .= raw_HF_eigenvector[:, energy_order, k_index]
    HF_flavor_index[:, k_index] .= cld.(energy_order, Nb)
end









HartreeMatrix=result.hamiltonian_work.hartree_hamiltonian
FockMatrix=result.hamiltonian_work.fock_hamiltonian

single_eigenvalue=projected_basis.energy
single_eigenvector=projected_basis.spinor
selected_folding_index=projected_basis.selected_folding_index
wave=momentum_mesh.momenta
occupied_count_by_k=result.filling_work.occupied_count_by_k
occupied_count_by_k_and_flavor =
    result.filling_work.occupied_count_by_k_and_flavor

scratch_dir=ENV["SCRATCH"]
savedir=joinpath(scratch_dir,"MWC_v3/data_output$(Int(filepos))")
mkpath(savedir)

savepath=joinpath(savedir,"$(args[1])n$(args[2])Area$(args[3])er$(args[4])uD$(args[5])flavor$(args[6])population$(args[7])lattice$(args[8])theta$(args[9])d$(args[10])NL$(args[11])Nq$(args[12])Nb$(args[13])cutoff$(args[14])try$(args[15])uDsign$(args[16])band.jld2")

jldsave(savepath;
    single_Ham=single_Ham, single_eigenvalue=single_eigenvalue,
    single_eigenvector=single_eigenvector, arguments=args,
    densitymatrix=densitymatrix, energy=energy, eout=eout,
    HFeigenvalue=HFeigenvalue, HF_eigenvector=HF_eigenvector,HF_flavor_index=HF_flavor_index,
    HartreeMatrix=HartreeMatrix, FockMatrix=FockMatrix,
    folding_coordinates=folding_coordinates,
    selected_folding_index=selected_folding_index, wave=wave,
    Area=Area, a1=a1, a2=a2, b1=b1, b2=b2,
    n=n, nreal=nreal, nreal_by_flavor=nreal_by_flavor,
    filling_real=filling_real, filling_real_by_flavor=filling_real_by_flavor,
    carrier_count=carrier_count, signed_carrier_count=signed_carrier_count,
    signed_carrier_count_by_flavor=signed_carrier_count_by_flavor,
    requested_carrier_density_by_flavor=requested_carrier_density_by_flavor,
    carrier_fraction_by_flavor=carrier_fraction_by_flavor,
    occupied_state_count=occupied_state_count,
    occupied_state_count_by_flavor=occupied_state_count_by_flavor,
    occupied_count_by_k=occupied_count_by_k,
    occupied_count_by_k_and_flavor=occupied_count_by_k_and_flavor,
    reference_occupation=reference_occupation,
    reference_occupied_state_count=reference_occupied_state_count,
    reference_occupied_state_count_by_flavor=reference_occupied_state_count_by_flavor,
    valley_by_flavor=valley_by_flavor,
    carrier_population_code=carrier_population_code,
    project_to_valence=project_to_valence,
    band_projection_code=band_projection_code,
    displacement_sign_code=displacement_sign_code,
    displacement_sign_by_flavor=displacement_sign_by_flavor
)

println("Saved to $savepath")