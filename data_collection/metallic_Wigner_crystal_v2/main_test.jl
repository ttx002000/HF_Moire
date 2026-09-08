

 #V2 implements flavor-dependent displacement field, and valence/conduction band projection
using Pkg
Pkg.activate(joinpath(@__DIR__,"../.."))
using JLD2

include("../../src/MWC_v2.jl")

args=parse.(Float64,ARGS)

n=args[1]                    # total density, nm^-2
δn=args[2]                   # requested density mismatch, nm^-2
m=Int(args[3])               # integer crystal filling
ϵr=args[4]
uD=args[5]
flavor_code=Int(args[6])
lattice_code=Int(args[7])
θ=args[8]/180*π
gate_distance=args[9]        # nm
NL=Int(args[10])
Nq=Int(args[11])
Nb=Int(args[12])
gcutoff=args[13]   
trytime=Int(args[14])  
displacement_sign_code = Int(args[15])
band_projection_code=Int(args[16])        # 1: conduction, 2: valence
filepos=Int(args[17])

stacking=1
seed_strength=1e-2
random_seed=rand(1:typemax(Int))
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




ncrystal = n - δn

ncrystal != 0 ||
    error("n-δn cannot be zero")

Area = m / ncrystal

Area > 0 ||
    error("m and n-δn must have the same sign")

a1,a2,b1,b2=get_lattice_vectors(lattice_code,Area,θ)

momentum_mesh=MomentumMesh(b1,b2,Nq,Nq)
Nk=momentum_mesh.momentum_count

filling = n * Area

signed_carrier_count =
    round(Int, filling * Nk)

reference_occupied_state_count =
    reference_occupation * Nb * Nf * Nk

occupied_state_count =
    reference_occupied_state_count +
    signed_carrier_count

carrier_count =
    abs(signed_carrier_count)

filling_real =
    signed_carrier_count / Nk

nreal =
    signed_carrier_count / (Nk * Area)

δnreal =
    nreal - ncrystal


0 < occupied_state_count < Nb * Nf * Nk ||
    error("occupied_state_count is outside the projected Hilbert space")

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
    projected_basis,
    momentum_mesh,
    single_Ham,
    form_factors,
    ϵr,
    gate_distance,
    occupied_state_count;
    reference_occupation=reference_occupation,
    seed_strength=seed_strength,
    random_seed=random_seed,
    maximum_iterations=maximum_iterations,
    density_tolerance=density_tolerance,
    energy_tolerance=energy_tolerance,
    verbose=true
)

println(result.solution)






densitymatrix=result.density_matrix
energy=result.solution.energy_per_carrier
eout=result.solution

HFeigenvalue=result.filling_work.eigenvalues
HF_eigenvector=result.filling_work.eigenvectors
HartreeMatrix=result.hamiltonian_work.hartree_hamiltonian
FockMatrix=result.hamiltonian_work.fock_hamiltonian

single_eigenvalue=projected_basis.energy
single_eigenvector=projected_basis.spinor
selected_folding_index=projected_basis.selected_folding_index
wave=momentum_mesh.momenta
occupied_count_by_k=result.filling_work.occupied_count_by_k
scratch_dir=ENV["SCRATCH"]
savedir=joinpath(scratch_dir,"MWC_v2/data_output$(Int(filepos))")
mkpath(savedir)

savepath=joinpath(savedir,"$(args[1])n$(args[2])dn$(args[3])m$(args[4])er$(args[5])uD$(args[6])flavor$(args[7])lattice$(args[8])theta$(args[9])d$(args[10])NL$(args[11])Nq$(args[12])Nb$(args[13])cutoff$(args[14])try$(args[15])uDsign$(args[16])band.jld2")
jldsave(savepath;single_Ham=single_Ham,single_eigenvalue=single_eigenvalue,
        single_eigenvector=single_eigenvector,arguments=args,densitymatrix=densitymatrix,
        energy=energy,eout=eout,HFeigenvalue=HFeigenvalue,HF_eigenvector=HF_eigenvector,
        HartreeMatrix=HartreeMatrix,FockMatrix=FockMatrix,folding_coordinates=folding_coordinates,
        selected_folding_index=selected_folding_index,wave=wave,Area=Area,
        a1=a1,a2=a2,b1=b1,b2=b2,occupied_state_count=occupied_state_count,
        signed_carrier_count=signed_carrier_count,carrier_count=carrier_count,
        nreal=nreal,δnreal=δnreal,
        filling_real=filling_real,valley_by_flavor=valley_by_flavor,occupied_count_by_k=occupied_count_by_k,
        reference_occupation=reference_occupation,
        reference_occupied_state_count=reference_occupied_state_count,
        project_to_valence=project_to_valence,
        band_projection_code=band_projection_code,
        displacement_sign_code=displacement_sign_code,
        displacement_sign_by_flavor=displacement_sign_by_flavor,)

println("Saved to $savepath")