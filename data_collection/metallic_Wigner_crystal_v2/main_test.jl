

 #V2 implements flavor-dependent displacement field, and valence/conduction band projection
using Pkg
Pkg.activate(joinpath(@__DIR__,"../.."))
using JLD2

include("../../src/MWC_v2.jl")

args=parse.(Float64,ARGS)

n=args[1]                         # signed carrier density, nm^-2
Area=args[2]                      # crystal unit-cell area, nm^2
ϵr=args[3]
uD=args[4]
flavor_code=Int(args[5])
lattice_code=Int(args[6])
θ=args[7]/180*π
gate_distance=args[8]             # nm
NL=Int(args[9])
Nq=Int(args[10])
Nb=Int(args[11])
gcutoff=args[12]                  # nm^-1
trytime=Int(args[13])
displacement_sign_code=Int(args[14])
band_projection_code=Int(args[15]) # 1: conduction, 2: valence
filepos=Int(args[16])

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

savepath=joinpath(savedir,"$(args[1])n$(args[2])Area$(args[3])er$(args[4])uD$(args[5])flavor$(args[6])lattice$(args[7])theta$(args[8])d$(args[9])NL$(args[10])Nq$(args[11])Nb$(args[12])cutoff$(args[13])try$(args[14])uDsign$(args[15])band.jld2")
jldsave(savepath;single_Ham=single_Ham,single_eigenvalue=single_eigenvalue,
        single_eigenvector=single_eigenvector,arguments=args,densitymatrix=densitymatrix,
        energy=energy,eout=eout,HFeigenvalue=HFeigenvalue,HF_eigenvector=HF_eigenvector,
        HartreeMatrix=HartreeMatrix,FockMatrix=FockMatrix,folding_coordinates=folding_coordinates,
        selected_folding_index=selected_folding_index,wave=wave,Area=Area,
        a1=a1,a2=a2,b1=b1,b2=b2,occupied_state_count=occupied_state_count,
        signed_carrier_count=signed_carrier_count,carrier_count=carrier_count,
        n=n,nreal=nreal,
        filling_real=filling_real,valley_by_flavor=valley_by_flavor,occupied_count_by_k=occupied_count_by_k,
        reference_occupation=reference_occupation,
        reference_occupied_state_count=reference_occupied_state_count,
        project_to_valence=project_to_valence,
        band_projection_code=band_projection_code,
        displacement_sign_code=displacement_sign_code,
        displacement_sign_by_flavor=displacement_sign_by_flavor,)

println("Saved to $savepath")