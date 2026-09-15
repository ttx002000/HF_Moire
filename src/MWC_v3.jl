using LinearAlgebra
using Random
using JLD2
# v3 fixes the carrier population of each flavor separately and forbids flavor coherence.

# Helper

function get_carrier_fraction_by_flavor(carrier_population_code::Int,
                                        valley_by_flavor::AbstractVector{<:Integer})
    flavor_count = length(valley_by_flavor)

    carrier_population_code == 1 &&
        return fill(1.0 / flavor_count, flavor_count)

    error("carrier_population_code must currently be 1 for equal flavor populations")
end


function get_valley_by_flavor(flavor_code::Int)
    flavor_code==1 && return [1]
    flavor_code==2 && return [-1]
    flavor_code==3 && return [1,-1]
    flavor_code==4 && return [1,1]
    flavor_code==5 && return [-1,-1]
    flavor_code==6 && return [1,1,-1,-1]
    error("flavor_code must be between 1 and 6")
end


function get_displacement_sign_by_flavor(
    displacement_sign_code::Int,
    flavor_count::Int
)
    if flavor_count == 1
        displacement_sign_code == 1 && return [1]
        displacement_sign_code == 2 && return [-1]

    elseif flavor_count == 2
        displacement_sign_code == 1 && return [1, 1]
        displacement_sign_code == 2 && return [1, -1]
        displacement_sign_code == 3 && return [-1, 1]
        displacement_sign_code == 4 && return [-1, -1]

    elseif flavor_count == 4
        displacement_sign_code == 1  && return [1, 1, 1, 1]
        displacement_sign_code == 2  && return [1, 1, 1, -1]
        displacement_sign_code == 3  && return [1, 1, -1, 1]
        displacement_sign_code == 4  && return [1, 1, -1, -1]
        displacement_sign_code == 5  && return [1, -1, 1, 1]
        displacement_sign_code == 6  && return [1, -1, 1, -1]
        displacement_sign_code == 7  && return [1, -1, -1, 1]
        displacement_sign_code == 8  && return [1, -1, -1, -1]
        displacement_sign_code == 9  && return [-1, 1, 1, 1]
        displacement_sign_code == 10 && return [-1, 1, 1, -1]
        displacement_sign_code == 11 && return [-1, 1, -1, 1]
        displacement_sign_code == 12 && return [-1, 1, -1, -1]
        displacement_sign_code == 13 && return [-1, -1, 1, 1]
        displacement_sign_code == 14 && return [-1, -1, 1, -1]
        displacement_sign_code == 15 && return [-1, -1, -1, 1]
        displacement_sign_code == 16 && return [-1, -1, -1, -1]
    end

    error(
        "Invalid displacement_sign_code=$displacement_sign_code " *
        "for flavor_count=$flavor_count"
    )
end





function get_lattice_vectors(lattice_code::Int,Area::Float64,θ::Float64)
    if lattice_code==1
        L=sqrt(Area)
        a1=L*[1.0,0.0]
        a2=L*[0.0,1.0]
    elseif lattice_code==2
        L=sqrt(2Area/√3)
        a1=L*[1.0,0.0]
        a2=L*[1/2,√3/2]
    else
        error("lattice_code must be 1 for square or 2 for triangular")
    end

    rotation=[cos(θ) -sin(θ);sin(θ) cos(θ)]
    a1=rotation*a1
    a2=rotation*a2
    reciprocal_vectors=2π*inv(hcat(a1,a2))'

    return a1,a2,reciprocal_vectors[:,1],reciprocal_vectors[:,2]
end









# Piece 1: Diagonalization and zero-temperature filling
BLAS.set_num_threads(1)

@inline function projected_band_flavor_index(projected_band_index::Int, flavor_index::Int, projected_band_count::Int)
    return projected_band_index + (flavor_index - 1) * projected_band_count
end


# The goal of this piece of code is to do the calculation for specifically Nq=1 case(v4 is also for spinless).
function get_Ham(k::Vector{Float64},uD::Float64,valley::Int64,stacking::Int,NL::Int)
 
   Ham=zeros(ComplexF64,2*NL,2*NL)
  Kac=4π/(3*0.246)*[1,0]*valley
  t0=3100
  t1=380
  t2=-21
  t3=290
  t4=141
  for layer in 1:NL-1
     Ham[2*layer-1:2*layer,2*layer+1:2*layer+2]=[t4*get_f((k+Kac)*stacking) t3*conj(get_f((k+Kac)*stacking));t1 t4*get_f((k+Kac)*stacking)]
  end

  if NL>2
   for layer in 1:NL-2
      Ham[2*layer-1:2*layer,2*layer+3:2*layer+4]=[0.0 t2/2;0.0 0.0]
   end
 end

  Ham=Ham+Ham'

  for layer in 1:NL
      Ham[2*layer-1:2*layer,2*layer-1:2*layer]=[uD*(layer-(NL+1)/2) -t0*get_f((k+Kac)*stacking);-t0*conj(get_f((k+Kac)*stacking)) uD*(layer-(NL+1)/2)]
  end

  
  
 return Ham
end



function get_f(k::Vector{Float64})
  delta1=1/√3*0.246*[0,1]
  delta2=1/√3*0.246*[√3/2,-1/2]
  delta3=1/√3*0.246*[-√3/2,-1/2]

  return exp(im*dot(k,delta1))+exp(im*dot(k,delta2))+exp(im*dot(k,delta3))
end



struct FillingDiagnostics
    occupied_state_count_by_flavor::Vector{Int}
    highest_occupied_energy_by_flavor::Vector{Float64}
    lowest_unoccupied_energy_by_flavor::Vector{Float64}
    chemical_potential_by_flavor::Vector{Float64}
    boundary_gap_by_flavor::Vector{Float64}
end


mutable struct HartreeFockFillingWork
    projected_band_count::Int
    flavor_count::Int
    eigenvalues::Matrix{Float64}
    eigenvectors::Array{ComplexF64,3}
    energy_sorted_level_indices::Vector{Int}
    occupied_count_by_k::Vector{Int}
    occupied_count_by_k_and_flavor::Matrix{Int}
    trial_density_matrix::Array{ComplexF64,3}
end


function HartreeFockFillingWork(projected_band_count::Int, flavor_count::Int, momentum_count::Int)
    hartree_fock_dimension = projected_band_count * flavor_count
    level_count_per_flavor = projected_band_count * momentum_count

    eigenvalues = zeros(Float64, hartree_fock_dimension, momentum_count)
    eigenvectors = zeros(ComplexF64, hartree_fock_dimension, hartree_fock_dimension, momentum_count)
    energy_sorted_level_indices = Vector{Int}(undef, level_count_per_flavor)
    occupied_count_by_k = zeros(Int, momentum_count)
    occupied_count_by_k_and_flavor = zeros(Int, momentum_count, flavor_count)
    trial_density_matrix = zeros(ComplexF64, hartree_fock_dimension, hartree_fock_dimension, momentum_count)

    return HartreeFockFillingWork(projected_band_count, flavor_count, eigenvalues, eigenvectors,
                                  energy_sorted_level_indices, occupied_count_by_k,
                                  occupied_count_by_k_and_flavor, trial_density_matrix)
end



function diagonalize_hartree_fock_hamiltonian!(filling_work::HartreeFockFillingWork,
                                               hartree_fock_hamiltonian::Array{ComplexF64,3})
    projected_band_count = filling_work.projected_band_count
    flavor_count = filling_work.flavor_count
    momentum_count = size(filling_work.eigenvalues, 2)
    hartree_fock_dimension = projected_band_count * flavor_count

    expected_size = (hartree_fock_dimension, hartree_fock_dimension, momentum_count)
    size(hartree_fock_hamiltonian) == expected_size ||
        throw(DimensionMismatch("The Hartree-Fock Hamiltonian has inconsistent dimensions."))

    fill!(filling_work.eigenvectors, 0)

    diagonalization_count = momentum_count * flavor_count

    Threads.@threads :static for diagonalization_index in 1:diagonalization_count
        k_index = mod1(diagonalization_index, momentum_count)
        flavor_index = fld(diagonalization_index - 1, momentum_count) + 1

        first_state_index = projected_band_flavor_index(1, flavor_index, projected_band_count)
        last_state_index = projected_band_flavor_index(projected_band_count, flavor_index,
                                                       projected_band_count) # This relies on the assumption that index within a flavor are close to each other. Which is delibrate, but I will keep it.
        flavor_state_indices = first_state_index:last_state_index

        flavor_hamiltonian = Hermitian(@view hartree_fock_hamiltonian[
            flavor_state_indices, flavor_state_indices, k_index
        ])
        flavor_eigen_decomposition = eigen(flavor_hamiltonian)

        filling_work.eigenvalues[flavor_state_indices, k_index] .=
            flavor_eigen_decomposition.values

        filling_work.eigenvectors[flavor_state_indices, flavor_state_indices, k_index] .=
            flavor_eigen_decomposition.vectors
    end

    return nothing
end




function construct_zero_temperature_density_matrix!(
    filling_work::HartreeFockFillingWork,
    occupied_state_count_by_flavor::AbstractVector{<:Integer}
)
    projected_band_count = filling_work.projected_band_count
    flavor_count = filling_work.flavor_count
    momentum_count = size(filling_work.eigenvalues, 2)
    level_count_per_flavor = projected_band_count * momentum_count

    length(occupied_state_count_by_flavor) == flavor_count ||
        throw(DimensionMismatch("occupied_state_count_by_flavor must contain one entry per flavor."))

    occupied_state_counts = Int.(occupied_state_count_by_flavor)

    all(occupied_state_count -> 0 <= occupied_state_count <= level_count_per_flavor,
        occupied_state_counts) ||
        throw(ArgumentError("Each flavor occupation must lie between 0 and $level_count_per_flavor."))

    fill!(filling_work.occupied_count_by_k, 0)
    fill!(filling_work.occupied_count_by_k_and_flavor, 0)
    fill!(filling_work.trial_density_matrix, 0)

    highest_occupied_energy_by_flavor = fill(NaN, flavor_count)
    lowest_unoccupied_energy_by_flavor = fill(NaN, flavor_count)
    chemical_potential_by_flavor = fill(NaN, flavor_count)
    boundary_gap_by_flavor = fill(NaN, flavor_count)

    for flavor_index in 1:flavor_count
        occupied_state_count = occupied_state_counts[flavor_index]

        first_state_index = projected_band_flavor_index(1, flavor_index, projected_band_count)
        last_state_index = projected_band_flavor_index(projected_band_count, flavor_index,
                                                       projected_band_count)
        flavor_state_indices = first_state_index:last_state_index

        flattened_eigenvalues =
            vec(@view filling_work.eigenvalues[flavor_state_indices, :])

        sortperm!(filling_work.energy_sorted_level_indices, flattened_eigenvalues;
                  alg=Base.Sort.MergeSort)

        for filling_index in 1:occupied_state_count
            linear_level_index = filling_work.energy_sorted_level_indices[filling_index]
            k_index = fld(linear_level_index - 1, projected_band_count) + 1

            filling_work.occupied_count_by_k_and_flavor[k_index, flavor_index] += 1
            filling_work.occupied_count_by_k[k_index] += 1
        end

        for k_index in 1:momentum_count
            occupied_count = filling_work.occupied_count_by_k_and_flavor[k_index, flavor_index]
            occupied_count == 0 && continue

            occupied_column_indices = first_state_index:(first_state_index + occupied_count - 1)
            occupied_eigenvectors = @view filling_work.eigenvectors[
                flavor_state_indices, occupied_column_indices, k_index
            ]
            density_matrix_block = @view filling_work.trial_density_matrix[
                flavor_state_indices, flavor_state_indices, k_index
            ]

            mul!(density_matrix_block, occupied_eigenvectors, adjoint(occupied_eigenvectors))
        end

        if occupied_state_count > 0
            highest_occupied_index =
                filling_work.energy_sorted_level_indices[occupied_state_count]

            highest_occupied_energy_by_flavor[flavor_index] =
                flattened_eigenvalues[highest_occupied_index]
        end

        if occupied_state_count < level_count_per_flavor
            lowest_unoccupied_index =
                filling_work.energy_sorted_level_indices[occupied_state_count + 1]

            lowest_unoccupied_energy_by_flavor[flavor_index] =
                flattened_eigenvalues[lowest_unoccupied_index]
        end

        if 0 < occupied_state_count < level_count_per_flavor
            highest_occupied_energy = highest_occupied_energy_by_flavor[flavor_index]
            lowest_unoccupied_energy = lowest_unoccupied_energy_by_flavor[flavor_index]

            chemical_potential_by_flavor[flavor_index] =
                (highest_occupied_energy + lowest_unoccupied_energy) / 2

            boundary_gap_by_flavor[flavor_index] =
                lowest_unoccupied_energy - highest_occupied_energy
        end
    end

    return FillingDiagnostics(occupied_state_counts, highest_occupied_energy_by_flavor,
                              lowest_unoccupied_energy_by_flavor, chemical_potential_by_flavor,
                              boundary_gap_by_flavor)
end

# Piece 2: Momentum mesh

struct MomentumMesh
    mesh_size_1::Int
    mesh_size_2::Int
    momentum_count::Int
    reciprocal_vector_1::Vector{Float64}
    reciprocal_vector_2::Vector{Float64}
    momentum_step_1::Vector{Float64}
    momentum_step_2::Vector{Float64}
    unit_cell_area::Float64
    mesh_integer_1::Vector{Int}
    mesh_integer_2::Vector{Int}
    momenta::Matrix{Float64}
end


function MomentumMesh(reciprocal_vector_1::AbstractVector{<:Real}, reciprocal_vector_2::AbstractVector{<:Real},
                      mesh_size_1::Int, mesh_size_2::Int)
    length(reciprocal_vector_1) == 2 ||
        throw(ArgumentError("reciprocal_vector_1 must have two components."))

    length(reciprocal_vector_2) == 2 ||
        throw(ArgumentError("reciprocal_vector_2 must have two components."))

    mesh_size_1 >= 1 ||
        throw(ArgumentError("mesh_size_1 must be positive."))

    mesh_size_2 >= 1 ||
        throw(ArgumentError("mesh_size_2 must be positive."))

    reciprocal_vector_1_float = Float64.(reciprocal_vector_1)
    reciprocal_vector_2_float = Float64.(reciprocal_vector_2)

    reciprocal_cell_area = abs(reciprocal_vector_1_float[1] * reciprocal_vector_2_float[2] -
                               reciprocal_vector_1_float[2] * reciprocal_vector_2_float[1])

    reciprocal_area_scale = norm(reciprocal_vector_1_float) * norm(reciprocal_vector_2_float)
    singularity_tolerance = 100 * eps(Float64) * reciprocal_area_scale

    reciprocal_cell_area > singularity_tolerance ||
        throw(ArgumentError("The reciprocal-lattice vectors must be linearly independent."))

    momentum_step_1 = reciprocal_vector_1_float ./ mesh_size_1
    momentum_step_2 = reciprocal_vector_2_float ./ mesh_size_2
    unit_cell_area = (2π)^2 / reciprocal_cell_area
    momentum_count = mesh_size_1 * mesh_size_2

    mesh_integer_1 = Vector{Int}(undef, momentum_count)
    mesh_integer_2 = Vector{Int}(undef, momentum_count)
    momenta = Matrix{Float64}(undef, 2, momentum_count)

    for mesh_integer_2_value in 0:mesh_size_2-1
        for mesh_integer_1_value in 0:mesh_size_1-1
            k_index = mesh_integer_1_value + 1 + mesh_size_1 * mesh_integer_2_value

            mesh_integer_1[k_index] = mesh_integer_1_value
            mesh_integer_2[k_index] = mesh_integer_2_value

            momenta[1, k_index] = mesh_integer_1_value * momentum_step_1[1] +
                                  mesh_integer_2_value * momentum_step_2[1]

            momenta[2, k_index] = mesh_integer_1_value * momentum_step_1[2] +
                                  mesh_integer_2_value * momentum_step_2[2]
        end
    end

    return MomentumMesh(mesh_size_1, mesh_size_2, momentum_count,
                        reciprocal_vector_1_float, reciprocal_vector_2_float,
                        momentum_step_1, momentum_step_2, unit_cell_area,
                        mesh_integer_1, mesh_integer_2, momenta)
end


@inline function momentum_transfer_mesh_coordinates(momentum_mesh::MomentumMesh, k_index::Int, k2_index::Int,
                                                    folding_difference_1::Int, folding_difference_2::Int)
    @inbounds k_mesh_integer_1 = momentum_mesh.mesh_integer_1[k_index]
    @inbounds k_mesh_integer_2 = momentum_mesh.mesh_integer_2[k_index]
    @inbounds k2_mesh_integer_1 = momentum_mesh.mesh_integer_1[k2_index]
    @inbounds k2_mesh_integer_2 = momentum_mesh.mesh_integer_2[k2_index]

    transfer_mesh_coordinate_1 = k2_mesh_integer_1 - k_mesh_integer_1 +
                                 momentum_mesh.mesh_size_1 * folding_difference_1

    transfer_mesh_coordinate_2 = k2_mesh_integer_2 - k_mesh_integer_2 +
                                 momentum_mesh.mesh_size_2 * folding_difference_2

    return transfer_mesh_coordinate_1, transfer_mesh_coordinate_2
end


# Piece 3: Projected folded basis

struct ProjectedBasis
    folding_coordinates::Matrix{Int}
    selected_folding_index::Array{Int,3}
    spinor::Array{ComplexF64,4}
    energy::Array{Float64,3}
end


function get_folding_coordinates(momentum_mesh::MomentumMesh, reciprocal_cutoff::Real)
    reciprocal_cutoff > 0 ||
        throw(ArgumentError("reciprocal_cutoff must be positive."))

    shortest_reciprocal_length = min(norm(momentum_mesh.reciprocal_vector_1),
                                     norm(momentum_mesh.reciprocal_vector_2))

    integer_search_limit = ceil(Int, 2 * reciprocal_cutoff / shortest_reciprocal_length)
    folding_coordinate_pairs = Tuple{Int,Int}[]

    for folding_coordinate_2 in -integer_search_limit:integer_search_limit
        for folding_coordinate_1 in -integer_search_limit:integer_search_limit
            reciprocal_momentum = folding_coordinate_1 .* momentum_mesh.reciprocal_vector_1 .+
                                  folding_coordinate_2 .* momentum_mesh.reciprocal_vector_2

            if norm(reciprocal_momentum) <= reciprocal_cutoff
                push!(folding_coordinate_pairs, (folding_coordinate_1, folding_coordinate_2))
            end
        end
    end

    sort!(folding_coordinate_pairs;
          by=folding_coordinates -> norm(folding_coordinates[1] .* momentum_mesh.reciprocal_vector_1 .+
                                         folding_coordinates[2] .* momentum_mesh.reciprocal_vector_2))

    folding_coordinates = Matrix{Int}(undef, 2, length(folding_coordinate_pairs))

    for folding_index in eachindex(folding_coordinate_pairs)
        folding_coordinates[1, folding_index] = folding_coordinate_pairs[folding_index][1]
        folding_coordinates[2, folding_index] = folding_coordinate_pairs[folding_index][2]
    end

    return folding_coordinates
end


function get_projected_basis(
    momentum_mesh::MomentumMesh,
    folding_coordinates::Matrix{Int},
    projected_band_count::Int,
    layer_count::Int,
    displacement_field::Real,
    valley_by_flavor::AbstractVector{<:Integer},
    displacement_sign_by_flavor::AbstractVector{<:Integer},
    stacking::Int,
    project_to_valence::Bool
)

    size(folding_coordinates, 1) == 2 ||
        throw(DimensionMismatch("folding_coordinates must have dimensions 2 × folding_vector_count."))

    orbital_count = 2 * layer_count
    flavor_count = length(valley_by_flavor)
    momentum_count = momentum_mesh.momentum_count
    folding_vector_count = size(folding_coordinates, 2)
    candidate_count = layer_count * folding_vector_count

    projected_band_count <= candidate_count ||
        throw(ArgumentError("The reciprocal cutoff does not provide enough candidate conduction states."))

    folding_vector_norm_by_index = Vector{Float64}(undef, folding_vector_count)

    for folding_index in 1:folding_vector_count
        folding_coordinate_1 = folding_coordinates[1, folding_index]
        folding_coordinate_2 = folding_coordinates[2, folding_index]

        folding_momentum_1 = folding_coordinate_1 * momentum_mesh.reciprocal_vector_1[1] +
                             folding_coordinate_2 * momentum_mesh.reciprocal_vector_2[1]

        folding_momentum_2 = folding_coordinate_1 * momentum_mesh.reciprocal_vector_1[2] +
                             folding_coordinate_2 * momentum_mesh.reciprocal_vector_2[2]

        folding_vector_norm_by_index[folding_index] = hypot(folding_momentum_1, folding_momentum_2)
    end

    outermost_available_folding_norm = maximum(folding_vector_norm_by_index)

    # selected_folding_index[alpha_index, k_index, flavor_index]
    # selected_spinor[orbital_index, alpha_index, k_index, flavor_index]
    # selected_energy[alpha_index, k_index, flavor_index]

    selected_folding_index = Array{Int,3}(undef, projected_band_count, momentum_count, flavor_count)
    selected_spinor = Array{ComplexF64,4}(undef, orbital_count, projected_band_count,
                                          momentum_count, flavor_count)
    selected_energy = Array{Float64,3}(undef, projected_band_count, momentum_count, flavor_count)

    candidate_folding_index = Vector{Int}(undef, candidate_count)
    candidate_spinor = Matrix{ComplexF64}(undef, orbital_count, candidate_count)
    candidate_energy = Vector{Float64}(undef, candidate_count)
    candidate_order = Vector{Int}(undef, candidate_count)
    relative_momentum = Vector{Float64}(undef, 2)


        if project_to_valence
            continuum_band_indices = 1:layer_count
        else
            continuum_band_indices = layer_count+1:2*layer_count
        end

    for flavor_index in 1:flavor_count
        valley = Int(valley_by_flavor[flavor_index])

        flavor_displacement_field =
            displacement_field *
            displacement_sign_by_flavor[flavor_index]

        for k_index in 1:momentum_count
            candidate_index = 0

            for folding_index in 1:folding_vector_count
                folding_coordinate_1 = folding_coordinates[1, folding_index]
                folding_coordinate_2 = folding_coordinates[2, folding_index]

                relative_momentum[1] = momentum_mesh.momenta[1, k_index] +
                                       folding_coordinate_1 * momentum_mesh.reciprocal_vector_1[1] +
                                       folding_coordinate_2 * momentum_mesh.reciprocal_vector_2[1]

                relative_momentum[2] = momentum_mesh.momenta[2, k_index] +
                                       folding_coordinate_1 * momentum_mesh.reciprocal_vector_1[2] +
                                       folding_coordinate_2 * momentum_mesh.reciprocal_vector_2[2]

                continuum_eigen = eigen(Hermitian(get_Ham(relative_momentum,Float64(flavor_displacement_field),
                        valley,
                        stacking,
                        layer_count)))

               for continuum_band_index in continuum_band_indices
                    candidate_index += 1
                    candidate_folding_index[candidate_index] = folding_index
                    candidate_energy[candidate_index] = continuum_eigen.values[continuum_band_index]
                    candidate_spinor[:, candidate_index] .= continuum_eigen.vectors[:, continuum_band_index]
                end
            end

            if project_to_valence
                sortperm!(
                    candidate_order,
                    candidate_energy;
                    alg=Base.Sort.MergeSort,
                    rev=true
                )
            else
                sortperm!(
                    candidate_order,
                    candidate_energy;
                    alg=Base.Sort.MergeSort
                )
            end

            for projected_band_index in 1:projected_band_count
                selected_candidate_index = candidate_order[projected_band_index]
                selected_folding_index[projected_band_index, k_index, flavor_index] =
                    candidate_folding_index[selected_candidate_index]
                selected_energy[projected_band_index, k_index, flavor_index] =
                    candidate_energy[selected_candidate_index]
                selected_spinor[:, projected_band_index, k_index, flavor_index] .=
                    candidate_spinor[:, selected_candidate_index]
            end
        end
    end

    largest_selected_folding_norm = 0.0

    for selected_folding_list_index in selected_folding_index
        selected_folding_norm = folding_vector_norm_by_index[selected_folding_list_index]
        largest_selected_folding_norm = max(largest_selected_folding_norm, selected_folding_norm)
    end

    reciprocal_norm_scale = max(outermost_available_folding_norm,
                                norm(momentum_mesh.reciprocal_vector_1),
                                norm(momentum_mesh.reciprocal_vector_2))

    shell_tolerance = 100 * eps(Float64) * reciprocal_norm_scale
    outer_shell_clearance = outermost_available_folding_norm - largest_selected_folding_norm

    outer_shell_clearance > shell_tolerance ||
        throw(ArgumentError(
            "The projected basis reaches the outermost available folding-vector shell: " *
            "largest selected |G| = $largest_selected_folding_norm, " *
            "outermost available |G| = $outermost_available_folding_norm. " *
            "Increase reciprocal_cutoff and rebuild folding_coordinates."
        ))

    return ProjectedBasis(folding_coordinates, selected_folding_index, selected_spinor, selected_energy)
end


# Piece 4: Single-particle Hamiltonian and screened Coulomb interaction

const COULOMB_PREFACTOR_MEV_NM = 9047.5636


function get_single_particle_hamiltonian(projected_basis::ProjectedBasis)
    projected_band_count, momentum_count, flavor_count = size(projected_basis.energy)
    hartree_fock_dimension = projected_band_count * flavor_count

    single_particle_hamiltonian =
        zeros(ComplexF64, hartree_fock_dimension, hartree_fock_dimension, momentum_count)

    for flavor_index in 1:flavor_count
        for k_index in 1:momentum_count
            for projected_band_index in 1:projected_band_count
                state_index = projected_band_flavor_index(projected_band_index, flavor_index,
                                                           projected_band_count)

                single_particle_hamiltonian[state_index, state_index, k_index] =
                    projected_basis.energy[projected_band_index, k_index, flavor_index]
            end
        end
    end

    return single_particle_hamiltonian
end


function coulomb_interaction(momentum_mesh::MomentumMesh, transfer_mesh_coordinate_1::Int,
                             transfer_mesh_coordinate_2::Int, dielectric_constant::Real,
                             gate_distance::Real)
    dielectric_constant > 0 ||
        throw(ArgumentError("dielectric_constant must be positive."))

    gate_distance >= 0 ||
        throw(ArgumentError("gate_distance cannot be negative."))

    momentum_transfer_1 = transfer_mesh_coordinate_1 * momentum_mesh.momentum_step_1[1] +
                          transfer_mesh_coordinate_2 * momentum_mesh.momentum_step_2[1]

    momentum_transfer_2 = transfer_mesh_coordinate_1 * momentum_mesh.momentum_step_1[2] +
                          transfer_mesh_coordinate_2 * momentum_mesh.momentum_step_2[2]

    momentum_transfer_norm = sqrt(momentum_transfer_1^2 + momentum_transfer_2^2)

    if transfer_mesh_coordinate_1 == 0 && transfer_mesh_coordinate_2 == 0
        return COULOMB_PREFACTOR_MEV_NM * gate_distance / dielectric_constant
    end

    return COULOMB_PREFACTOR_MEV_NM *
           tanh(momentum_transfer_norm * gate_distance) /
           (dielectric_constant * momentum_transfer_norm)
end


# Piece 5: Form factors

struct FormFactorWork
    overlaps_for_k::Matrix{ComplexF64}
end

function FormFactorWork(projected_basis::ProjectedBasis)
    projected_band_count, momentum_count, flavor_count = size(projected_basis.energy)
    overlaps_for_k = zeros(ComplexF64, projected_band_count, projected_band_count * momentum_count)
    return FormFactorWork(overlaps_for_k)
end


struct FormFactorThreadWorkspacePool
    form_factor_work_by_thread::Vector{FormFactorWork}
end

function FormFactorThreadWorkspacePool(projected_basis::ProjectedBasis)
    maximum_thread_id = Threads.maxthreadid()
    form_factor_work_by_thread = [FormFactorWork(projected_basis) for thread_index in 1:maximum_thread_id]
    return FormFactorThreadWorkspacePool(form_factor_work_by_thread)
end


function calculate_form_factors_for_k!(form_factor_work::FormFactorWork, projected_basis::ProjectedBasis,
                                       k_index::Int, flavor_index::Int)
    orbital_count, projected_band_count, momentum_count, flavor_count = size(projected_basis.spinor)

    spinors_at_k = @view projected_basis.spinor[:, :, k_index, flavor_index]
    spinors_at_all_k2 = reshape(@view(projected_basis.spinor[:, :, :, flavor_index]),
                                orbital_count, projected_band_count * momentum_count)

    mul!(form_factor_work.overlaps_for_k, adjoint(spinors_at_k), spinors_at_all_k2)
    return nothing
end


function precompute_form_factors(projected_basis::ProjectedBasis)
    projected_band_count, momentum_count, flavor_count = size(projected_basis.energy)

    # form_factors[alpha_index, gamma_index, k2_index, k_index, flavor_index]
    #     = <u_{flavor,alpha}(k) | u_{flavor,gamma}(k2)>
    form_factors = Array{ComplexF64,5}(undef, projected_band_count, projected_band_count,
                                      momentum_count, momentum_count, flavor_count)

    form_factor_thread_workspace_pool = FormFactorThreadWorkspacePool(projected_basis)
    form_factor_calculation_count = momentum_count * flavor_count

    Threads.@threads :static for calculation_index in 1:form_factor_calculation_count
        k_index = mod1(calculation_index, momentum_count)
        flavor_index = fld(calculation_index - 1, momentum_count) + 1
        current_thread_id = Threads.threadid()

        thread_form_factor_work =
            form_factor_thread_workspace_pool.form_factor_work_by_thread[current_thread_id]

        calculate_form_factors_for_k!(thread_form_factor_work, projected_basis, k_index, flavor_index)

        stored_form_factors_for_k =
            reshape(@view(form_factors[:, :, :, k_index, flavor_index]),
                    projected_band_count, projected_band_count * momentum_count)

        copyto!(stored_form_factors_for_k, thread_form_factor_work.overlaps_for_k)
    end

    return form_factors
end

# Piece 6: Hartree Hamiltonian

@inline function selected_folding_coordinates(projected_basis::ProjectedBasis, projected_band_index::Int,
                                              k_index::Int, flavor_index::Int)
    folding_list_index = projected_basis.selected_folding_index[projected_band_index, k_index, flavor_index]

    return projected_basis.folding_coordinates[1, folding_list_index],
           projected_basis.folding_coordinates[2, folding_list_index]
end


struct HartreeWork
    reciprocal_vector_coordinates::Vector{Tuple{Int,Int}}
    reciprocal_vector_index::Array{Int,4}
    interaction_by_reciprocal_vector::Vector{Float64}
    density_by_reciprocal_vector::Vector{ComplexF64}
end


function HartreeWork(projected_basis::ProjectedBasis, momentum_mesh::MomentumMesh,
                     dielectric_constant::Real, gate_distance::Real)
    projected_band_count, momentum_count, flavor_count = size(projected_basis.energy)

    reciprocal_vector_index =
        Array{Int,4}(undef, projected_band_count, projected_band_count, momentum_count, flavor_count)

    reciprocal_vector_lookup = Dict{Tuple{Int,Int},Int}()
    reciprocal_vector_coordinates = Tuple{Int,Int}[]

    for flavor_index in 1:flavor_count, k_index in 1:momentum_count,
        alpha_index in 1:projected_band_count, beta_index in 1:projected_band_count

        alpha_folding_coordinates =
            selected_folding_coordinates(projected_basis, alpha_index, k_index, flavor_index)

        beta_folding_coordinates =
            selected_folding_coordinates(projected_basis, beta_index, k_index, flavor_index)

        hartree_reciprocal_coordinates =
            (beta_folding_coordinates[1] - alpha_folding_coordinates[1],
             beta_folding_coordinates[2] - alpha_folding_coordinates[2])

        if haskey(reciprocal_vector_lookup, hartree_reciprocal_coordinates)
            hartree_vector_index = reciprocal_vector_lookup[hartree_reciprocal_coordinates]
        else
            push!(reciprocal_vector_coordinates, hartree_reciprocal_coordinates)
            hartree_vector_index = length(reciprocal_vector_coordinates)
            reciprocal_vector_lookup[hartree_reciprocal_coordinates] = hartree_vector_index
        end

        reciprocal_vector_index[alpha_index, beta_index, k_index, flavor_index] =
            hartree_vector_index
    end

    reciprocal_vector_count = length(reciprocal_vector_coordinates)
    interaction_by_reciprocal_vector = Vector{Float64}(undef, reciprocal_vector_count)

    for hartree_vector_index in 1:reciprocal_vector_count
        hartree_reciprocal_coordinates = reciprocal_vector_coordinates[hartree_vector_index]

        transfer_mesh_coordinate_1 =
            momentum_mesh.mesh_size_1 * hartree_reciprocal_coordinates[1]

        transfer_mesh_coordinate_2 =
            momentum_mesh.mesh_size_2 * hartree_reciprocal_coordinates[2]

        interaction_by_reciprocal_vector[hartree_vector_index] =
            coulomb_interaction(momentum_mesh, transfer_mesh_coordinate_1, transfer_mesh_coordinate_2,
                                dielectric_constant, gate_distance)
    end

    density_by_reciprocal_vector = zeros(ComplexF64, reciprocal_vector_count)

    return HartreeWork(reciprocal_vector_coordinates, reciprocal_vector_index,
                       interaction_by_reciprocal_vector, density_by_reciprocal_vector)
end

function construct_hartree_hamiltonian!(hartree_hamiltonian,density_matrix,form_factors,
                                        hartree_work,momentum_mesh,reference_occupation)
    projected_band_count = size(form_factors, 1)
    momentum_count = size(form_factors, 3)
    flavor_count = size(form_factors, 5)
    inverse_total_area = 1 / (momentum_count * momentum_mesh.unit_cell_area)

    fill!(hartree_hamiltonian, 0)
    fill!(hartree_work.density_by_reciprocal_vector, 0)

    for flavor_index in 1:flavor_count, k_index in 1:momentum_count,
        alpha_index in 1:projected_band_count, beta_index in 1:projected_band_count

        alpha_f_index =
            projected_band_flavor_index(alpha_index, flavor_index, projected_band_count)

        beta_f_index =
            projected_band_flavor_index(beta_index, flavor_index, projected_band_count)

        hartree_vector_index =
            hartree_work.reciprocal_vector_index[alpha_index, beta_index, k_index, flavor_index]

        form_factor = form_factors[alpha_index, beta_index, k_index, k_index, flavor_index]

        relative_density_element=density_matrix[alpha_f_index,beta_f_index,k_index]

        if alpha_f_index==beta_f_index
            relative_density_element-=reference_occupation
        end

        hartree_work.density_by_reciprocal_vector[hartree_vector_index] +=
            conj(form_factor)*relative_density_element
    end

    for flavor_index in 1:flavor_count, k_index in 1:momentum_count,
        alpha_index in 1:projected_band_count, beta_index in 1:projected_band_count

        alpha_f_index =
            projected_band_flavor_index(alpha_index, flavor_index, projected_band_count)

        beta_f_index =
            projected_band_flavor_index(beta_index, flavor_index, projected_band_count)

        hartree_vector_index =
            hartree_work.reciprocal_vector_index[alpha_index, beta_index, k_index, flavor_index]

        form_factor = form_factors[alpha_index, beta_index, k_index, k_index, flavor_index]

        hartree_hamiltonian[alpha_f_index, beta_f_index, k_index] =
            hartree_work.interaction_by_reciprocal_vector[hartree_vector_index] *
            hartree_work.density_by_reciprocal_vector[hartree_vector_index] *
            form_factor * inverse_total_area
    end

    return nothing
end

# Piece 7: Fock Hamiltonian

struct FockWork
    folding_difference_coordinates::Vector{Tuple{Int,Int}}
    folding_difference_index_by_pair::Matrix{Int}
    form_factor_index_pairs::Matrix{Vector{Tuple{Int,Int}}}
    active_folding_difference_indices::Vector{Int}
    folding_difference_is_active::BitVector
end


function FockWork(projected_basis::ProjectedBasis)
    flavor_count = size(projected_basis.energy, 3)
    folding_vector_count = size(projected_basis.folding_coordinates, 2)

    folding_difference_coordinates = Tuple{Int,Int}[]
    folding_difference_lookup = Dict{Tuple{Int,Int},Int}()

    folding_difference_index_by_pair =
        Matrix{Int}(undef, folding_vector_count, folding_vector_count)

    for bra_folding_index in 1:folding_vector_count, ket_folding_index in 1:folding_vector_count
        folding_difference =
            (projected_basis.folding_coordinates[1, ket_folding_index] -
             projected_basis.folding_coordinates[1, bra_folding_index],
             projected_basis.folding_coordinates[2, ket_folding_index] -
             projected_basis.folding_coordinates[2, bra_folding_index])

        if haskey(folding_difference_lookup, folding_difference)
            folding_difference_index = folding_difference_lookup[folding_difference]
        else
            push!(folding_difference_coordinates, folding_difference)
            folding_difference_index = length(folding_difference_coordinates)
            folding_difference_lookup[folding_difference] = folding_difference_index
        end

        folding_difference_index_by_pair[bra_folding_index, ket_folding_index] =
            folding_difference_index
    end

    folding_difference_count = length(folding_difference_coordinates)

    form_factor_index_pairs =
        Matrix{Vector{Tuple{Int,Int}}}(undef, folding_difference_count, flavor_count)

    for folding_difference_index in 1:folding_difference_count, flavor_index in 1:flavor_count
        form_factor_index_pairs[folding_difference_index, flavor_index] = Tuple{Int,Int}[]
    end

    active_folding_difference_indices = Int[]
    folding_difference_is_active = falses(folding_difference_count)

    return FockWork(folding_difference_coordinates, folding_difference_index_by_pair,
                    form_factor_index_pairs, active_folding_difference_indices,
                    folding_difference_is_active)
end



struct FockThreadWorkspacePool
    fock_work_by_thread::Vector{FockWork}
end


function FockThreadWorkspacePool(projected_basis::ProjectedBasis)
    maximum_thread_id = Threads.maxthreadid()

    fock_work_by_thread =
        Vector{FockWork}(undef, maximum_thread_id)

    for thread_index in 1:maximum_thread_id
        fock_work_by_thread[thread_index] =
            FockWork(projected_basis)
    end

    return FockThreadWorkspacePool(fock_work_by_thread)
end


function prepare_fock_form_factor_pairs!(fock_work::FockWork, projected_basis::ProjectedBasis,
                                         k_index::Int, k2_index::Int)
    projected_band_count = size(projected_basis.energy, 1)
    flavor_count = size(projected_basis.energy, 3)

    for folding_difference_index in fock_work.active_folding_difference_indices
        fock_work.folding_difference_is_active[folding_difference_index] = false

        for flavor_index in 1:flavor_count
            empty!(fock_work.form_factor_index_pairs[folding_difference_index, flavor_index])
        end
    end

    empty!(fock_work.active_folding_difference_indices)

    for f_index in 1:flavor_count, alpha_index in 1:projected_band_count,
        gamma_index in 1:projected_band_count

        alpha_folding_index =
            projected_basis.selected_folding_index[alpha_index, k_index, f_index]

        gamma_folding_index =
            projected_basis.selected_folding_index[gamma_index, k2_index, f_index]

        folding_difference_index =
            fock_work.folding_difference_index_by_pair[alpha_folding_index, gamma_folding_index]

        if !fock_work.folding_difference_is_active[folding_difference_index]
            fock_work.folding_difference_is_active[folding_difference_index] = true
            push!(fock_work.active_folding_difference_indices, folding_difference_index)
        end

        push!(fock_work.form_factor_index_pairs[folding_difference_index, f_index],
              (alpha_index, gamma_index))
    end

    return nothing
end






function construct_fock_hamiltonian!(fock_hamiltonian, density_matrix, form_factors,
                                     fock_thread_workspace_pool, projected_basis, momentum_mesh,
                                     dielectric_constant, gate_distance, reference_occupation)

    projected_band_count, momentum_count, flavor_count = size(projected_basis.energy)
    hartree_fock_dimension = projected_band_count * flavor_count
    inverse_total_area = 1 / (momentum_count * momentum_mesh.unit_cell_area)

    fill!(fock_hamiltonian, 0)

    Threads.@threads :static for k_index in 1:momentum_count
        current_thread_id = Threads.threadid()
        thread_fock_work = fock_thread_workspace_pool.fock_work_by_thread[current_thread_id]

        for k2_index in 1:momentum_count
            prepare_fock_form_factor_pairs!(thread_fock_work, projected_basis, k_index, k2_index)

            for folding_difference_index in thread_fock_work.active_folding_difference_indices
                folding_difference =
                    thread_fock_work.folding_difference_coordinates[folding_difference_index]

                transfer_mesh_coordinate_1, transfer_mesh_coordinate_2 =
                    momentum_transfer_mesh_coordinates(momentum_mesh, k_index, k2_index,
                                                       folding_difference[1], folding_difference[2])

                coulomb_over_area =
                    inverse_total_area * coulomb_interaction(momentum_mesh,
                                                             transfer_mesh_coordinate_1,
                                                             transfer_mesh_coordinate_2,
                                                             dielectric_constant, gate_distance)

                for flavor_index in 1:flavor_count
                    form_factor_index_pairs =
                        thread_fock_work.form_factor_index_pairs[
                            folding_difference_index, flavor_index
                        ]

                    isempty(form_factor_index_pairs) && continue

                    for (alpha_index, gamma_index) in form_factor_index_pairs,
                        (delta_index, beta_index) in form_factor_index_pairs

                        alpha_f_index =
                            projected_band_flavor_index(alpha_index, flavor_index,
                                                        projected_band_count)

                        delta_f_index =
                            projected_band_flavor_index(delta_index, flavor_index,
                                                        projected_band_count)

                        alpha_f_index > delta_f_index && continue

                        gamma_f_index =
                            projected_band_flavor_index(gamma_index, flavor_index,
                                                        projected_band_count)

                        beta_f_index =
                            projected_band_flavor_index(beta_index, flavor_index,
                                                        projected_band_count)

                        left_form_factor =
                            form_factors[alpha_index, gamma_index, k2_index,
                                         k_index, flavor_index]

                        right_form_factor =
                            form_factors[delta_index, beta_index, k2_index,
                                         k_index, flavor_index]

                        density_element =
                            density_matrix[gamma_f_index, beta_f_index, k2_index]

                        if gamma_index == beta_index
                            density_element -= reference_occupation
                        end

                        fock_hamiltonian[alpha_f_index, delta_f_index, k_index] -=
                            coulomb_over_area * left_form_factor *
                            density_element * conj(right_form_factor)
                    end
                end
            end
        end

        for column_state_index in 1:hartree_fock_dimension
            fock_hamiltonian[column_state_index, column_state_index, k_index] =
                real(fock_hamiltonian[column_state_index, column_state_index, k_index])

            for row_state_index in 1:column_state_index-1
                fock_hamiltonian[column_state_index, row_state_index, k_index] =
                    conj(fock_hamiltonian[row_state_index, column_state_index, k_index])
            end
        end
    end

    return nothing
end

# Piece 8: Full Hartree-Fock Hamiltonian

struct HartreeFockHamiltonianWork
    hartree_hamiltonian::Array{ComplexF64,3}
    fock_hamiltonian::Array{ComplexF64,3}
    hartree_fock_hamiltonian::Array{ComplexF64,3}
end

function HartreeFockHamiltonianWork(projected_basis::ProjectedBasis)
    projected_band_count, momentum_count, flavor_count = size(projected_basis.energy)
    hartree_fock_dimension = projected_band_count * flavor_count
    hamiltonian_size = (hartree_fock_dimension, hartree_fock_dimension, momentum_count)

    hartree_hamiltonian = zeros(ComplexF64, hamiltonian_size)
    fock_hamiltonian = zeros(ComplexF64, hamiltonian_size)
    hartree_fock_hamiltonian = zeros(ComplexF64, hamiltonian_size)

    return HartreeFockHamiltonianWork(hartree_hamiltonian, fock_hamiltonian, hartree_fock_hamiltonian)
end


function construct_hartree_fock_hamiltonian!(hamiltonian_work,density_matrix,
                                             single_particle_hamiltonian,form_factors,
                                             hartree_work,fock_thread_workspace_pool,
                                             projected_basis,momentum_mesh,
                                             dielectric_constant,gate_distance,reference_occupation)
    construct_hartree_hamiltonian!(hamiltonian_work.hartree_hamiltonian,density_matrix,
                               form_factors,hartree_work,momentum_mesh,reference_occupation)

    construct_fock_hamiltonian!(hamiltonian_work.fock_hamiltonian,density_matrix,form_factors,
                            fock_thread_workspace_pool,projected_basis,momentum_mesh,
                            dielectric_constant,gate_distance,reference_occupation)

    copyto!(
        hamiltonian_work.hartree_fock_hamiltonian,
        single_particle_hamiltonian
    )

    hamiltonian_work.hartree_fock_hamiltonian .+=
        hamiltonian_work.hartree_hamiltonian

    hamiltonian_work.hartree_fock_hamiltonian .+=
        hamiltonian_work.fock_hamiltonian

    return nothing
end


# Piece 9: Optimal damping algorithm

struct OptimalDampingWork
    trial_hamiltonian_work::HartreeFockHamiltonianWork
    trial_density_difference::Array{ComplexF64,3}
end

function OptimalDampingWork(projected_basis::ProjectedBasis)
    projected_band_count, momentum_count, flavor_count = size(projected_basis.energy)
    hartree_fock_dimension = projected_band_count * flavor_count

    trial_hamiltonian_work = HartreeFockHamiltonianWork(projected_basis)
    trial_density_difference = zeros(ComplexF64, hartree_fock_dimension, hartree_fock_dimension, momentum_count)

    return OptimalDampingWork(trial_hamiltonian_work, trial_density_difference)
end

function hartree_fock_energy(density_matrix,single_particle_hamiltonian,
                             hamiltonian_work,reference_occupation)
    relative_density_matrix=copy(density_matrix)
    hartree_fock_dimension,_,momentum_count=size(relative_density_matrix)

    for k_index in 1:momentum_count, state_index in 1:hartree_fock_dimension
        relative_density_matrix[state_index,state_index,k_index]-=reference_occupation
    end

    single_particle_energy=real(dot(single_particle_hamiltonian,relative_density_matrix))
    hartree_energy=0.5*real(dot(hamiltonian_work.hartree_hamiltonian,relative_density_matrix))
    fock_energy=0.5*real(dot(hamiltonian_work.fock_hamiltonian,relative_density_matrix))

    return single_particle_energy+hartree_energy+fock_energy
end

function calculate_oda_energy_coefficients!(optimal_damping_work, current_density_matrix,
                                            trial_density_matrix, current_hamiltonian_work)
    trial_hamiltonian_work = optimal_damping_work.trial_hamiltonian_work
    trial_density_difference = optimal_damping_work.trial_density_difference

    copyto!(trial_density_difference, trial_density_matrix)
    trial_density_difference .-= current_density_matrix

    energy_linear_coefficient = real(dot(current_hamiltonian_work.hartree_fock_hamiltonian, trial_density_difference))

    hartree_energy_curvature = real(dot(trial_hamiltonian_work.hartree_hamiltonian, trial_density_difference) -
                                    dot(current_hamiltonian_work.hartree_hamiltonian, trial_density_difference))

    fock_energy_curvature = real(dot(trial_hamiltonian_work.fock_hamiltonian, trial_density_difference) -
                                 dot(current_hamiltonian_work.fock_hamiltonian, trial_density_difference))

    energy_quadratic_coefficient = 0.5 * (hartree_energy_curvature + fock_energy_curvature)

    return energy_linear_coefficient, energy_quadratic_coefficient
end

function optimal_damping_parameter(energy_linear_coefficient, energy_quadratic_coefficient)
    if energy_quadratic_coefficient > 0
        return clamp(-energy_linear_coefficient / (2 * energy_quadratic_coefficient), 0.0, 1.0)
    end

    full_step_energy_change = energy_linear_coefficient + energy_quadratic_coefficient
    return full_step_energy_change < 0 ? 1.0 : 0.0
end

function apply_oda_update!(current_density_matrix, current_hamiltonian_work,
                           optimal_damping_work, single_particle_hamiltonian, damping_parameter)
    trial_hamiltonian_work = optimal_damping_work.trial_hamiltonian_work
    trial_density_difference = optimal_damping_work.trial_density_difference

    current_density_matrix .+= damping_parameter .* trial_density_difference

    current_hamiltonian_work.hartree_hamiltonian .*= 1 - damping_parameter
    current_hamiltonian_work.hartree_hamiltonian .+= damping_parameter .* trial_hamiltonian_work.hartree_hamiltonian

    current_hamiltonian_work.fock_hamiltonian .*= 1 - damping_parameter
    current_hamiltonian_work.fock_hamiltonian .+= damping_parameter .* trial_hamiltonian_work.fock_hamiltonian

    copyto!(current_hamiltonian_work.hartree_fock_hamiltonian, single_particle_hamiltonian)
    current_hamiltonian_work.hartree_fock_hamiltonian .+= current_hamiltonian_work.hartree_hamiltonian
    current_hamiltonian_work.hartree_fock_hamiltonian .+= current_hamiltonian_work.fock_hamiltonian

    return nothing
end

# Piece 10: Self-consistent Hartree-Fock loop with ODA

function solve_hartree_fock_with_oda!(
    density_matrix,
    filling_work,
    current_hamiltonian_work,
    optimal_damping_work,
    single_particle_hamiltonian,
    form_factors,
    hartree_work,
    fock_thread_workspace_pool,
    projected_basis,
    momentum_mesh,
    dielectric_constant,
    gate_distance,
    occupied_state_count_by_flavor;
    reference_occupation=0,
    maximum_iterations=500,
    density_tolerance=1e-7,
    energy_tolerance=1e-8,
    verbose=true
)
    maximum_iterations >= 1 ||
        throw(ArgumentError("maximum_iterations must be positive."))
    total_occupied_state_count = sum(occupied_state_count_by_flavor)
    reference_occupied_state_count = reference_occupation * length(filling_work.eigenvalues)
    carrier_count = abs(total_occupied_state_count - reference_occupied_state_count)
    carrier_count > 0 || throw(ArgumentError("carrier_count must be positive."))

    construct_hartree_fock_hamiltonian!(
        current_hamiltonian_work,
        density_matrix,
        single_particle_hamiltonian,
        form_factors,
        hartree_work,
        fock_thread_workspace_pool,
        projected_basis,
        momentum_mesh,
        dielectric_constant,
        gate_distance,
        reference_occupation
    )

    current_energy =
       hartree_fock_energy(density_matrix,single_particle_hamiltonian,
                    current_hamiltonian_work,reference_occupation)

    has_converged = false
    iteration_count = 0
    density_residual = Inf
    energy_change_per_carrier = Inf
    damping_parameter = 0.0
    filling_diagnostics = nothing

    for iteration_index in 1:maximum_iterations
        iteration_start=time()
        iteration_count = iteration_index
      
        diagonalize_hartree_fock_hamiltonian!(
            filling_work,
            current_hamiltonian_work.hartree_fock_hamiltonian
        )

        filling_diagnostics =
        construct_zero_temperature_density_matrix!(
            filling_work,
            occupied_state_count_by_flavor
        )
        trial_density_matrix =
            filling_work.trial_density_matrix

        construct_hartree_fock_hamiltonian!(
            optimal_damping_work.trial_hamiltonian_work,
            trial_density_matrix,
            single_particle_hamiltonian,
            form_factors,
            hartree_work,
            fock_thread_workspace_pool,
            projected_basis,
            momentum_mesh,
            dielectric_constant,
            gate_distance,
            reference_occupation
        )

        energy_linear_coefficient, energy_quadratic_coefficient =
            calculate_oda_energy_coefficients!(
                optimal_damping_work,
                density_matrix,
                trial_density_matrix,
                current_hamiltonian_work
            )

        density_residual=sum(abs2,optimal_damping_work.trial_density_difference)/momentum_mesh.momentum_count

        oda_parameter = optimal_damping_parameter(energy_linear_coefficient,
                                                energy_quadratic_coefficient)

        damping_parameter = oda_parameter < 1e-8 ? 0.02 : oda_parameter

        previous_energy = current_energy

        apply_oda_update!(
            density_matrix,
            current_hamiltonian_work,
            optimal_damping_work,
            single_particle_hamiltonian,
            damping_parameter
        )

        current_energy =
            hartree_fock_energy(density_matrix,single_particle_hamiltonian,
                    current_hamiltonian_work,reference_occupation)

        energy_change_per_carrier=abs(current_energy-previous_energy)/carrier_count

        iteration_time=time()-iteration_start
        if verbose
            println("iteration ",iteration_index,
            ": energy/carrier = ",current_energy/carrier_count,
            ", energy change/carrier = ",energy_change_per_carrier,
            ", density residual = ",density_residual,
            ", ODA parameter = ",damping_parameter,
            ", time = ",iteration_time)
        end

       if density_residual<=density_tolerance && energy_change_per_carrier<=energy_tolerance
            has_converged = true
            break
        end
    end

    has_converged ||
    error("Hartree-Fock calculation did not converge after $maximum_iterations iterations.")

    return (
        converged=has_converged,
        iteration_count=iteration_count,
        energy=current_energy,
        energy_per_carrier=current_energy/carrier_count,
        density_residual=density_residual,
        energy_change_per_carrier=energy_change_per_carrier,
        damping_parameter=damping_parameter,
        filling_diagnostics=filling_diagnostics
    )
end


# Piece 11: Initialization and driver

function initialize_density_matrix!(density_matrix::Array{ComplexF64,3},
                                    filling_work::HartreeFockFillingWork,
                                    single_particle_hamiltonian::Array{ComplexF64,3},
                                    occupied_state_count_by_flavor::AbstractVector{<:Integer})

    random_initialization_hamiltonian = zeros(ComplexF64, size(single_particle_hamiltonian))

    projected_band_count = filling_work.projected_band_count
    flavor_count = filling_work.flavor_count
    momentum_count = size(single_particle_hamiltonian, 3)

    random_complex_matrix = zeros(ComplexF64, projected_band_count, projected_band_count)
    random_hermitian_matrix = similar(random_complex_matrix)

    for flavor_index in 1:flavor_count, k_index in 1:momentum_count
        first_state_index =
            projected_band_flavor_index(1, flavor_index, projected_band_count)

        last_state_index =
            projected_band_flavor_index(projected_band_count, flavor_index, projected_band_count)

        flavor_state_indices = first_state_index:last_state_index

        randn!(random_complex_matrix)
        random_hermitian_matrix .=
            0.5 .* (random_complex_matrix .+ adjoint(random_complex_matrix))

        @views random_initialization_hamiltonian[
            flavor_state_indices, flavor_state_indices, k_index
        ] .= random_hermitian_matrix
    end

    diagonalize_hartree_fock_hamiltonian!(filling_work, random_initialization_hamiltonian)
    construct_zero_temperature_density_matrix!(filling_work, occupied_state_count_by_flavor)
    copyto!(density_matrix, filling_work.trial_density_matrix)

    return nothing
end

function initialize_density_matrix_from_seed!(density_matrix::Array{ComplexF64,3}, 
    projected_basis::ProjectedBasis, occupied_state_count_by_flavor::AbstractVector{<:Integer},
     seed_folder::Union{Nothing,AbstractString}, carrier_population_code::Int, 
     displacement_sign_code::Int, valley_by_flavor::AbstractVector{<:Integer}; reference_occupation::Int=0)   

    seed_folder === nothing && return false
    isdir(seed_folder) || return false

    seed_paths = filter(path -> isfile(path) && endswith(lowercase(path), ".jld2"), readdir(seed_folder; join=true))
    isempty(seed_paths) && return false

    seed_path = rand(seed_paths)
    seed_data = JLD2.load(seed_path)

    seed_density_matrix = seed_data["densitymatrix"]
    seed_spinor = seed_data["single_eigenvector"]
    seed_selected_folding_index = seed_data["selected_folding_index"]
    seed_folding_coordinates = seed_data["folding_coordinates"]
    seed_reference_occupation = Int(seed_data["reference_occupation"])

    seed_reference_occupation == reference_occupation || error("The seed uses a different reference_occupation.")

    seed_orbital_count, seed_band_count, seed_momentum_count, seed_flavor_count = size(seed_spinor)
    orbital_count, projected_band_count, momentum_count, flavor_count = size(projected_basis.spinor)

    seed_orbital_count == orbital_count || throw(DimensionMismatch("The seed and current runs have different NL."))
    seed_momentum_count == momentum_count || throw(DimensionMismatch("The seed and current runs have different Nq."))
    seed_flavor_count == flavor_count || throw(DimensionMismatch("The seed and current runs have different flavor counts."))
    seed_data["carrier_population_code"] == carrier_population_code || error("The seed uses a different carrier_population_code.")
    seed_data["displacement_sign_code"] == displacement_sign_code || error("The seed uses a different displacement_sign_code.")
    seed_data["valley_by_flavor"] == valley_by_flavor || error("The seed uses a different valley_by_flavor.")

    expected_seed_density_size = (seed_band_count * flavor_count, seed_band_count * flavor_count, momentum_count)
    size(seed_density_matrix) == expected_seed_density_size || throw(DimensionMismatch("The saved density matrix has inconsistent dimensions."))

    fill!(density_matrix, 0)

    basis_overlap = zeros(ComplexF64, projected_band_count, seed_band_count)
    current_relative_count_by_flavor = zeros(Float64, flavor_count)

    for flavor_index in 1:flavor_count, k_index in 1:momentum_count
        fill!(basis_overlap, 0)

        for projected_band_index in 1:projected_band_count, seed_band_index in 1:seed_band_count
            current_folding_index = projected_basis.selected_folding_index[projected_band_index, k_index, flavor_index]
            seed_folding_index = seed_selected_folding_index[seed_band_index, k_index, flavor_index]

            same_folding_vector = projected_basis.folding_coordinates[1, current_folding_index] == seed_folding_coordinates[1, seed_folding_index] &&
                                  projected_basis.folding_coordinates[2, current_folding_index] == seed_folding_coordinates[2, seed_folding_index]

            if same_folding_vector
                current_spinor = @view projected_basis.spinor[:, projected_band_index, k_index, flavor_index]
                old_spinor = @view seed_spinor[:, seed_band_index, k_index, flavor_index]
                basis_overlap[projected_band_index, seed_band_index] = dot(current_spinor, old_spinor)
            end
        end

        seed_first_state_index = projected_band_flavor_index(1, flavor_index, seed_band_count)
        seed_last_state_index = projected_band_flavor_index(seed_band_count, flavor_index, seed_band_count)
        seed_flavor_state_indices = seed_first_state_index:seed_last_state_index

        seed_relative_density_matrix = copy(@view seed_density_matrix[seed_flavor_state_indices, seed_flavor_state_indices, k_index])

        for seed_band_index in 1:seed_band_count
            seed_relative_density_matrix[seed_band_index, seed_band_index] -= reference_occupation
        end

        mapped_relative_density_matrix = basis_overlap * seed_relative_density_matrix * adjoint(basis_overlap)
        mapped_relative_density_matrix = 0.5 .* (mapped_relative_density_matrix + adjoint(mapped_relative_density_matrix))

        first_state_index = projected_band_flavor_index(1, flavor_index, projected_band_count)
        last_state_index = projected_band_flavor_index(projected_band_count, flavor_index, projected_band_count)
        flavor_state_indices = first_state_index:last_state_index
        density_matrix_block = @view density_matrix[flavor_state_indices, flavor_state_indices, k_index]

        copyto!(density_matrix_block, mapped_relative_density_matrix)
        current_relative_count_by_flavor[flavor_index] += real(tr(mapped_relative_density_matrix))
    end

    states_per_flavor = projected_band_count * momentum_count

    for flavor_index in 1:flavor_count
        target_relative_count = occupied_state_count_by_flavor[flavor_index] - reference_occupation * states_per_flavor
        current_relative_count = current_relative_count_by_flavor[flavor_index]

        if abs(current_relative_count) < 1e-12
            abs(target_relative_count) < 1e-12 || return false
            population_rescaling = 0.0
        else
            population_rescaling = target_relative_count / current_relative_count
        end

        first_state_index = projected_band_flavor_index(1, flavor_index, projected_band_count)
        last_state_index = projected_band_flavor_index(projected_band_count, flavor_index, projected_band_count)
        flavor_state_indices = first_state_index:last_state_index

        for k_index in 1:momentum_count
            density_matrix_block = @view density_matrix[flavor_state_indices, flavor_state_indices, k_index]
            density_matrix_block .*= population_rescaling

            for projected_band_index in 1:projected_band_count
                density_matrix_block[projected_band_index, projected_band_index] += reference_occupation
            end
        end
    end

    println("Using initialization seed: $seed_path")
    return true
end








function run_hartree_fock_with_oda(
    projected_basis::ProjectedBasis,
    momentum_mesh::MomentumMesh,
    single_particle_hamiltonian::Array{ComplexF64,3},
    form_factors::Array{ComplexF64,5},
    dielectric_constant::Real,
    gate_distance::Real,
    occupied_state_count_by_flavor::AbstractVector{<:Integer};
    reference_occupation::Int=0,
    maximum_iterations::Int=5000000,
    density_tolerance::Real=1e-7,
    energy_tolerance::Real=1e-8,
    verbose::Bool=true,
    seed_folder::Union{Nothing,AbstractString}=nothing,
    carrier_population_code::Int,
    displacement_sign_code::Int,
    valley_by_flavor::AbstractVector{<:Integer},
)
    projected_band_count, momentum_count, flavor_count =
        size(projected_basis.energy)

    hartree_fock_dimension =
        projected_band_count * flavor_count

    density_matrix =
        zeros(
            ComplexF64,
            hartree_fock_dimension,
            hartree_fock_dimension,
            momentum_count
        )

    filling_work =
        HartreeFockFillingWork(
            projected_band_count,
            flavor_count,
            momentum_count
        )

if rand() > 0.5 || !initialize_density_matrix_from_seed!(density_matrix, projected_basis, occupied_state_count_by_flavor, seed_folder, carrier_population_code, displacement_sign_code, valley_by_flavor; reference_occupation=reference_occupation)
    initialize_density_matrix!(density_matrix, filling_work, single_particle_hamiltonian, occupied_state_count_by_flavor)
end

    hartree_work =
        HartreeWork(
            projected_basis,
            momentum_mesh,
            dielectric_constant,
            gate_distance
        )

    fock_thread_workspace_pool =
        FockThreadWorkspacePool(projected_basis)

    hamiltonian_work =
        HartreeFockHamiltonianWork(projected_basis)

    optimal_damping_work =
        OptimalDampingWork(projected_basis)

    solution =
        solve_hartree_fock_with_oda!(
            density_matrix,
            filling_work,
            hamiltonian_work,
            optimal_damping_work,
            single_particle_hamiltonian,
            form_factors,
            hartree_work,
            fock_thread_workspace_pool,
            projected_basis,
            momentum_mesh,
            dielectric_constant,
            gate_distance,
            occupied_state_count_by_flavor;
            reference_occupation=reference_occupation,
            maximum_iterations=maximum_iterations,
            density_tolerance=density_tolerance,
            energy_tolerance=energy_tolerance,
            verbose=verbose
        )

    diagonalize_hartree_fock_hamiltonian!(
        filling_work,
        hamiltonian_work.hartree_fock_hamiltonian
    )

    final_filling_diagnostics =
        construct_zero_temperature_density_matrix!(
            filling_work,
            occupied_state_count_by_flavor
        )

    solution =
        merge(
            solution,
            (filling_diagnostics=final_filling_diagnostics,)
        )

    return (
        solution=solution,
        density_matrix=density_matrix,
        hamiltonian_work=hamiltonian_work,
        filling_work=filling_work
    )
end