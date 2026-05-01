using LinearAlgebra
struct RMGBandSpec
    name::Symbol
    spin::Int
    valley::Int
    band::Int
    cnp_shift::Float64
    z::Float64
end

struct TMDBandSpec
    name::Symbol
    spin::Int
    valley::Int
    band::Int
    mass::Float64
    energy_offset::Float64
    z::Float64
end

struct ProjectedSelection
    rmg_bands::Vector{RMGBandSpec}
    tmd_bands::Vector{TMDBandSpec}
end


struct ActiveBand
    name::Symbol
    material::Symbol
    spin::Int
    valley::Int
    band::Int
    flavor::Int
    energy_shift::Float64
    mass::Float64
    z::Float64
end




function get_f(k::Vector{Float64})
  delta1=1/√3*0.246*[0,1]
  delta2=1/√3*0.246*[√3/2,-1/2]
  delta3=1/√3*0.246*[-√3/2,-1/2]

  return exp(im*dot(k,delta1))+exp(im*dot(k,delta2))+exp(im*dot(k,delta3))
end



function Hamiltonian(k::Vector{Float64},uD::Float64,valley::Int64,stacking::Int,NL::Int)
 
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



struct ProjectedModel
    active::Vector{ActiveBand}

    k_set::Vector{Vector{Float64}}
    k_index::Vector{Vector{Int}}
    Area::Float64

    # active-band single-particle Hamiltonian
    # size: Nactive × Nactive × Nk
    single_matrix::Array{ComplexF64,3}

    # bookkeeping
    active_flavor::Vector{Int}              # active index a -> flavor id
    local_band_of_active::Vector{Int}       # active index a -> local band inside flavor
    active_in_flavor::Vector{Vector{Int}}   # flavor f -> global active indices


    # flavor f occupies this contiguous global active-index range
    flavor_ranges::Vector{UnitRange{Int}}
    # formfactor[f][n1,n2,ik,jk]
    # n1,n2 are local band indices inside flavor f
    formfactor::Vector{Array{ComplexF64,4}}
    z_flavor::Vector{Float64}
    # interaction storage
    # Vmat[f1,f2,ik,jk] = V_{f1,f2}(k_i - k_j)
    Vmat::Array{Float64,4}

    # Vzero[f1,f2] = V_{f1,f2}(q=0), used only for Hartree
    Vzero::Matrix{Float64}


    # fock_kernel[f,g][n1,n4,n3,n2] is an Nk × Nk matrix.
    # It includes V_fg(k1-k3) and projected form factors.
    # It does NOT include 1 / Area.
    fock_kernel::Matrix{Array{Matrix{ComplexF64},4}}

    # mostly for debugging
    flavor_key_to_id::Dict{Tuple{Symbol,Symbol,Int,Int},Int}
end


flavor_key(spec::RMGBandSpec) = (:RMG, spec.name, spec.spin, spec.valley)
flavor_key(spec::TMDBandSpec) = (:TMD, spec.name, spec.spin, spec.valley)

function get_flavor_id!(
    flavor_key_to_id::Dict{Tuple{Symbol,Symbol,Int,Int},Int},
    key::Tuple{Symbol,Symbol,Int,Int},
)
    if haskey(flavor_key_to_id, key)
        return flavor_key_to_id[key]
    else
        id = length(flavor_key_to_id) + 1
        flavor_key_to_id[key] = id
        return id
    end
end



function build_k_grid(radius::Float64, num_points::Int)
    kx_grid = collect(range(-radius/2, stop=radius/2, length=num_points))
    ky_grid = collect(range(-radius/2, stop=radius/2, length=num_points))

    Area = 4π^2 / ((kx_grid[2] - kx_grid[1]) * (ky_grid[2] - ky_grid[1]))

    k_set = Vector{Float64}[]
    k_index = Vector{Int}[]

    for ix in eachindex(kx_grid), iy in eachindex(ky_grid)
        push!(k_set, [kx_grid[ix], ky_grid[iy]])
        push!(k_index, [ix, iy])
    end

    return k_set, k_index, Area
end







function build_active_bands(selection::ProjectedSelection)
    active_unsorted = ActiveBand[]
    flavor_key_to_id = Dict{Tuple{Symbol,Symbol,Int,Int},Int}()

    # RMG bands
    for spec in selection.rmg_bands
        key = flavor_key(spec)
        flavor = get_flavor_id!(flavor_key_to_id, key)

            push!(active_unsorted, ActiveBand(
            spec.name,
            :RMG,
            spec.spin,
            spec.valley,
            spec.band,
            flavor,
            spec.cnp_shift,
            0.0,
            spec.z,
        ))
    end

    # TMD bands
    for spec in selection.tmd_bands
        key = flavor_key(spec)
        flavor = get_flavor_id!(flavor_key_to_id, key)

                push!(active_unsorted, ActiveBand(
            spec.name,
            :TMD,
            spec.spin,
            spec.valley,
            spec.band,
            flavor,
            spec.energy_offset,
            spec.mass,
            spec.z,
        ))
    end

    active = sort(active_unsorted, by = x -> x.flavor)

    Nactive = length(active)
    Nflavor = length(flavor_key_to_id)

    active_flavor = zeros(Int, Nactive)
    active_in_flavor = [Int[] for _ in 1:Nflavor]

    for a in 1:Nactive
        f = active[a].flavor
        active_flavor[a] = f
        push!(active_in_flavor[f], a)
    end

    local_band_of_active = zeros(Int, Nactive)

    for f in 1:Nflavor
        for (local_band, a) in enumerate(active_in_flavor[f])
            local_band_of_active[a] = local_band
        end
    end

    flavor_ranges = Vector{UnitRange{Int}}(undef, Nflavor)
    for f in 1:Nflavor
        inds = active_in_flavor[f]
        @assert !isempty(inds)
        flavor_ranges[f] = first(inds):last(inds)
    end

    return active,
           active_flavor,
           local_band_of_active,
           active_in_flavor,
           flavor_ranges,
           flavor_key_to_id
end


















function flavor_z_values(
    active::Vector{ActiveBand},
    active_in_flavor::Vector{Vector{Int}},
)
    Nflavor = length(active_in_flavor)
    z_flavor = zeros(Float64, Nflavor)

    for f in 1:Nflavor
        inds = active_in_flavor[f]
        @assert !isempty(inds)

        z0 = active[first(inds)].z
        for a in inds
            @assert active[a].z == z0 "All active bands in the same flavor must have the same z. flavor=$f"
        end

        z_flavor[f] = z0
    end

    return z_flavor
end







function build_fock_kernel(
    active_in_flavor::Vector{Vector{Int}},
    flavor_ranges::Vector{UnitRange{Int}},
    formfactor::Vector{Array{ComplexF64,4}},
    Vmat::Array{Float64,4},
    k_set::Vector{Vector{Float64}},
)
    Nflavor = length(active_in_flavor)
    Nk = length(k_set)

    fock_kernel = Matrix{Array{Matrix{ComplexF64},4}}(undef, Nflavor, Nflavor)

    for f in 1:Nflavor
        rf = flavor_ranges[f]
        nf = length(rf)
        Ff = formfactor[f]

        for g in 1:Nflavor
            rg = flavor_ranges[g]
            ng = length(rg)
            Fg = formfactor[g]

            kernels_fg = Array{Matrix{ComplexF64},4}(undef, nf, ng, nf, ng)

            Threads.@threads for n1 in 1:nf
                for n4 in 1:ng
                    for n3 in 1:nf
                        for n2 in 1:ng
                            K = zeros(ComplexF64, Nk, Nk)

                            @inbounds for ik in 1:Nk
                                for jk in 1:Nk
                                    K[ik, jk] =
                                        Vmat[f, g, ik, jk] *
                                        Ff[n1, n3, ik, jk] *
                                        Fg[n2, n4, jk, ik]
                                end
                            end

                            kernels_fg[n1, n4, n3, n2] = K
                        end
                    end
                end
            end

            fock_kernel[f, g] = kernels_fg
        end
    end

    return fock_kernel
end







function Coulomb_projected_z(
    q::Vector{Float64},
    dz::Float64,
)
    nq = norm(q)

    if nq == 0.0
        # Neutralized unscreened q=0 finite piece:
        # V_ll'(0) = -9047.5636 * |z_l - z_l'|
        return -9047.5636 * abs(dz)
    else
        # Unscreened finite-q interlayer Coulomb:
        # V_ll'(q) = 9047.5636 / |q| * exp(-|q| |z_l-z_l'|)
        return 9047.5636 / nq * exp(-abs(dz) * nq)
    end
end


function build_Vmat_z(
    z_flavor::Vector{Float64},
    k_set::Vector{Vector{Float64}};
    epsilon_r::Float64 = 1.0,
)
    Nflavor = length(z_flavor)
    Nk = length(k_set)

    Vmat = zeros(Float64, Nflavor, Nflavor, Nk, Nk)
    Vzero = zeros(Float64, Nflavor, Nflavor)

    for f1 in 1:Nflavor
        for f2 in 1:Nflavor
            dz = abs(z_flavor[f1] - z_flavor[f2])

            Vzero[f1, f2] = Coulomb_projected_z(
                [0.0, 0.0],
                dz,
            ) / epsilon_r

            Threads.@threads for ik in 1:Nk
                for jk in 1:Nk
                    Vmat[f1, f2, ik, jk] =
                        Coulomb_projected_z(
                            k_set[ik] - k_set[jk],
                            dz,
                        ) / epsilon_r
                end
            end
        end
    end

    return Vmat, Vzero
end





function build_projected_single_particle(
    radius::Float64,
    num_points::Int,
    uD::Float64,
    NL::Int,
    selection::ProjectedSelection;
    stacking::Int = 1,
)
    k_set, k_index, Area = build_k_grid(radius, num_points)

    active,
    active_flavor,
    local_band_of_active,
    active_in_flavor,
    flavor_ranges,
    flavor_key_to_id =
        build_active_bands(selection)

    Nactive = length(active)
    Nk = length(k_set)
    Nflavor = length(active_in_flavor)

    single_matrix = zeros(ComplexF64, Nactive, Nactive, Nk)

    # U_by_flavor[f] has size:
    #
    # RMG flavor:
    #   2NL × number_of_active_bands_in_flavor × Nk
    #
    # TMD flavor:
    #   number_of_active_bands_in_flavor × number_of_active_bands_in_flavor × Nk
    #
    # The TMD one is just a trivial orthonormal microscopic basis.
    U_by_flavor = Vector{Array{ComplexF64,3}}(undef, Nflavor)

    for f in 1:Nflavor
        a0 = active_in_flavor[f][1]
        nb = length(active_in_flavor[f])

        if active[a0].material == :RMG
            U_by_flavor[f] = zeros(ComplexF64, 2 * NL, nb, Nk)

        elseif active[a0].material == :TMD
            U_by_flavor[f] = zeros(ComplexF64, nb, nb, Nk)

        else
            error("Unknown material $(active[a0].material).")
        end
    end

    for a in 1:Nactive
        ab = active[a]
        f = ab.flavor
        local_band = local_band_of_active[a]

        if ab.material == :RMG
            for ik in 1:Nk
                H = Hamiltonian(k_set[ik], uD, ab.valley, stacking, NL)
                eig = eigen(Hermitian(H))

                single_matrix[a, a, ik] =
                    real(eig.values[ab.band]) + ab.energy_shift

                U_by_flavor[f][:, local_band, ik] .= eig.vectors[:, ab.band]
            end

        elseif ab.material == :TMD
            for ik in 1:Nk
                k = k_set[ik]

                single_matrix[a, a, ik] =
                    dot(k, k) / (2 * ab.mass) * 76.1996 + ab.energy_shift

                U_by_flavor[f][local_band, local_band, ik] = 1.0 + 0.0im
            end

        else
            error("Unknown active material $(ab.material).")
        end
    end

       return active,
           k_set,
           k_index,
           Area,
           single_matrix,
           active_flavor,
           local_band_of_active,
           active_in_flavor,
           flavor_ranges,
           U_by_flavor,
           flavor_key_to_id
end





function build_projected_formfactors(
    U_by_flavor::Vector{Array{ComplexF64,3}};
    atol::Float64 = 1e-12,
)
    Nflavor = length(U_by_flavor)
    formfactor = Vector{Array{ComplexF64,4}}(undef, Nflavor)

    for f in 1:Nflavor
        U = U_by_flavor[f]
        _, nb, Nk = size(U)

        Ff = zeros(ComplexF64, nb, nb, Nk, Nk)

        Threads.@threads for ik in 1:Nk
            for jk in 1:Nk

                if ik == jk
                    # Exact same-momentum form factor:
                    # <u_n(k)|u_m(k)> = δ_nm
                    for n1 in 1:nb
                        for n2 in 1:nb
                            Ff[n1, n2, ik, jk] = (n1 == n2) ? (1.0 + 0.0im) : (0.0 + 0.0im)
                        end
                    end

                else
                    # General finite-q form factor:
                    # <u_n1(k_i)|u_n2(k_j)>
                    for n1 in 1:nb
                        for n2 in 1:nb
                            val = dot(
                                view(U, :, n1, ik),
                                view(U, :, n2, jk),
                            )

                            Ff[n1, n2, ik, jk] = val
                        end
                    end
                end
            end
        end

        formfactor[f] = Ff
    end

    return formfactor
end


function build_projected_model(
    radius::Float64,
    num_points::Int,
    uD::Float64,
    NL::Int,
    selection::ProjectedSelection,
    epsilon_r::Float64;
    stacking::Int = 1
)
    active,
    k_set,
    k_index,
    Area,
    single_matrix,
    active_flavor,
    local_band_of_active,
    active_in_flavor,
    flavor_ranges,
    U_by_flavor,
    flavor_key_to_id =
        build_projected_single_particle(
            radius,
            num_points,
            uD,
            NL,
            selection;
            stacking = stacking,
        )

    formfactor = build_projected_formfactors(U_by_flavor)

    z_flavor = flavor_z_values(active, active_in_flavor)

        Vmat, Vzero = build_Vmat_z(
            z_flavor,
            k_set;
            epsilon_r = epsilon_r,
        )

    fock_kernel = build_fock_kernel(
    active_in_flavor,
    flavor_ranges,
    formfactor,
    Vmat,
    k_set,
)

        return ProjectedModel(
        active,
        k_set,
        k_index,
        Area,
        single_matrix,
        active_flavor,
        local_band_of_active,
        active_in_flavor,
        flavor_ranges,
        formfactor,
        z_flavor,
        Vmat,
        Vzero,
        fock_kernel,
        flavor_key_to_id,
    )
end




mutable struct HFWorkProjected
    Fock_matrix::Array{ComplexF64,3}          # dimension × dimension × Nk
    Hartree_matrix::Matrix{ComplexF64}        # dimension × dimension

    density_matrix_new::Array{ComplexF64,3}
    output_density_matrix::Array{ComplexF64,3}
    DeltaMatrix::Array{ComplexF64,3}

    HF_eigenvalues::Matrix{Float64}           # dimension × Nk
    HF_eigenvectors::Array{ComplexF64,3}      # dimension × dimension × Nk

    fermifactor::Vector{Float64}
    evals_flat::Vector{Float64}
    occ_matrix::Matrix{Float64}               # dimension × Nk

    H_scratch::Vector{Matrix{ComplexF64}}
    Vocc_scratch::Vector{Matrix{ComplexF64}}

    # one thread-local buffer for Ff * Pfg
   

    DIIS_input_density_matrix::Vector{Array{ComplexF64,3}}
    DIIS_input_DeltaMatrix::Vector{Array{ComplexF64,3}}
    diis_head::Int
    diis_len::Int

    flavor_density::Vector{Float64}
end



function HFWorkProjected(model::ProjectedModel, DIIS_size::Int)
    dimension = length(model.active)
    Nk = length(model.k_set)
    Nflavor = length(model.active_in_flavor)

    Z3 = zeros(ComplexF64, dimension, dimension, Nk)
    Z2 = zeros(ComplexF64, dimension, dimension)

    nscratch = Threads.maxthreadid()

    return HFWorkProjected(
        copy(Z3),  # Fock_matrix
        copy(Z2),  # Hartree_matrix

        copy(Z3),  # density_matrix_new
        copy(Z3),  # output_density_matrix
        copy(Z3),  # DeltaMatrix

        zeros(Float64, dimension, Nk),  # HF_eigenvalues
        copy(Z3),                       # HF_eigenvectors

        zeros(Float64, dimension * Nk), # fermifactor
        zeros(Float64, dimension * Nk), # evals_flat
        zeros(Float64, dimension, Nk),  # occ_matrix

        [zeros(ComplexF64, dimension, dimension) for _ in 1:nscratch], # H_scratch
        [zeros(ComplexF64, dimension, dimension) for _ in 1:nscratch], # Vocc_scratch

      

        [zeros(ComplexF64, dimension, dimension, Nk) for _ in 1:DIIS_size],
        [zeros(ComplexF64, dimension, dimension, Nk) for _ in 1:DIIS_size],
        1,
        0,

        zeros(Float64, Nflavor),
    )
end



function Construct_projector!(
    work::HFWorkProjected,
    model::ProjectedModel,
    density_matrix::Array{ComplexF64,3},
    energy_input::Float64,
    BG_density_matrix::Array{ComplexF64,3},
    target_density::Float64,
    bg_particle_density::Float64,
    temp::Float64,
    mixing::Float64,
)
    k_set = model.k_set
    Area = model.Area
    single_matrix = model.single_matrix

    dimension = length(model.active)
    Nk = length(k_set)
    Nflavor = length(model.active_in_flavor)

    @assert size(density_matrix) == (dimension, dimension, Nk)
    @assert size(BG_density_matrix) == (dimension, dimension, Nk)

    fill!(work.Fock_matrix, 0.0 + 0.0im)
    fill!(work.Hartree_matrix, 0.0 + 0.0im)
    fill!(work.flavor_density, 0.0)

    # ------------------------------------------------------------
    # Fock term
    #
    # Eq. 275:
    #
    # Fock_fg(k1)[n1,n4] =
    #   sum_{n3,n2} sum_{k3}
    #       V_fg(k1-k3)
    #       F_f[n1,n3,k1,k3]
    #       P_fg[n3,n2,k3]
    #       F_g[n2,n4,k3,k1]
    #
    # The k3 sum is done by BLAS:
    #
    #   outvec += K * densvec
    #
    # No 1/Area here. Apply it once at the end.
    # ------------------------------------------------------------

    tic = time()

    Threads.@threads for f in 1:Nflavor
        for g in 1:Nflavor
            rf = model.flavor_ranges[f]
            rg = model.flavor_ranges[g]

            nf = length(rf)
            ng = length(rg)

            kernels_fg = model.fock_kernel[f, g]

            for n1 in 1:nf
                a = rf[n1]

                for n4 in 1:ng
                    b = rg[n4]

                    outvec = @view work.Fock_matrix[a, b, :]

                    for n3 in 1:nf
                        c = rf[n3]

                        for n2 in 1:ng
                            d = rg[n2]

                            K = kernels_fg[n1, n4, n3, n2]

                            # outvec[ik] += sum_jk K[ik,jk] * density_matrix[c,d,jk]
                            mul!(
                                outvec,
                                K,
                                @view(density_matrix[c, d, :]),
                                1.0,
                                1.0,
                            )
                        end
                    end
                end
            end
        end
    end

    work.Fock_matrix .*= 1 / Area

    toc = time()
    println("Focktime", toc - tic)
    #=
    # Optional diagnostic. This does not modify Fock_matrix.
    max_fock_nonherm = 0.0
    @inbounds for ik in 1:Nk
        for a in 1:dimension
            for b in 1:dimension
                max_fock_nonherm = max(
                    max_fock_nonherm,
                    abs(work.Fock_matrix[a, b, ik] - conj(work.Fock_matrix[b, a, ik])),
                )
            end
        end
    end

    if max_fock_nonherm > 1e-8
        println("Warning: Fock_matrix is not Hermitian. max violation = ", max_fock_nonherm)
    end
    =#

    # ------------------------------------------------------------
    # Hartree term
    #
    # For now we use the simplified projected Hartree:
    #
    #   flavor_density[g] = sum_k tr(P_gg(k))
    #
    #   Hartree_f = sum_g Vzero[f,g] * flavor_density[g] / Area
    #
    # Same-momentum projected form factor is identity inside each
    # flavor block, so Hartree is diagonal in the active-band basis.
    # ------------------------------------------------------------

    tic = time()

    for g in 1:Nflavor
        rg = model.flavor_ranges[g]

        total_g = 0.0
        for ik in 1:Nk
            total_g += real(tr(@view density_matrix[rg, rg, ik]))
        end

        work.flavor_density[g] = total_g
    end

    for f in 1:Nflavor
        rf = model.flavor_ranges[f]

        hartree_f = 0.0
        for g in 1:Nflavor
            hartree_f += model.Vzero[f, g] * work.flavor_density[g]
        end

        for a in rf
            work.Hartree_matrix[a, a] = hartree_f + 0.0im
        end
    end

    work.Hartree_matrix .*= 1 / Area

    toc = time()
    println("Hartreetime", toc - tic)

    # ------------------------------------------------------------
    # Diagonalize HF Hamiltonian
    #
    # H_HF(k) = single(k) + Hartree - Fock(k)
    # ------------------------------------------------------------

    tic = time()

    old_blas_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)

    Threads.@threads for ik in 1:Nk
        tid = Threads.threadid()
        Htmp = work.H_scratch[tid]

        copy!(Htmp, work.Hartree_matrix)
        @views Htmp .-= work.Fock_matrix[:, :, ik]
        @views Htmp .+= single_matrix[:, :, ik]

        F = eigen!(Hermitian(Htmp))

        @views copy!(work.HF_eigenvectors[:, :, ik], F.vectors)
        @views copy!(work.HF_eigenvalues[:, ik], real(F.values))
    end

    BLAS.set_num_threads(old_blas_threads)

    # ------------------------------------------------------------
    # Finite-temperature occupations
    # ------------------------------------------------------------

    work.evals_flat .= vec(work.HF_eigenvalues)

    val_s = minimum(work.evals_flat)-2*temp
    val_e = maximum(work.evals_flat)+2*temp

    fermi_level, renormalized_density = find_FL_iterative!(
        work.fermifactor,
        work.evals_flat,
        target_density,
        val_s,
        val_e,
        temp,
        Area,
        bg_particle_density,
    )

    @inbounds for ik in 1:Nk, band in 1:dimension
        work.occ_matrix[band, ik] =
            fermi_occ((work.HF_eigenvalues[band, ik] - fermi_level) / temp)
    end


    
    # ------------------------------------------------------------
    # Construct new density matrix
    #
    # density_matrix_new(k) = V(k) * occ(k) * V(k)'
    # Then subtract background.
    # ------------------------------------------------------------
     
    Threads.@threads for ik in 1:Nk
        tid = Threads.threadid()
        Vocc = work.Vocc_scratch[tid]

        @views copy!(Vocc, work.HF_eigenvectors[:, :, ik])

        @inbounds for band in 1:dimension
            occ = work.occ_matrix[band, ik]
            @views Vocc[:, band] .*= occ
        end

        mul!(
            @view(work.density_matrix_new[:, :, ik]),
            Vocc,
            adjoint(@view(work.HF_eigenvectors[:, :, ik])),
            1.0,
            0.0,
        )

        @views work.density_matrix_new[:, :, ik] .-= BG_density_matrix[:, :, ik]
    end

    toc = time()
    println("constructing DMtime", toc - tic)



    # ------------------------------------------------------------
    # Mixing and diagnostics
    # ------------------------------------------------------------

    tic = time()

    @. work.output_density_matrix =
        mixing * density_matrix + (1.0 - mixing) * work.density_matrix_new

    @. work.DeltaMatrix = work.density_matrix_new - density_matrix

    eout = real(sum(abs2, work.DeltaMatrix) / Nk)

    energy = 0.0
    @inbounds for ik in 1:Nk
        energy += real(tr(
            density_matrix[:, :, ik] *
            (
                work.Hartree_matrix / 2 -
                work.Fock_matrix[:, :, ik] / 2 +
                single_matrix[:, :, ik]
            )
        )) / Nk
    end

    energy_change = real(energy - energy_input)

    toc = time()
    println("othertime", toc - tic)

    return eout,
           energy_change,
           real(energy),
           fermi_level,
           renormalized_density
end


function get_initial_proj(
    model::ProjectedModel,
    NL::Int;
    init_amp::Float64 = 0.01,
)
    dimension = length(model.active)
    Nk = length(model.k_set)

    BG_density_matrix = zeros(ComplexF64, dimension, dimension, Nk)
    initial_density_matrix = zeros(ComplexF64, dimension, dimension, Nk)

    for a in 1:dimension
        band = model.active[a]

        filled_background = false

        if band.material == :RMG
            filled_background = band.band <= NL
        elseif band.material == :TMD
            filled_background = band.mass < 0
        else
            error("Unknown material: $(band.material)")
        end

        if filled_background
            for ik in 1:Nk
                BG_density_matrix[a, a, ik] = 1.0 + 0.0im
            end
        end
    end

    for ik in 1:Nk
        A = randn(dimension, dimension) .+ im .* randn(dimension, dimension)
        initial_density_matrix[:, :, ik] .= init_amp .* (A + A')
    end

    return initial_density_matrix, BG_density_matrix
end




@inline function diis_push!(
    DIIS_input_density_matrix::Vector{Array{ComplexF64,3}},
    DIIS_input_DeltaMatrix::Vector{Array{ComplexF64,3}},
    input_density_matrix::Array{ComplexF64,3},
    DeltaMatrix::Array{ComplexF64,3},
    diis_head::Int,
    diis_len::Int,
)
    copy!(DIIS_input_density_matrix[diis_head], input_density_matrix)
    copy!(DIIS_input_DeltaMatrix[diis_head], DeltaMatrix)

    diis_head = (diis_head == length(DIIS_input_density_matrix)) ? 1 : diis_head + 1
    diis_len = min(diis_len + 1, length(DIIS_input_density_matrix))

    return diis_head, diis_len
end


function iteration(
    initial_density_matrix::Array{ComplexF64,3},
    BG_density_matrix::Array{ComplexF64,3},
    model::ProjectedModel,
    target_density::Float64,
    temp::Float64;
    DIIS_size::Int = 5,
    mixing::Float64 = 0.5
)
    dimension = length(model.active)
    Nk = length(model.k_set)

    @assert size(initial_density_matrix) == (dimension, dimension, Nk)
    @assert size(BG_density_matrix) == (dimension, dimension, Nk)

    bg_particle_density = background_density(BG_density_matrix, model.Area)

    eout = 1.0
    itcount = 0
    bad_count = 0
    energy = 0.0
    energy_change = 0.0
    fermi_level = 0.0
    renormalized_density = 0.0

    input_density_matrix = copy(initial_density_matrix)

    work = HFWorkProjected(model, DIIS_size)

    eout_hist = Float64[]
    PLATEAU_N = 10
    PLATEAU_FRAC = 0.10
    E_EPS = 1e-30

    diis_fire_once = false
    diis_cooldown = 0
    DIIS_COOLDOWN = 10

    while ((eout > 1e-14) || (bad_count < DIIS_size + 2) || (abs(energy_change) > 1e-8))

        if eout < 1e-14
            bad_count += 1
        else
            bad_count = 0
        end


        

        dmk_used = input_density_matrix

        tic = time()

        use_diis =
            ((itcount > 60 && abs(eout) > 1e-2) ||
             (itcount > 50 && abs(eout) < 1e-6) ||
             diis_fire_once) &&
            (work.diis_len >= 4)

        if use_diis
            dmk = implement_DIIS(
                work.DIIS_input_density_matrix,
                work.DIIS_input_DeltaMatrix,
                model.k_set,
                work.diis_len,
            )

            if dmk === nothing
                A = randn(dimension, dimension, Nk) .+ im .* randn(dimension, dimension, Nk)
                dmk = similar(input_density_matrix)

                @inbounds for ik in 1:Nk
                    dmk[:, :, ik] .= 0.01 .* (A[:, :, ik] + A[:, :, ik]')
                end

                itcount = 0
                work.diis_head = 1
                work.diis_len = 0
                empty!(eout_hist)
                diis_fire_once = false
                diis_cooldown = 0

                println("DIIS failed, random restart")
            end

            dmk_used = dmk

            eout, energy_change, energy, fermi_level, renormalized_density =
                Construct_projector!(
                    work,
                    model,
                    dmk,
                    energy,
                    BG_density_matrix,
                    target_density,
                    bg_particle_density,
                    temp,
                    mixing,
                )

            println("using DIIS")

            if diis_fire_once
                diis_fire_once = false
                empty!(eout_hist)
                diis_cooldown = DIIS_COOLDOWN
            end
        else
            eout, energy_change, energy, fermi_level, renormalized_density =
                Construct_projector!(
                    work,
                    model,
                    input_density_matrix,
                    energy,
                    BG_density_matrix,
                    target_density,
                    bg_particle_density,
                    temp,
                    mixing,
                )
        end

        work.diis_head, work.diis_len = diis_push!(
            work.DIIS_input_density_matrix,
            work.DIIS_input_DeltaMatrix,
            dmk_used,
            work.DeltaMatrix,
            work.diis_head,
            work.diis_len,
        )

        copy!(input_density_matrix, work.output_density_matrix)

        itcount += 1
     



        println(
            time() - tic,
            " eout=$eout",
            " energy_change=$energy_change",
            " energy=$energy",
            " fermi_level=$fermi_level",
            " density=$renormalized_density",
            " itcount=$itcount",
        )
        flush(stdout)

   

        if diis_cooldown > 0
            diis_cooldown -= 1
        end

        push!(eout_hist, eout)
        if length(eout_hist) > PLATEAU_N
            popfirst!(eout_hist)
        end

        if !diis_fire_once && diis_cooldown == 0 && length(eout_hist) == PLATEAU_N
            e0 = eout_hist[1]
            e1 = eout_hist[end]

            rel_change = abs(e1 - e0) / max(abs(e0), E_EPS)

            if rel_change < PLATEAU_FRAC
                diis_fire_once = true
                println("Plateau detected: |Δe|/|e| ≈ $(rel_change). Will fire DIIS once.")
            end
        end
    end



    return work.HF_eigenvalues,
           work.HF_eigenvectors,
           energy,
           work.DIIS_input_density_matrix,
           fermi_level,
           work.Hartree_matrix,
           work.Fock_matrix,
           eout,
           renormalized_density
end


function implement_DIIS(
    DIIS_input_projector::Vector{Array{ComplexF64,3}},
    DIIS_input_DeltaMatrix::Vector{Array{ComplexF64,3}},
    k_set::Vector{Vector{Float64}},
    DIIS_size::Int,
)
    Bmatrix = zeros(ComplexF64, DIIS_size + 1, DIIS_size + 1)

    for ja in 1:DIIS_size
        Bmatrix[ja, DIIS_size + 1] = 1
        Bmatrix[DIIS_size + 1, ja] = 1
    end

    for ja in 1:DIIS_size
        for jb in 1:DIIS_size
            for jk in eachindex(k_set)
                Bmatrix[ja, jb] += real(tr(
                    DIIS_input_DeltaMatrix[ja][:, :, jk]' *
                    DIIS_input_DeltaMatrix[jb][:, :, jk]
                ))
            end
        end
    end

    inB = safe_inverse(Bmatrix)

    if inB !== nothing
        rhs = zeros(Float64, DIIS_size + 1)
        rhs[end] = 1.0

        coeff = inB * rhs

        dmk = zero(DIIS_input_projector[1])

        @inbounds for ja in 1:DIIS_size
            BLAS.axpy!(coeff[ja], DIIS_input_projector[ja], dmk)
            BLAS.axpy!(coeff[ja], DIIS_input_DeltaMatrix[ja], dmk)
        end

        return dmk
    else
        return nothing
    end
end

function safe_inverse(A)
    try
        return inv(A)
    catch e
        if isa(e, SingularException)
            println("Matrix is singular, doing pseudoinverse.")
            return pinv(A, 1e-8)
        else
            return nothing
        end
    end
end



function background_density(BG_density_matrix::Array{ComplexF64,3}, Area::Float64)
    Nk = size(BG_density_matrix, 3)

    total = 0.0
    for ik in 1:Nk
        total += real(tr(@view BG_density_matrix[:, :, ik]))
    end

    return total / Area
end





@inline function fermi_occ(x::Float64)::Float64
    if x > 40.0
        return 0.0
    elseif x < -40.0
        return 1.0
    else
        return 1.0 / (exp(x) + 1.0)
    end
end


function find_FL_iterative!(
    fermifactor::Vector{Float64},
    quasi_particle_energy::Vector{Float64},
    target_density::Float64,
    val_s::Float64,
    val_e::Float64,
    temp::Float64,
    Area::Float64,
    bg_particle_density::Float64;
    maxiter::Int = 300,
)
    lo = val_s
    hi = val_e
    try_FL = 0.5 * (lo + hi)

    stan = target_density == 0.0 ? 1e-9 : abs(1e-8 * target_density)

    fl = 0.0
    for _ in 1:maxiter
        try_FL = 0.5 * (lo + hi)

        @inbounds for i in eachindex(quasi_particle_energy)
            fermifactor[i] = fermi_occ((quasi_particle_energy[i] - try_FL) / temp)
        end

        fl = sum(fermifactor) / Area - bg_particle_density

        if abs(fl - target_density) < stan
            return try_FL, fl
        elseif fl > target_density
            hi = try_FL
        else
            lo = try_FL
        end
    end

    return try_FL, fl
end





function get_projected_selection(
    sel::Int,
    NL::Int,
    deltaE::Float64,
    m_TMD::Float64,
    z_TMD::Float64;
)
    z_RMG = 0.0

    if sel == 1
        return ProjectedSelection(
            [
                RMGBandSpec(:RMG_s1_vp, 1, +1, NL+1, 0.0, z_RMG),
                RMGBandSpec(:RMG_s1_vm, 1, -1, NL+1, 0.0, z_RMG),
                RMGBandSpec(:RMG_s2_vp, 2, +1, NL+1, 0.0, z_RMG),
                RMGBandSpec(:RMG_s2_vm, 2, -1, NL+1, 0.0, z_RMG),
            ],
            [
                TMDBandSpec(:TMD_s1, 1, 0, 1, m_TMD, deltaE, z_TMD),
            ],
        )
    else
        error("For now only sel == 1 is implemented.")
    end
end