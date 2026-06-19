


using LinearAlgebra


struct StripeSPBookkeeping
    # single-particle band data
    band_energies::Vector{Float64}
    band_vecs::Matrix{ComplexF64}

    # stripe single-particle Hamiltonian
    single_particle_matrix::Array{ComplexF64,3}  # (Nn, Nn, Nkred)

    # full momentum grid
    kvecs_cart::Vector{Vector{Float64}}
    k_int::Vector{NTuple{2,Int}}
    kint_to_gid::Dict{NTuple{2,Int},Int}

    # stripe ordering data
    Qint::NTuple{2,Int}
    Qvec_cart::Vector{Float64}

    k_equiv_id::Vector{Vector{Int}}
    k_equiv_n::Vector{Vector{Int}}
    k_repre_id::Vector{Int}

    nvals::Vector{Int}
    gid_padded::Matrix{Int}
    valid_ns::Vector{Vector{Int}}

    # interaction/form factor data
    formfactors::Array{ComplexF64,4}
    VF_matrix::Array{ComplexF64,4}

    # geometry
    Area::Float64
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






function get_single_particle(
    radius::Float64,
    num_points::Int,
    uD::Float64,
    NL::Int,
    grid_angle::Float64,
    Qvec_index::Vector{Int};
    band_index::Int,
    stacking::Int = 1,
)

  valley_index=1
  if iseven(num_points)
      error("Use odd num_points so that the grid contains k=0.")
  end


  kx_grid = collect(range(-radius / 2, stop = radius / 2, length = num_points))
  ky_grid = collect(range(-radius / 2, stop = radius / 2, length = num_points))

  dkx = kx_grid[2] - kx_grid[1]
  dky = ky_grid[2] - ky_grid[1]

  Rθ = [
      cos(grid_angle/180*π) -sin(grid_angle/180*π)
      sin(grid_angle/180*π)  cos(grid_angle/180*π)
  ]

  T1 = Rθ * [dkx, 0.0]
  T2 = Rθ * [0.0, dky]

  Area = 4π^2 / abs(T1[1] * T2[2] - T1[2] * T2[1])

  Qint = (Qvec_index[1], Qvec_index[2])
  Qvec_cart = Qint[1] * T1 + Qint[2] * T2

   if Qint==(0,0)
     error("this is wrong Q")
   end

    norb = 2 * NL
    Nk = num_points^2

    if band_index < 1 || band_index > norb
        error("band_index must be between 1 and $(norb).")
    end

    band_energies = Vector{Float64}(undef, Nk)
    band_vecs = Matrix{ComplexF64}(undef, norb, Nk)

    kvecs_cart = Vector{Vector{Float64}}(undef, Nk)
    k_int = Vector{NTuple{2,Int}}(undef, Nk)



    gid = 0

    for ix in eachindex(kx_grid), iy in eachindex(ky_grid)
        gid += 1

        nx = ix - (num_points + 1) ÷ 2
        ny = iy - (num_points + 1) ÷ 2

        k = nx * T1 + ny * T2
        Ham = Hamiltonian(k, uD, valley_index, stacking, NL)

        FFF = eigen(Hermitian(Ham))

        band_energies[gid] = real(FFF.values[band_index])
        uvec=FFF.vectors[:, band_index]
        if (uD>=0 && band_index==NL+1)|| (uD<0 && band_index==NL)
             agg=angle(uvec[2*NL])
             uvec=uvec*exp(-im*agg)
        elseif  (uD>=0 && band_index==NL)|| (uD<0 && band_index==NL+1)
             agg=angle(uvec[1])
             uvec=uvec*exp(-im*agg)
        else
            error("not implemented gauge")

        end


        band_vecs[:, gid] .= uvec

        kvecs_cart[gid] = k
       k_int[gid] = (nx, ny)
    end

    kint_to_gid = Dict(k_int[i] => i for i in eachindex(k_int))

    # Large enough search range along Q.
    ncut = Int(ceil(radius * sqrt(2) / norm(Qvec_cart))) + 5

    have_visited = falses(length(k_int))

    # k_equiv_id[class_id] = list of global momentum gids in one equivalence class
    k_equiv_id = Vector{Vector{Int}}()

    # k_equiv_m[class_id] = the corresponding raw integer m relative to the seed point
    k_equiv_m = Vector{Vector{Int}}()

    for ja in eachindex(k_int)

        if have_visited[ja] == false

            temp_equiv = Int[]
            temp_m = Int[]

            seed = k_int[ja]

            for mm in -ncut:ncut
                kkk = (seed[1] + mm * Qint[1], seed[2] + mm * Qint[2])

                if haskey(kint_to_gid, kkk)
                    gid_now = kint_to_gid[kkk]

                    push!(temp_equiv, gid_now)
                    push!(temp_m, mm)
                    have_visited[gid_now] = true
                end
            end

            # Sort along Q by raw m.
            perm = sortperm(temp_m)
            temp_equiv = temp_equiv[perm]
            temp_m = temp_m[perm]

            push!(k_equiv_id, temp_equiv)
            push!(k_equiv_m, temp_m)
        end
    end

    # Choose representative.
    #
    # For memory-efficient padded 4D form factors, I recommend choosing
    # the first point along the finite Q-chain, so n = 0,1,...,L-1.
    # This avoids having different classes like [-10,2] and [3,10].
    k_repre_id = zeros(Int, length(k_equiv_id))
    k_equiv_n = [zeros(Int, length(k_equiv_id[ja])) for ja in eachindex(k_equiv_id)]

    for ja in eachindex(k_equiv_id)
        # representative = first point along Q
        k_repre_id[ja] = k_equiv_id[ja][1]

        m0 = k_equiv_m[ja][1]

        for jb in eachindex(k_equiv_id[ja])
            # Now n is always 0,1,2,... inside each stripe.
            k_equiv_n[ja][jb] = k_equiv_m[ja][jb] - m0
        end
    end

    max_n = maximum(maximum.(k_equiv_n))
    min_n = minimum(minimum.(k_equiv_n))
    @assert min_n == 0

    nvals = collect(min_n:max_n)
    Nn = length(nvals)
    Nkred = length(k_equiv_id)

    # gid_padded[n_index, k_reduced_index] gives global momentum gid.
    # 0 means invalid.
    gid_padded = zeros(Int, Nn, Nkred)
    valid_ns = [Int[] for _ in 1:Nkred]

    for kred in 1:Nkred
        for jj in eachindex(k_equiv_id[kred])
            gid_now = k_equiv_id[kred][jj]
            n_now = k_equiv_n[kred][jj]

            n_index = n_now - min_n + 1

            gid_padded[n_index, kred] = gid_now
            push!(valid_ns[kred], n_index)
        end
    end
    

    single_particle_matrix = zeros(ComplexF64, Nn, Nn, Nkred)

    for k_red in 1:Nkred
        for n_index in valid_ns[k_red]
            gid = gid_padded[n_index, k_red]
            single_particle_matrix[n_index, n_index, k_red] = band_energies[gid]
        end
    end


   

    formfactors = zeros(ComplexF64, Nn, Nn, Nkred, Nkred)

    @inbounds for k1_red in 1:Nkred
        for k2_red in 1:Nkred
            for n1_index in valid_ns[k1_red]
                gid1 = gid_padded[n1_index, k1_red]

                for n2_index in valid_ns[k2_red]
                    gid2 = gid_padded[n2_index, k2_red]

                    formfactors[n1_index, n2_index, k1_red, k2_red] =
                        dot(band_vecs[:, gid1], band_vecs[:, gid2])
                end
            end
        end
    end




    VF_matrix = zeros(ComplexF64, Nn, Nn, Nkred, Nkred)

    @inbounds for k1_red in 1:Nkred
        for k2_red in 1:Nkred
            for n1_index in valid_ns[k1_red]
                gid1 = gid_padded[n1_index, k1_red]

                for n3_index in valid_ns[k2_red]
                    gid3 = gid_padded[n3_index, k2_red]

                    dq1 = k_int[gid1][1] - k_int[gid3][1]
                    dq2 = k_int[gid1][2] - k_int[gid3][2]

                    VF_matrix[n1_index, n3_index, k1_red, k2_red] =
                        Coulomb([dq1, dq2], T1, T2) *
                        formfactors[n1_index, n3_index, k1_red, k2_red]
                end
            end
        end
    end




   return StripeSPBookkeeping(
    band_energies,
    band_vecs,

    single_particle_matrix,

    kvecs_cart,
    k_int,
    kint_to_gid,

    Qint,
    Qvec_cart,

    k_equiv_id,
    k_equiv_n,
    k_repre_id,

    nvals,
    gid_padded,
    valid_ns,

    formfactors,
    VF_matrix,

    Area,
)


       

end


function Coulomb(k::Vector{Int},T1::Vector{Float64},T2::Vector{Float64})::Float64
   D=25
   return k==[0,0] ? D*9047.5636 : tanh(norm([T1 T2]*k*D))/norm(k[1]*T1+k[2]*T2)*9047.5636
end


function construct_n_index_pairs(nvals::Vector{Int})

    Nn = length(nvals)
    valid_n_pairs = NTuple{4,Int}[]

    for n1 in 1:Nn
        for n2 in 1:Nn
            for n3 in n2:Nn

                n4 = n1 + n2 - n3

                if 1 <= n4 <= Nn
                    push!(valid_n_pairs, (n1, n2, n3, n4))
                end
            end
        end
    end

    return valid_n_pairs
end



mutable struct StripeHFWork
    Fock_matrix::Array{ComplexF64,3}          # (Nn, Nn, Nkred)
    Hartree_matrix::Array{ComplexF64,3}       # (Nn, Nn, Nkred)
    HF_matrix::Array{ComplexF64,3}            # (Nn, Nn, Nkred)

    density_matrix_new::Array{ComplexF64,3}
    output_density_matrix::Array{ComplexF64,3}
    DeltaMatrix::Array{ComplexF64,3}

    HF_eigenvalues::Matrix{Float64}           # (Nn, Nkred)
    HF_eigenvectors::Array{ComplexF64,3}      # (Nn, Nn, Nkred)

    fermifactor::Vector{Float64}
    evals_flat::Vector{Float64}
    occ_matrix::Matrix{Float64}               # (Nn, Nkred)

    H_scratch::Vector{Matrix{ComplexF64}}
    Vocc_scratch::Vector{Matrix{ComplexF64}}
end


function StripeHFWork(book)

    Nn, Nkred = size(book.gid_padded)

    Z3 = zeros(ComplexF64, Nn, Nn, Nkred)
    nscratch = Threads.maxthreadid()

    return StripeHFWork(
        copy(Z3),                         # Fock_matrix
        copy(Z3),                         # Hartree_matrix
        copy(Z3),                         # HF_matrix

        copy(Z3),                         # density_matrix_new
        copy(Z3),                         # output_density_matrix
        copy(Z3),                         # DeltaMatrix

        zeros(Float64, Nn, Nkred),        # HF_eigenvalues
        copy(Z3),                         # HF_eigenvectors

        zeros(Float64, Nn * Nkred),       # fermifactor
        zeros(Float64, Nn * Nkred),       # evals_flat
        zeros(Float64, Nn, Nkred),        # occ_matrix

        [zeros(ComplexF64, Nn, Nn) for _ in 1:nscratch],
        [zeros(ComplexF64, Nn, Nn) for _ in 1:nscratch],
    )
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
    quasi_particle_energy::AbstractVector{Float64},
    target_density::Float64,
    val_s::Float64,
    val_e::Float64,
    temp::Float64,
    Area::Float64;
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

        fl = sum(fermi_occ((quasi_particle_energy[i] - try_FL) / temp)
                 for i in eachindex(quasi_particle_energy)) / Area

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




function construct_projector!(
    work::StripeHFWork,
    book::StripeSPBookkeeping,
    valid_n_pairs::Vector{NTuple{4,Int}},
    density_matrix::Array{ComplexF64,3},
    energy_input::Float64,
    target_density::Float64,
    temp::Float64,
    mixing::Float64,
    ϵr::Float64,
    band_index::Int,
    NL::Int,
)

    Nn, Nkred = size(book.gid_padded)

    fill!(work.Fock_matrix, 0)
    fill!(work.Hartree_matrix, 0)
    fill!(work.HF_matrix, 0)

    fill!(work.density_matrix_new, 0)
    fill!(work.output_density_matrix, 0)
    fill!(work.DeltaMatrix, 0)

    fill!(work.HF_eigenvalues, Inf)
    fill!(work.HF_eigenvectors, 0)
    fill!(work.occ_matrix, 0)

    # ============================================================
    # Fock construction goes here.
    #
    # Fill:
    #     work.Fock_matrix
    #
    # Available:
    #     book.formfactors
    #     book.VF_matrix
    #     book.gid_padded
    #     book.valid_ns
    #     valid_n_pairs
    #     density_matrix
    #     ϵr
    #     book.Area
    # ============================================================


       Threads.@threads for k2_red in 1:Nkred

        dim2 = length(book.valid_ns[k2_red])

        for (n1, n2, n3, n4) in valid_n_pairs

            if n3 <= dim2

                temp_sum = sum(
                    @views book.VF_matrix[n1, n3, :, k2_red] .*
                           book.formfactors[n2, n4, k2_red, :] .*
                           density_matrix[n4, n1, :]
                )

                work.Fock_matrix[n2, n3, k2_red] += temp_sum 

            end
        end

        for n2 in 1:dim2
            work.Fock_matrix[n2, n2, k2_red] =
                complex(real(work.Fock_matrix[n2, n2, k2_red]), 0.0)

            for n3 in (n2 + 1):dim2
                work.Fock_matrix[n3, n2, k2_red] =
                    conj(work.Fock_matrix[n2, n3, k2_red])
            end
        end
    end

    work.Fock_matrix.*=(2.0 / (ϵr * book.Area))



        # ============================================================
    # Hartree construction.
    #
    # First construct:
    #
    #   hartree_sum[n1,n4] =
    #       sum_k1 VF(n1,k1; n4,k1) * P(k1)[n4,n1]
    #
    # Then:
    #
    #   Hartree(k2)[n2,n3] =
    #       2/A * F(n2,k2; n3,k2)
    #       * sum_{n1,n4} hartree_sum[n1,n4]
    #       subject to n1+n2=n3+n4.
    #
    # valid_n_pairs stores (n1,n2,n3,n4), with n2 <= n3.
    # ============================================================

        # Hartree construction

    hartree_sum = zeros(ComplexF64, Nn, Nn)

    for n1 in 1:Nn, n4 in 1:Nn
        hartree_sum[n1, n4] = sum(
            book.VF_matrix[n1, n4, k1_red, k1_red] *
            density_matrix[n4, n1, k1_red]
            for k1_red in 1:Nkred
        )
    end

    Threads.@threads for k2_red in 1:Nkred
        dim2 = length(book.valid_ns[k2_red])

        for (n1, n2, n3, n4) in valid_n_pairs
            if n3 <= dim2
                work.Hartree_matrix[n2, n3, k2_red] +=
                    book.formfactors[n2, n3, k2_red, k2_red] *
                    hartree_sum[n1, n4] 
            end
        end

        for n2 in 1:dim2
            work.Hartree_matrix[n2, n2, k2_red] =
                real(work.Hartree_matrix[n2, n2, k2_red])

            for n3 in n2+1:dim2
                work.Hartree_matrix[n3, n2, k2_red] =
                    conj(work.Hartree_matrix[n2, n3, k2_red])
            end
        end
    end


    work.Hartree_matrix.*=(2.0 / (ϵr * book.Area))



    old_blas_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)

    Threads.@threads for k_red in 1:Nkred
        tid = Threads.threadid()
        Htmp = work.H_scratch[tid]

        dimk = length(book.valid_ns[k_red])
        r = 1:dimk

        fill!(Htmp, 0)

        @views Htmp[r, r] .=
            book.single_particle_matrix[r, r, k_red] .+
            work.Hartree_matrix[r, r, k_red] .-
            work.Fock_matrix[r, r, k_red]

        FFF = eigen!(Hermitian(@view Htmp[r, r]))

        @views work.HF_eigenvalues[r, k_red] .= real.(FFF.values)
        @views work.HF_eigenvectors[r, r, k_red] .= FFF.vectors
        @views work.HF_matrix[r, r, k_red] .= Htmp[r, r]
    end

    BLAS.set_num_threads(old_blas_threads)

    # collect only physical eigenvalues
    eval_count = 0
    for k_red in 1:Nkred
        dimk = length(book.valid_ns[k_red])
        r = 1:dimk

        @views work.evals_flat[eval_count + 1 : eval_count + dimk] .=
            work.HF_eigenvalues[r, k_red]

        eval_count += dimk
    end

    evals_view = @view work.evals_flat[1:eval_count]

    val_s = minimum(evals_view)
    val_e = maximum(evals_view)
     




        full_density = eval_count / book.Area

        if band_index < NL + 1
            target_electron_density = full_density + target_density
        else
            target_electron_density = target_density
        end

        if target_electron_density < -1e-12 || target_electron_density > full_density + 1e-12
            error(
                "target_density inconsistent with this projected band. " *
                "target_density=$target_density, " *
                "target_electron_density=$target_electron_density, " *
                "full_density=$full_density"
            )
        end

        fermi_level, renormalized_electron_density = find_FL_iterative!(
            work.fermifactor,
            evals_view,
            target_electron_density,
            val_s,
            val_e,
            temp,
            book.Area,
        )

        if band_index < NL + 1
            renormalized_density = renormalized_electron_density - full_density
        else
            renormalized_density = renormalized_electron_density
        end
    for k_red in 1:Nkred
        dimk = length(book.valid_ns[k_red])

        for band in 1:dimk
            work.occ_matrix[band, k_red] =
                fermi_occ((work.HF_eigenvalues[band, k_red] - fermi_level) / temp)
        end
    end

    Threads.@threads for k_red in 1:Nkred
        tid = Threads.threadid()
        Vocc = work.Vocc_scratch[tid]

        dimk = length(book.valid_ns[k_red])
        r = 1:dimk

        fill!(Vocc, 0)

        @views Vocc[r, r] .= work.HF_eigenvectors[r, r, k_red]

        for band in 1:dimk
            @views Vocc[r, band] .*= work.occ_matrix[band, k_red]
        end

        mul!(
            @view(work.density_matrix_new[r, r, k_red]),
            @view(Vocc[r, r]),
            adjoint(@view(work.HF_eigenvectors[r, r, k_red])),
            1.0,
            0.0,
        )
    end

    @. work.output_density_matrix =
        mixing * density_matrix + (1.0 - mixing) * work.density_matrix_new

    @. work.DeltaMatrix =
        work.density_matrix_new - density_matrix

    eout = real(sum(abs2, work.DeltaMatrix) / eval_count)

    energy = 0.0

    for k_red in 1:Nkred
        dimk = length(book.valid_ns[k_red])
        r = 1:dimk

        @views energy_matrix =
            book.single_particle_matrix[r, r, k_red] .+
            0.5 .* work.Hartree_matrix[r, r, k_red] .-
            0.5 .* work.Fock_matrix[r, r, k_red]

        @views energy += real(tr(density_matrix[r, r, k_red] * energy_matrix))
    end

    energy = energy / eval_count
    energy_change = real(energy - energy_input)

    return eout, energy_change, real(energy), fermi_level, renormalized_density
end


@inline function diis_push_stripe!(
    DIIS_input_density_matrix::Vector{Array{ComplexF64,3}},
    DIIS_input_DeltaMatrix::Vector{Array{ComplexF64,3}},
    input_density_matrix::Array{ComplexF64,3},
    DeltaMatrix::Array{ComplexF64,3},
    diis_head::Int,
    diis_len::Int,
)

    copy!(DIIS_input_density_matrix[diis_head], input_density_matrix)
    copy!(DIIS_input_DeltaMatrix[diis_head], DeltaMatrix)

    if diis_head == length(DIIS_input_density_matrix)
        diis_head = 1
    else
        diis_head += 1
    end

    diis_len = min(diis_len + 1, length(DIIS_input_density_matrix))

    return diis_head, diis_len
end


function safe_inverse(A)
    try
        return inv(A)
    catch e
        if isa(e, SingularException)
            println("Matrix is singular, doing pseudoinverse.")
            return pinv(A, 1e-8)
        else
            rethrow(e)
        end
    end
end


function implement_DIIS_stripe(
    DIIS_input_density_matrix::Vector{Array{ComplexF64,3}},
    DIIS_input_DeltaMatrix::Vector{Array{ComplexF64,3}},
    diis_len::Int,
)

    Bmatrix = zeros(ComplexF64, diis_len + 1, diis_len + 1)

    for i in 1:diis_len
        Bmatrix[i, diis_len + 1] = 1
        Bmatrix[diis_len + 1, i] = 1
    end

    for i in 1:diis_len
        for j in 1:diis_len
            Bmatrix[i, j] = real(
                sum(conj.(DIIS_input_DeltaMatrix[i]) .* DIIS_input_DeltaMatrix[j])
            )
        end
    end

    inB = safe_inverse(Bmatrix)

    rhs = zeros(Float64, diis_len + 1)
    rhs[end] = 1.0

    coeff = inB * rhs

    dmk = zero(DIIS_input_density_matrix[1])

    for i in 1:diis_len
        BLAS.axpy!(coeff[i], DIIS_input_density_matrix[i], dmk)
        BLAS.axpy!(coeff[i], DIIS_input_DeltaMatrix[i], dmk)
    end

    return dmk
end




function iteration_stripe(
    initial_density_matrix::Array{ComplexF64,3},
    book::StripeSPBookkeeping,
    ϵr::Float64,
    target_density::Float64,
    temp::Float64,
    bandindex::Int,
    NL::Int;
    DIIS_size::Int = 5,
    mixing::Float64 = 0.5,
    eout_tol::Float64 = 1e-12,
    energy_tol::Float64 = 1e-8,
)

    valid_n_pairs = construct_n_index_pairs(book.nvals)

    work = StripeHFWork(book)

    input_density_matrix = copy(initial_density_matrix)

    DIIS_input_density_matrix =
        [similar(input_density_matrix) for _ in 1:DIIS_size]

    DIIS_input_DeltaMatrix =
        [similar(input_density_matrix) for _ in 1:DIIS_size]

    diis_head = 1
    diis_len = 0

    eout = 1.0
    energy = 0.0
    energy_change = 0.0
    fermi_level = 0.0
    renormalized_density = 0.0

    itcount = 0
    bad_count = 0

    eout_hist = Float64[]
    PLATEAU_N = 10
    PLATEAU_FRAC = 0.10
    E_EPS = 1e-30

    diis_fire_once = false
    diis_cooldown = 0
    DIIS_COOLDOWN = 10

    while true

        if eout < eout_tol
            bad_count += 1
        else
            bad_count = 0
        end

        if bad_count >= DIIS_size + 2 && abs(energy_change) < energy_tol
            break
        end
        dmk_used = input_density_matrix

        use_diis =
            (
                (itcount > 60 && abs(eout) > 1e-2) ||
                (itcount > 50 && abs(eout) < 1e-6) ||
                diis_fire_once
            ) && (diis_len >= 4)

        tic = time()

        if use_diis
            dmk = implement_DIIS_stripe(
                DIIS_input_density_matrix,
                DIIS_input_DeltaMatrix,
                diis_len,
            )

            dmk_used = dmk

            eout, energy_change, energy, fermi_level, renormalized_density =
                construct_projector!(
                    work,
                    book,
                    valid_n_pairs,
                    dmk,
                    energy,
                    target_density,
                    temp,
                    mixing,
                    ϵr,
                    bandindex,
                    NL
                )

            println("using DIIS")

            if diis_fire_once
                diis_fire_once = false
                empty!(eout_hist)
                diis_cooldown = DIIS_COOLDOWN
            end
        else
            eout, energy_change, energy, fermi_level, renormalized_density =
                construct_projector!(
                    work,
                    book,
                    valid_n_pairs,
                    input_density_matrix,
                    energy,
                    target_density,
                    temp,
                    mixing,
                    ϵr,
                    bandindex,
                    NL
                )
        end

        diis_head, diis_len = diis_push_stripe!(
            DIIS_input_density_matrix,
            DIIS_input_DeltaMatrix,
            dmk_used,
            work.DeltaMatrix,
            diis_head,
            diis_len,
        )

        copy!(input_density_matrix, work.output_density_matrix)

        itcount += 1

        println(
            "time=", time() - tic,
            " eout=", eout,
            " energy=", energy,
            " energy_change=", energy_change,
            " fermi_level=", fermi_level,
            " density=", renormalized_density,
            " itcount=", itcount,
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
           input_density_matrix,
           fermi_level,
           work.Hartree_matrix,
           work.Fock_matrix,
           eout,
           renormalized_density
end



function random_hermitian_density_matrix(book; random_amp::Float64 = 0.1)


    Nn, Nkred = size(book.gid_padded)
    density_matrix = zeros(ComplexF64, Nn, Nn, Nkred)

    for k_red in 1:Nkred
        dimk = length(book.valid_ns[k_red])
        r = 1:dimk

        A = randn(ComplexF64, dimk, dimk)
        H = random_amp .* (A + A') ./ 2

        @views density_matrix[r, r, k_red] .= H
    end

    return density_matrix
end
