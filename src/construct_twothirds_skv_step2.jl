using EllipticFunctions
using LinearAlgebra
using Plots
using JLD2,StaticArrays
using Combinatorics
BLAS.set_num_threads(1)




function get_spinor_norm(NL::Int,q::Vector{Float64})
    qh=q[1]+im*q[2]
    aGr=0.246
    t0=3100
    t1=380
    vf=√3/2*aGr*t0
    spinor=zeros(ComplexF64,NL)
    for ja in 1:NL
      spinor[ja]=qh^(ja-1)*(vf/t1)^(ja-1)
    end
  return 1/norm(spinor)
end


function get_spinor(NL::Int,q::Vector{Float64})
    qh=q[1]+im*q[2]
    aGr=0.246
    t0=3100
    t1=380
    vf=√3/2*aGr*t0
    spinor=zeros(ComplexF64,NL)
    for ja in 1:NL
      spinor[ja]=qh^(ja-1)*(vf/t1)^(ja-1)
    end
  return spinor/norm(spinor)
end


struct sigma_related
    pv::ComplexF64
    ev::ComplexF64
end


# ------------------------------------------------------------
# Slaterorb = occupied-orbital pair (p,q), with p < q
# basisorb  = basis-orbital pair (a,b)
# ------------------------------------------------------------


# C[I,J] = (-1)^(p+q+k+l) det S^(pq|kl)
#
# Here:
#   Slaterorb_pairs[I] = (p,q)
#   Slaterorb_pairs[J] = (k,l)
#
# This is the singular-safe second-cofactor matrix.

struct QResolvedPrecomp
    Tn1::Vector{Int}
    Tn2::Vector{Int}
    basisorb_total_groups::Vector{Vector{Tuple{Int,Int}}}
    Ω::Matrix{ComplexF64}
end

function QResolvedPrecomp(Tgrid, spinor_set)
    Tn1, Tn2 = Tgrid_components(Tgrid)
    basisorb_total_groups = make_basisorb_total_groups(Tn1, Tn2)
    Ω = spinor_overlap_matrix(spinor_set)

    return QResolvedPrecomp(Tn1, Tn2, basisorb_total_groups, Ω)
end



function Slaterorb_second_minor_matrix(S::AbstractMatrix)
    N = size(S, 1)

    Slaterorb_pairs = [(p, q) for p in 1:N-1 for q in p+1:N]
    nSlaterorb = length(Slaterorb_pairs)

    Slaterorb_keep = [
        [r for r in 1:N if r != p && r != q]
        for (p, q) in Slaterorb_pairs
    ]

    C = zeros(ComplexF64, nSlaterorb, nSlaterorb)

    @inbounds for I in 1:nSlaterorb, J in 1:nSlaterorb
        p, q = Slaterorb_pairs[I]
        k, l = Slaterorb_pairs[J]

        sgn = isodd(p + q + k + l) ? -1 : 1
        C[I, J] = sgn * det(S[Slaterorb_keep[I], Slaterorb_keep[J]])
    end

    return C, Slaterorb_pairs
end


# Convert Tgrid::Vector{Vector{Int}} into two flat Int arrays.
# This avoids repeatedly doing Tgrid[a][1], Tgrid[a][2] in hot loops.
function Tgrid_components(Tgrid)
    M = length(Tgrid)

    Tn1 = Vector{Int}(undef, M)
    Tn2 = Vector{Int}(undef, M)

    @inbounds for a in 1:M
        Tn1[a] = Tgrid[a][1]
        Tn2[a] = Tgrid[a][2]
    end

    return Tn1, Tn2
end




# Group basisorb pairs (a,b) by total momentum T_a + T_b.
#
# Each returned group is a Vector{Tuple{Int,Int}} of basisorb pairs.
#
# If (a,b) and (c,d) are in the same group, then
#
#       T_a + T_b = T_c + T_d
#
# so momentum conservation is automatic.
function make_basisorb_total_groups(Tn1::AbstractVector{<:Integer},
                                    Tn2::AbstractVector{<:Integer})
    M = length(Tn1)

    group_dict = Dict{Tuple{Int, Int}, Vector{Tuple{Int, Int}}}()

    @inbounds for a in 1:M, b in 1:M
        a == b && continue

        total_key = (Tn1[a] + Tn1[b], Tn2[a] + Tn2[b])

        push!(get!(group_dict, total_key, Tuple{Int, Int}[]), (a, b))
    end

    return collect(values(group_dict))
end



# For a block of basisorb pairs (a,b), build the Slater-amplitude matrix
#
#   A[I, col] = X[a,p] X[b,q] - X[a,q] X[b,p]
#
# where
#
#   Slaterorb_pairs[I] = (p,q)
#   basisorb_group[col] = (a,b)
#
# Therefore A is nSlaterorb × nbasisorb_in_this_block.
function Slaterorb_amplitude_matrix(basisorb_group, Xmatrix, Slaterorb_pairs)
    nSlaterorb = length(Slaterorb_pairs)
    nbasisorb = length(basisorb_group)

    A = zeros(ComplexF64, nSlaterorb, nbasisorb)

    @inbounds for col in 1:nbasisorb
        a, b = basisorb_group[col]

        for I in 1:nSlaterorb
            p, q = Slaterorb_pairs[I]

            A[I, col] =
                Xmatrix[a, p] * Xmatrix[b, q] -
                Xmatrix[a, q] * Xmatrix[b, p]
        end
    end

    return A
end



# Precompute spinor overlaps:
#
#   Ω[a,c] = spinor_set[a]' * spinor_set[c]
#
# This is small for M ~ 600.
function spinor_overlap_matrix(spinor_set)
    M = length(spinor_set)

    Ω = zeros(ComplexF64, M, M)

    @inbounds for c in 1:M
        for a in 1:M
            Ω[a, c] = dot(spinor_set[a], spinor_set[c])
        end
    end

    return Ω
end



# ------------------------------------------------------------
# Main function:
#
# Extract G(q), where q = T_a - T_c.
#
# The returned object is:
#
#   Gq[(n1,n2)]
#
# such that for any bare interaction V(q),
#
#   E = 1/2 * sum_q V(q) * G(q)
#
# No S^{-1} is used.
# No V_{abcd} tensor is constructed.
# ------------------------------------------------------------


function q_resolved_pair_weight(
    Dmatrix,
    Fmatrix,
    precomp::QResolvedPrecomp;
    skip_q0::Bool = false,
)
    M, Nelectron = size(Dmatrix)

    @assert size(Fmatrix, 1) == M
    @assert size(Fmatrix, 2) == Nelectron
    
  

    # Cross-overlap matrix between the two Slater determinants.
    # This may be singular. We never invert it.
    S = Dmatrix' * Fmatrix

    # C is the singular-safe Slaterorb cofactor matrix.
    C, Slaterorb_pairs = Slaterorb_second_minor_matrix(S)

    Tn1 = precomp.Tn1
    Tn2 = precomp.Tn2
    @assert length(Tn1) == M
    basisorb_total_groups = precomp.basisorb_total_groups
    Ω = precomp.Ω

    # Thread-local dictionaries. We merge them after the threaded loop.
    partial_Gq = Vector{Dict{Tuple{Int, Int}, ComplexF64}}(
        undef,
        length(basisorb_total_groups),
    )

    for ig in eachindex(basisorb_total_groups)
        basisorb_group = basisorb_total_groups[ig]

        # AD[:,P] corresponds to basisorb P = (a,b) in D.
        AD = Slaterorb_amplitude_matrix(
            basisorb_group,
            Dmatrix,
            Slaterorb_pairs,
        )

        # AF[:,Q] corresponds to basisorb Q = (c,d) in F.
        AF = Dmatrix === Fmatrix ? AD :
             Slaterorb_amplitude_matrix(
                 basisorb_group,
                 Fmatrix,
                 Slaterorb_pairs,
             )

        # manybody_block[P,Q] =
        #
        #   A_D(a,b)^dagger C A_F(c,d)
        #
        # where:
        #   basisorb_group[P] = (a,b)
        #   basisorb_group[Q] = (c,d)
        manybody_block = AD' * (C * AF)

        local_Gq = Dict{Tuple{Int, Int}, ComplexF64}()

        @inbounds for P in eachindex(basisorb_group)
            a, b = basisorb_group[P]

            for Q in eachindex(basisorb_group)
                c, d = basisorb_group[Q]

                qkey = (Tn1[a] - Tn1[c], Tn2[a] - Tn2[c])

                if skip_q0 && qkey == (0, 0)
                    continue
                end

                form_factor = Ω[a, c] * Ω[b, d]

                local_Gq[qkey] =
                    get(local_Gq, qkey, 0.0 + 0.0im) +
                    form_factor * manybody_block[P, Q]
            end
        end

        partial_Gq[ig] = local_Gq
    end

    # Merge all thread-local G(q) dictionaries.
    Gq = Dict{Tuple{Int, Int}, ComplexF64}()

    for local_Gq in partial_Gq
        for (qkey, value) in local_Gq
            Gq[qkey] = get(Gq, qkey, 0.0 + 0.0im) + value
        end
    end

    return Gq,S
end


function Slaterorb_first_minor_matrix(S::AbstractMatrix)
    N = size(S, 1)
    @assert size(S, 2) == N

    C = zeros(ComplexF64, N, N)

    if N == 1
        C[1, 1] = 1.0 + 0.0im
        return C
    end

    keep = [[r for r in 1:N if r != p] for p in 1:N]

    @inbounds for p in 1:N, k in 1:N
        sgn = isodd(p + k) ? -1 : 1
        C[p, k] = sgn * det(S[keep[p], keep[k]])
    end

    return C
end


function onebody_transition_rho(Dmatrix, Fmatrix)
    # rho[a,b] = <Phi_D| c_b^dagger c_a |Phi_F>

    S = Dmatrix' * Fmatrix

    # C[p,k] = (-1)^(p+k) det S^(p|k)
    C = Slaterorb_first_minor_matrix(S)

    rho = Fmatrix * transpose(C) * Dmatrix'

    return rho, S
end




function main_func(args::Vector{Float64})

    flux1=args[1]
    flux2=args[2]
    NL=Int(args[3])
    moiream=args[4]
    N1=Int(args[5])
    N2=Int(args[6])
    N1f=Int(args[7])
    N2f=Int(args[8])
    xi00_re=args[9]
    xi00_im=args[10]
    grid_cutoff=args[11]
    type=Int(args[12])
    file_pos=Int(args[13])


    scratch_dir = ENV["SCRATCH"]
   






    topo_sectors = 0:2
    Nvec = length(topo_sectors)

    spinor_set_all =
        Vector{Vector{Vector{ComplexF64}}}(undef, Nvec)

    params_all = Vector{Any}(undef, Nvec)

    overall_mag_list_all =
        Vector{Vector{ComplexF64}}(undef, Nvec)

    orbital_basis_orthogonal_list_all =
        Vector{Vector{Matrix{ComplexF64}}}(undef, Nvec)

    for (s, topo_sec_s) in enumerate(topo_sectors)
        filepath_s=joinpath(scratch_dir, "constrcut_twothirds_skv/data_output$(Int(args[13]))/$(args[1])f1$(args[2])f2$(args[3])NL$(args[4])am$(args[5])N1$(args[6])N2$(args[7])N1f$(args[8])N2f$(args[9])xir$(args[10])xii$(args[11])grid$(topo_sec_s)topo$(args[12])type.jld2")
 


        isfile(filepath_s) || error("No file exists: $filepath_s")

        st_s = load(filepath_s)

        spinor_set_all[s] =
            st_s["spinor_set"]

        params_all[s] =
            st_s["params"]

        overall_mag_list_all[s] =
            st_s["overall_mag_list"]

        orbital_basis_orthogonal_list_all[s] =
            st_s["orbital_basis_orthogonal_list"]
    end



    spinor_set_ref = spinor_set_all[1]

    for s in 2:Nvec
        @assert length(spinor_set_all[s]) == length(spinor_set_ref)

        for i in eachindex(spinor_set_ref)
            @assert isapprox(
                spinor_set_all[s][i],
                spinor_set_ref[i];
                atol = 1e-10,
                rtol = 1e-10,
            )
        end
    end

    spinor_set = spinor_set_ref


    for s in 2:Nvec
    @assert params_all[s].Tgrid == params_all[1].Tgrid
    end

    Tgrid = params_all[1].Tgrid




  
    
    precomp = QResolvedPrecomp(Tgrid, spinor_set)

  
  

        # ------------------------------------------------------------
        # Raw transition objects in the three-state space
        # ------------------------------------------------------------

        deno_raw = zeros(ComplexF64, Nvec, Nvec)

        M_basis = size(
            orbital_basis_orthogonal_list_all[1][1],
            1,
        )

        Nelectron = size(
            orbital_basis_orthogonal_list_all[1][1],
            2,
        )

        # rho_raw[m,n,a,b] = <Psi_a|c_n^dagger c_m|Psi_b>
        rho_raw = zeros(
            ComplexF64,
            M_basis,
            M_basis,
            Nvec,
            Nvec,
        )

        # Gq_raw[q][a,b] = <Psi_a|G(q)|Psi_b>
        Gq_raw = Dict{
            Tuple{Int,Int},
            Matrix{ComplexF64}
        }()

    for s in 1:Nvec
        @assert length(overall_mag_list_all[s]) ==
                length(orbital_basis_orthogonal_list_all[s])

        for Dmatrix in orbital_basis_orthogonal_list_all[s]
            @assert size(Dmatrix) == (M_basis, Nelectron)
        end
    end
     

            # ------------------------------------------------------------
        # deno_raw[a,b] = <Psi_a|Psi_b>
        # ------------------------------------------------------------

        for a in 1:Nvec
            coeff_a = overall_mag_list_all[a]
            dets_a = orbital_basis_orthogonal_list_all[a]

            for b in 1:Nvec
                coeff_b = overall_mag_list_all[b]
                dets_b = orbital_basis_orthogonal_list_all[b]

                deno_ab = 0.0 + 0.0im

                for ia in eachindex(coeff_a)
                    Dmatrix = dets_a[ia]

                    for ib in eachindex(coeff_b)
                        Fmatrix = dets_b[ib]

                        coeff_ab =
                            conj(coeff_a[ia]) *
                            coeff_b[ib]

                        deno_ab +=
                            coeff_ab *
                            det(Dmatrix' * Fmatrix)
                    end
                end

                deno_raw[a, b] = deno_ab
            end
        end
        @assert isapprox(
            deno_raw,
            deno_raw';
            atol = 1e-9,
            rtol = 1e-9,
        )
        # Remove only numerical non-Hermiticity.
        deno_raw .= (deno_raw + deno_raw') / 2

        println("Raw many-body overlap matrix:")
        display(deno_raw)




        
        # ------------------------------------------------------------
        # rho_raw[:,:,a,b] = <Psi_a|c^dagger c|Psi_b>
        # ------------------------------------------------------------

        for a in 1:Nvec
            coeff_a = overall_mag_list_all[a]
            dets_a = orbital_basis_orthogonal_list_all[a]

            for b in 1:Nvec
                coeff_b = overall_mag_list_all[b]
                dets_b = orbital_basis_orthogonal_list_all[b]

                rho_ab =
                    zeros(ComplexF64, M_basis, M_basis)

                for ia in eachindex(coeff_a)
                    Dmatrix = dets_a[ia]

                    for ib in eachindex(coeff_b)
                        Fmatrix = dets_b[ib]

                        coeff_ab =
                            conj(coeff_a[ia]) *
                            coeff_b[ib]

                        rho_ij, _ =
                            onebody_transition_rho(
                                Dmatrix,
                                Fmatrix,
                            )

                        rho_ab .+= coeff_ab .* rho_ij
                    end
                end

                @views rho_raw[:, :, a, b] .= rho_ab
            end
        end


   for a in 1:Nvec, b in 1:Nvec
        @assert isapprox(
            rho_raw[:, :, a, b]',
            rho_raw[:, :, b, a];
            atol = 1e-9,
            rtol = 1e-9,
        )
    end

   println("Julia threads = ", Threads.nthreads())


        # ------------------------------------------------------------
        # Gq_raw[q][a,b] = <Psi_a|G(q)|Psi_b>
        # ------------------------------------------------------------

        println("Julia threads = ", Threads.nthreads())

        print_lock = ReentrantLock()

        for a in 1:Nvec
            coeff_a = overall_mag_list_all[a]
            dets_a = orbital_basis_orthogonal_list_all[a]

            for b in 1:Nvec
                coeff_b = overall_mag_list_all[b]
                dets_b = orbital_basis_orthogonal_list_all[b]

                local_Gq_list =
                    Matrix{
                        Dict{Tuple{Int,Int},ComplexF64}
                    }(
                        undef,
                        length(coeff_a),
                        length(coeff_b),
                    )

                Threads.@threads for ia in eachindex(coeff_a)
                    lock(print_lock) do
                        println(
                            "state pair ($a,$b), determinant " *
                            "$ia/$(length(coeff_a)), " *
                            "thread $(Threads.threadid())",
                        )
                        flush(stdout)
                    end

                    Dmatrix = dets_a[ia]

                    for ib in eachindex(coeff_b)
                        Fmatrix = dets_b[ib]

                        local_Gq_list[ia, ib], _ =
                            q_resolved_pair_weight(
                                Dmatrix,
                                Fmatrix,
                                precomp;
                                skip_q0 = false,
                            )
                    end
                end

                # Sum the Slater-determinant contributions into the
                # (a,b) entry of every q-resolved 3×3 matrix.
                for ia in eachindex(coeff_a)
                    for ib in eachindex(coeff_b)
                        coeff_ab =
                            conj(coeff_a[ia]) *
                            coeff_b[ib]

                        for (qkey, value) in local_Gq_list[ia, ib]
                            Gmatrix = get!(
                                Gq_raw,
                                qkey,
                            ) do
                                zeros(
                                    ComplexF64,
                                    Nvec,
                                    Nvec,
                                )
                            end

                            Gmatrix[a, b] +=
                                coeff_ab * value
                        end
                    end
                end
            end
        end


    orthogonalize_matrix =
        metric_gram_schmidt(deno_raw)


    deno_orth =
        orthogonalize_matrix' *
        deno_raw *
        orthogonalize_matrix

    @assert isapprox(
        deno_orth,
        Matrix{ComplexF64}(I, Nvec, Nvec);
        atol = 1e-10,
        rtol = 1e-10,
    )

    println("Orthogonalization matrix:")
    display(orthogonalize_matrix)

    println("Overlap after orthogonalization:")
    display(deno_orth)


    # ------------------------------------------------------------
    # Transform G(q) into the orthonormal state basis
    # ------------------------------------------------------------

    Gq_orth = Dict{
        Tuple{Int,Int},
        Matrix{ComplexF64}
    }()

    for (qkey, Gmatrix_raw) in Gq_raw
        Gq_orth[qkey] =
            orthogonalize_matrix' *
            Gmatrix_raw *
            orthogonalize_matrix
    end


        # ------------------------------------------------------------
    # Transform rho into the orthonormal many-body basis
    #
    # For each fixed one-body index pair (m,n),
    #
    # rho_orth[m,n,:,:] =
    #     A' * rho_raw[m,n,:,:] * A
    # ------------------------------------------------------------

    rho_orth = zeros(
        ComplexF64,
        M_basis,
        M_basis,
        Nvec,
        Nvec,
    )

    for m in 1:M_basis
        for n in 1:M_basis
            rho_state_raw =
                Matrix(@view rho_raw[m, n, :, :])

            @views rho_orth[m, n, :, :] .=
                orthogonalize_matrix' *
                rho_state_raw *
                orthogonalize_matrix
        end
    end




        for s in 1:Nvec
        @assert isapprox(
            tr(@view rho_orth[:, :, s, s]),
            Nelectron;
            atol = 1e-8,
            rtol = 1e-8,
        )
       end


        for (qkey, Gq_matrix) in Gq_orth
        minus_qkey = (-qkey[1], -qkey[2])

        if haskey(Gq_orth, minus_qkey)
            @assert isapprox(
                Gq_matrix',
                Gq_orth[minus_qkey];
                atol = 1e-8,
                rtol = 1e-8,
            )
        end
    end



    savepath=joinpath(scratch_dir, "constrcut_twothirds_skv/data_output$(Int(args[13]))/onetwobody/onetwobody_$(args[1])f1$(args[2])f2$(args[3])NL$(args[4])am$(args[5])N1$(args[6])N2$(args[7])N1f$(args[8])N2f$(args[9])xir$(args[10])xii$(args[11])grid$(args[12])type.jld2")
    jldsave(
        savepath;

        # --------------------------------------------------------
        # Raw transition data in the original three-state basis
        # --------------------------------------------------------

        deno_raw = deno_raw,
        rho_raw = rho_raw,
        Gq_raw = Gq_raw,

        # --------------------------------------------------------
        # Orthogonalization matrix
        #
        # |Psi_orth[s]> =
        #     sum_a |Psi_raw[a]> * orthogonalize_matrix[a,s]
        # --------------------------------------------------------

        orthogonalize_matrix = orthogonalize_matrix,

        # --------------------------------------------------------
        # Full transition data in the orthonormal three-state basis
        # Includes both diagonal and off-diagonal matrix elements
        # --------------------------------------------------------

        deno_orth = deno_orth,
        rho_orth = rho_orth,
        Gq_orth = Gq_orth,

        # --------------------------------------------------------
        # Underlying state information
        # --------------------------------------------------------

        spinor_set = spinor_set,
        params_all = params_all,
        overall_mag_list_all = overall_mag_list_all,
        orbital_basis_orthogonal_list_all =
            orbital_basis_orthogonal_list_all,

        topo_sectors = collect(topo_sectors),
    )

  return nothing


end


# ------------------------------------------------------------
# Metric Gram-Schmidt using deno_raw as the inner product
# ------------------------------------------------------------

function metric_gram_schmidt(
    overlap::AbstractMatrix{<:Complex},
)
    N = size(overlap, 1)

    @assert size(overlap, 2) == N
    @assert abs(det(overlap))>10^(-8)

    A = zeros(ComplexF64, N, N)

    for s in 1:N
        # Begin with raw state s.
        v = zeros(ComplexF64, N)
        v[s] = 1.0 + 0.0im

        # Remove projections onto the previously constructed states.
        for t in 1:s-1
            q = @view A[:, t]

            projection =
                dot(q, overlap * v)

            v .-= projection .* q
        end

        norm2 =
            real(dot(v, overlap * v))

        @assert norm2 > 0 """
        The raw many-body states are linearly dependent.
        Gram-Schmidt failed at state $s.
        norm² = $norm2
        """

        A[:, s] .= v ./ sqrt(norm2)
    end

    return A
end

