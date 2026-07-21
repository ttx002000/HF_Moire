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
    topo_sec=Int(args[12])
    type=Int(args[13])
    file_pos=Int(args[14])


    scratch_dir = ENV["SCRATCH"]
    filepath=joinpath(scratch_dir, "constrcut_twothirds_skv/data_output$(Int(args[14]))/$(args[1])f1$(args[2])f2$(args[3])NL$(args[4])am$(args[5])N1$(args[6])N2$(args[7])N1f$(args[8])N2f$(args[9])xir$(args[10])xii$(args[11])grid$(args[12])topo$(args[13])type.jld2")
 
    if isfile(filepath)
        st=load(filepath)
    else
        error("no file exists")
    end



    spinor_set=st["spinor_set"]
    params=st["params"]
    orbital_basis_orthogonal_list=st["orbital_basis_orthogonal_list"]
    overall_mag_list=st["overall_mag_list"]
    possible_config=st["possible_config"]
    vac_fac_list=st["vac_fac_list"]
    orbital_Rmatrix_list=st["orbital_Rmatrix_list"]
     orbital_basis_norm_list=st["orbital_basis_norm_list"]
      orbital_basis_orthogonal_list=st["orbital_basis_orthogonal_list"]
      M_eta_list=st["M_eta_list"]
    Mkkmatrix_list=st["Mkkmatrix_list"]
    N_coeff_list=st["N_coeff_list"]



    Tgrid=params.Tgrid
    
    precomp = QResolvedPrecomp(Tgrid, spinor_set)

  
  

    deno=0.0+0.0*im


    for ja in eachindex(overall_mag_list)
        for  jb in eachindex(overall_mag_list)
            coeff=overall_mag_list[ja]'*overall_mag_list[jb]
            Dmatrix=orbital_basis_orthogonal[ja]
            Fmatrix=orbital_basis_orthogonal[jb]
            deno+=det(Dmatrix'*Fmatrix)*coeff


        end
    end




    final_G_q=Dict{Tuple{Int, Int}, ComplexF64}()

    local_Gq_list=Matrix{Dict{Tuple{Int,Int},ComplexF64}}(undef,length(overall_mag_list),length(overall_mag_list))
    print_lock = ReentrantLock()s
    Threads.@threads for ja in eachindex(overall_mag_list)

    
        for  jb in eachindex(overall_mag_list)
                lock(print_lock) do
                        println("starting $(ja)/$(length(overall_mag_list)), $(jb)/$(length(overall_mag_list))")
                        flush(stdout)
                    end
            Dmatrix=orbital_basis_orthogonal[ja]
            Fmatrix=orbital_basis_orthogonal[jb]
            local_Gq_list[ja,jb],_=q_resolved_pair_weight(
                                                    Dmatrix,
                                                    Fmatrix,
                                                    precomp;
                                                    skip_q0 = false,
                                                )


        end
    end
   
    for ja in eachindex(overall_mag_list), jb in eachindex(overall_mag_list)
        coeff=overall_mag_list[ja]'*overall_mag_list[jb]
        for (qkey,vals) in local_Gq_list[ja,jb]
                final_G_q[qkey] = (get(final_G_q, qkey, 0.0 + 0.0im) + vals*coeff/deno)
        end

    end




        # ------------------------------------------------------------
    # One-body density matrix:
    #
    # rho_matrix[a,b] = <Psi| c_b^dagger c_a |Psi> / <Psi|Psi>
    # ------------------------------------------------------------

    M_basis = size(orbital_basis_orthogonal[1], 1)
    rho_num = zeros(ComplexF64, M_basis, M_basis)
  
    for ja in eachindex(overall_mag_list), jb in eachindex(overall_mag_list)
        coeff = overall_mag_list[ja]' * overall_mag_list[jb]

        Dmatrix = orbital_basis_orthogonal[ja]
        Fmatrix = orbital_basis_orthogonal[jb]

        rho_ij, _ = onebody_transition_rho(Dmatrix, Fmatrix)

        rho_num .+= coeff .* rho_ij
     
    end

    rho_matrix = rho_num ./ deno





    savepath=joinpath(scratch_dir, "constrcut_twothirds_skv/data_output$(Int(args[14]))/onetwobody/onetwobody_$(args[1])f1$(args[2])f2$(args[3])NL$(args[4])am$(args[5])N1$(args[6])N2$(args[7])N1f$(args[8])N2f$(args[9])xir$(args[10])xii$(args[11])grid$(args[12])topo$(args[13])type.jld2")
 
    jldsave(savepath,
     final_G_q=final_G_q, rho_matrix=rho_matrix,deno=deno,
     params=params,overall_mag_list=overall_mag_list,spinor_set=spinor_set,
            possible_config=possible_config, vac_fac_list= vac_fac_list,orbital_Rmatrix_list=orbital_Rmatrix_list,
               orbital_basis_norm_list=orbital_basis_norm_list, orbital_basis_orthogonal_list=orbital_basis_orthogonal_list,M_eta_list=M_eta_list,
               Mkkmatrix_list=Mkkmatrix_list,N_coeff_list=N_coeff_list)

  return nothing


end