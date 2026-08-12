using EllipticFunctions
using LinearAlgebra
using Plots
using JLD2,StaticArrays
using Combinatorics
using SparseArrays
using Base.Threads
BLAS.set_num_threads(1)
#I think I forgot to orthonormalize those state

function c_dot(z1::ComplexF64,z2::ComplexF64)
    return 1/2*(z1*conj(z2)+z2*conj(z1))
end
function c_cross(z1::ComplexF64,z2::ComplexF64)
   return 1/2*im*(z1*conj(z2)-z2*conj(z1))
end

function wrapparallel(k::Complex, b1::Complex, b2::Complex)
    A = [real(b1) real(b2);
         imag(b1) imag(b2)]
    rhs = [real(k), imag(k)]
    uv = A \ rhs
    u, v = uv[1], uv[2]

    frac(t) = t - round(t)   # matches Mathematica t - Round[t]
    return frac(u)*b1 + frac(v)*b2
end


function wrap_n1n2(k::Complex, b1::Complex, b2::Complex)
    A = [real(b1) real(b2);
         imag(b1) imag(b2)]
    rhs = [real(k), imag(k)]
    uv = A \ rhs
    u, v = uv[1], uv[2]

    frac(t) = t - round(t)   # matches Mathematica t - Round[t]
    return frac(u)*b1 + frac(v)*b2, Int(round(u)), Int(round(v))
end





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



function torusSigma(z,L1,L2,Lb,bigA)
    w_z,n1,n2=wrap_n1n2(z,L1/Lb,L2/Lb)

    pv=(-1)^(n1+n2)*wsigma(w_z, omega=(L1/(2*Lb), L2/(2*Lb)))*exp(-bigA*w_z^2)
    ev=(L1*conj(L1)*n1^2+L2*conj(L2)*n2^2)/(4*Lb^2)+(conj(L2)*n2+conj(L1)*n1)/(2*Lb)*w_z+n1*n2*conj(L2)*L1/(2*Lb^2)
    
    return sigma_related(pv,ev)

end


function modsigma(z,a1,a2,lb,α)
    w_z,n1,n2=wrap_n1n2(z,a1/lb,a2/lb)
    pv=(-1)^(n1+n2)*wsigma(w_z, omega=(a1/(2*lb), a2/(2*lb)))*exp(-α*w_z^2)
    ev=(a1*conj(a1)*n1^2+a2*conj(a2)*n2^2)/(4*lb^2)+(conj(a2)*n2+conj(a1)*n1)/(2*lb)*w_z+n1*n2*conj(a2)*a1/(2*lb^2)

    return sigma_related(pv,ev)
end





function landaulevel(z,k,a1,a2,lb,α,b1,b2)
    ss=modsigma(z/lb+im*lb*(k-b1/2-b2/2),a1,a2,lb,α)
    ev_add=im/2*(conj(k-b1/2-b2/2))*z-z*conj(z)/(4*lb^2)-lb^2/4*conj(k)*k+lb^2*k*(conj(b1+b2)/4)

    return sigma_related(ss.pv,ss.ev+ev_add)
end


function complex_vec(vv::Vector{Float64})
   return vv[1]+im*vv[2]
end



struct GqPrecomp
    qkeys::Vector{Tuple{Int,Int}}
    qindex::Dict{Tuple{Int,Int},Int}
    q_pairs::Vector{Vector{Tuple{Int,Int}}}
    minus_q_index::Vector{Int}

    Ω::Matrix{ComplexF64}

    Slaterorb_pairs::Vector{Tuple{Int,Int}}
    Slaterorb_keep::Vector{Vector{Int}}
    Slaterorb_first_keep::Vector{Vector{Int}}
end

function GqPrecomp(Tgrid, spinor_set, Nelectron)
    M = length(Tgrid)
    @assert length(spinor_set) == M

    # --------------------------------------------------------
    # Ω[a,c] = <s_a|s_c>
    # --------------------------------------------------------

    Ω = zeros(ComplexF64, M, M)

    for a in 1:M, c in 1:M
        Ω[a,c] = dot(spinor_set[a], spinor_set[c])
    end

    # --------------------------------------------------------
    # Group every ordered basis pair (a,c) by
    #
    # q = T_a - T_c.
    #
    # Across all q groups there are exactly M^2 entries.
    # --------------------------------------------------------

    qpair_dict =
        Dict{Tuple{Int,Int},Vector{Tuple{Int,Int}}}()

    for a in 1:M, c in 1:M
        qkey = (
            Tgrid[a][1] - Tgrid[c][1],
            Tgrid[a][2] - Tgrid[c][2],
        )

        push!(
            get!(qpair_dict, qkey, Tuple{Int,Int}[]),
            (a,c),
        )
    end

    qkeys = sort!(collect(keys(qpair_dict)))
    q_pairs = [qpair_dict[qkey] for qkey in qkeys]

    qindex = Dict{Tuple{Int,Int},Int}()

    for iq in eachindex(qkeys)
        qindex[qkeys[iq]] = iq
    end

    minus_q_index = Vector{Int}(undef, length(qkeys))

    for iq in eachindex(qkeys)
        q1, q2 = qkeys[iq]
        minus_q_index[iq] = qindex[(-q1,-q2)]
    end

    # --------------------------------------------------------
    # Occupied-orbital pairs (p,r), p < r.
    #
    # These depend only on Nelectron, so don't rebuild them
    # for every pair of Slater determinants.
    # --------------------------------------------------------

    Slaterorb_pairs = [
        (p,r)
        for p in 1:Nelectron-1
        for r in p+1:Nelectron
    ]

    Slaterorb_keep = [
        [x for x in 1:Nelectron if x != p && x != r]
        for (p,r) in Slaterorb_pairs
    ]

    Slaterorb_first_keep = [
    [x for x in 1:Nelectron if x != p]
    for p in 1:Nelectron]


    return GqPrecomp(
        qkeys,
        qindex,
        q_pairs,
        minus_q_index,
        Ω,
        Slaterorb_pairs,
        Slaterorb_keep,
        Slaterorb_first_keep,
    )
end




function Slaterorb_first_minor_matrix(
    S::AbstractMatrix,
    precomp::GqPrecomp,
)
    N = size(S, 1)
    @assert size(S, 2) == N

    keep = precomp.Slaterorb_first_keep

    C = zeros(ComplexF64, N, N)

    if N == 1
        C[1,1] = 1.0 + 0.0im
        return C
    end

    for p in 1:N, m in 1:N
        sgn = isodd(p + m) ? -1 : 1
        C[p,m] =
            sgn * det(S[keep[p], keep[m]])
    end

    return C
end




function Slaterorb_second_minor_matrix(
    S::AbstractMatrix,
    precomp::GqPrecomp,
)
    pairs = precomp.Slaterorb_pairs
    keep = precomp.Slaterorb_keep

    C = zeros(ComplexF64, length(pairs), length(pairs))

    for I in eachindex(pairs), J in eachindex(pairs)
        p, r = pairs[I]
        m, n = pairs[J]

        sgn = isodd(p + r + m + n) ? -1 : 1
        C[I,J] = sgn * det(S[keep[I], keep[J]])
    end

    return C
end

mutable struct SlaterPairWorkspace
    Xq::Array{ComplexF64,3}
    Gq::Vector{ComplexF64}
end


function SlaterPairWorkspace(Nelectron, Nq)
    Xq = zeros(ComplexF64, Nelectron, Nelectron, Nq)
    Gq = zeros(ComplexF64, Nq)

    return SlaterPairWorkspace(Xq, Gq)
end


function build_Xq!(
    workspace::SlaterPairWorkspace,
    Dmatrix,
    Fmatrix,
    precomp::GqPrecomp,
)
    Xq = workspace.Xq
    fill!(Xq, 0.0 + 0.0im)

    Nelectron = size(Dmatrix, 2)

    for iq in eachindex(precomp.qkeys)
        X = @view Xq[:,:,iq]

        for (a,c) in precomp.q_pairs[iq]
            ff = precomp.Ω[a,c]

            for p in 1:Nelectron
                left = conj(Dmatrix[a,p]) * ff

                for m in 1:Nelectron
                    X[p,m] += left * Fmatrix[c,m]
                end
            end
        end
    end

    return nothing
end

function slater_pair_Gq!(
    workspace::SlaterPairWorkspace,
    Dmatrix,
    Fmatrix,
    precomp::GqPrecomp;
    check_q0 = true,
)
    M, Nelectron = size(Dmatrix)
    @assert size(Fmatrix) == (M, Nelectron)

    S = Dmatrix' * Fmatrix
    overlap = det(S)

    C = Slaterorb_second_minor_matrix(S, precomp)
    C1 = Slaterorb_first_minor_matrix(S, precomp)

    build_Xq!(workspace, Dmatrix, Fmatrix, precomp)

    Xq = workspace.Xq
    Gq = workspace.Gq
    fill!(Gq, 0.0 + 0.0im)

    pairs = precomp.Slaterorb_pairs

    for iq in eachindex(precomp.qkeys)
        iminus = precomp.minus_q_index[iq]

        Xplus = @view Xq[:,:,iq]
        Xminus = @view Xq[:,:,iminus]

        value = 0.0 + 0.0im

        for I in eachindex(pairs)
            p, r = pairs[I]

            for J in eachindex(pairs)
                m, n = pairs[J]

                pair_value =
                    Xplus[p,m] * Xminus[r,n] -
                    Xplus[p,n] * Xminus[r,m] -
                    Xplus[r,m] * Xminus[p,n] +
                    Xplus[r,n] * Xminus[p,m]

                value += C[I,J] * pair_value
            end
        end

        Gq[iq] = value
    end

    if check_q0
        iq0 = precomp.qindex[(0,0)]
        expected_q0 = Nelectron * (Nelectron - 1) * overlap

        @assert isapprox(
            Gq[iq0],
            expected_q0;
            atol = 1e-7,
            rtol = 1e-7,
        )
    end

    return overlap, C1, C
end




function orthogonalize_Gq(
    Gq_raw_array,
    reconstruct_matrix,
)
    Nvec, Nvec2, Nq = size(Gq_raw_array)
    @assert Nvec == Nvec2
    @assert size(reconstruct_matrix) == (Nvec, Nvec)

    Gq_orth_array = zeros(ComplexF64, Nvec, Nvec, Nq)

    Threads.@threads :static for iq in 1:Nq
        Gq_raw = @view Gq_raw_array[:,:,iq]

        Gq_orth_array[:,:,iq] .=
            reconstruct_matrix' * Gq_raw * reconstruct_matrix
    end

    return Gq_orth_array
end


function orthogonalize_rho(
    rho_raw,
    reconstruct_matrix,
)
    M_basis, M_basis2, Nvec, Nvec2 = size(rho_raw)

    @assert M_basis == M_basis2
    @assert Nvec == Nvec2
    @assert size(reconstruct_matrix) == (Nvec, Nvec)

    rho_orth = zeros(
        ComplexF64,
        M_basis,
        M_basis,
        Nvec,
        Nvec,
    )

    Threads.@threads :static for m in 1:M_basis
        for n in 1:M_basis
            rho_state = @view rho_raw[m,n,:,:]

            @views rho_orth[m,n,:,:] .=
                reconstruct_matrix' *
                rho_state *
                reconstruct_matrix
        end
    end

    return rho_orth
end





function metric_gram_schmidt(overlap)
    N = size(overlap, 1)

    @assert size(overlap, 2) == N
    println(det(overlap),"determinant overlap")
    @assert abs(det(overlap)) > 1e-9

    A = zeros(ComplexF64, N, N)

    for s in 1:N
        v = zeros(ComplexF64, N)
        v[s] = 1.0 + 0.0im

        for t in 1:s-1
            q = @view A[:,t]
            projection = dot(q, overlap * v)
            v .-= projection .* q
        end

        norm2 = real(dot(v, overlap * v))

        @assert norm2 > 0 """
        The raw many-body states are linearly dependent.
        Gram-Schmidt failed at state $s.
        norm² = $norm2
        """

        A[:,s] .= v ./ sqrt(norm2)
    end

    return A
end


function calculate_block_overlap_Gq(
    orbital_basis_orthogonal_list,
    config_range,
    precomp::GqPrecomp,
)
    Nconfig_block = length(config_range)
    Nelectron = size(orbital_basis_orthogonal_list[first(config_range)], 2)
    Nq = length(precomp.qkeys)

    slater_overlap_block = zeros(ComplexF64, Nconfig_block, Nconfig_block)
    Gq_slater_block = zeros(ComplexF64, Nconfig_block, Nconfig_block, Nq)

    Norbital_pair = length(precomp.Slaterorb_pairs)

    second_minor_block = zeros(
        ComplexF64,
        Norbital_pair,
        Norbital_pair,
        Nconfig_block,
        Nconfig_block,
    )

    workspace = SlaterPairWorkspace(Nelectron, Nq)

    for local_i in 1:Nconfig_block
        global_i = config_range[local_i]
        Dmatrix = orbital_basis_orthogonal_list[global_i]

        for local_j in local_i:Nconfig_block
            global_j = config_range[local_j]
            Fmatrix = orbital_basis_orthogonal_list[global_j]

            overlap, _, C2 = slater_pair_Gq!(
                workspace,
                Dmatrix,
                Fmatrix,
                precomp,
            )
            @views second_minor_block[:,:,local_i,local_j] .= C2

            slater_overlap_block[local_i,local_j] = overlap

            for iq in 1:Nq
                Gq_slater_block[local_i,local_j,iq] = workspace.Gq[iq]
            end

            if local_i == local_j
                nothing
            else
                
                slater_overlap_block[local_j,local_i] = conj(overlap)
                @views second_minor_block[:,:,local_j,local_i] .= C2'

                for iq in 1:Nq
                    minus_iq = precomp.minus_q_index[iq]

                    Gq_slater_block[local_j,local_i,iq] =
                        conj(workspace.Gq[minus_iq])
                end
            end
        end
    end

    return (
        slater_overlap_block,
        Gq_slater_block,
        second_minor_block,
    )
end





function construct_manybody_rho_blockwise(
    orbital_basis_orthogonal_list,
    overall_mag_matrix,
    blocks,
    block_config_ranges,
    precomp::GqPrecomp,
)
    Nconfig = length(orbital_basis_orthogonal_list)
    M_basis, Nelectron = size(orbital_basis_orthogonal_list[1])
    Nvec = size(overall_mag_matrix, 2)
    Nblock = length(blocks)

    rho_raw = zeros(ComplexF64, M_basis, M_basis, Nvec, Nvec)

    first_minor_matrix = zeros(
        ComplexF64,
        Nelectron,
        Nelectron,
        Nconfig,
        Nconfig,
    )

    for block_a_index in 1:Nblock
        states_a = blocks[block_a_index]
        config_range_a = block_config_ranges[block_a_index]

        for block_b_index in block_a_index:Nblock
            states_b = blocks[block_b_index]
            config_range_b = block_config_ranges[block_b_index]

            rho_block = zeros(
                ComplexF64,
                M_basis,
                M_basis,
                length(states_a),
                length(states_b),
            )

            for global_i in config_range_a
                Dmatrix = orbital_basis_orthogonal_list[global_i]

                for global_j in config_range_b
                    Fmatrix = orbital_basis_orthogonal_list[global_j]

                    orbital_overlap = Dmatrix' * Fmatrix

                    first_minor =
                        Slaterorb_first_minor_matrix(
                            orbital_overlap,
                            precomp,
                        )

                    @views first_minor_matrix[:,:,global_i,global_j] .=
                        first_minor

                    if block_a_index == block_b_index
                        nothing
                    else
                        @views first_minor_matrix[:,:,global_j,global_i] .=
                            first_minor'
                    end

                    rho_slater_pair =
                        Fmatrix *
                        transpose(first_minor) *
                        Dmatrix'

                    for local_state_a in eachindex(states_a)
                        state_a = states_a[local_state_a]
                        coefficient_a =
                            overall_mag_matrix[global_i,state_a]

                        for local_state_b in eachindex(states_b)
                            state_b = states_b[local_state_b]
                            coefficient_b =
                                overall_mag_matrix[global_j,state_b]

                            weight =
                                conj(coefficient_a) *
                                coefficient_b

                            rho_state_pair =
                                @view rho_block[
                                    :,
                                    :,
                                    local_state_a,
                                    local_state_b,
                                ]

                            rho_state_pair .+=
                                weight .* rho_slater_pair
                        end
                    end
                end
            end

            for local_state_a in eachindex(states_a)
                state_a = states_a[local_state_a]

                for local_state_b in eachindex(states_b)
                    state_b = states_b[local_state_b]

                    rho_ab =
                        @view rho_block[
                            :,
                            :,
                            local_state_a,
                            local_state_b,
                        ]

                    rho_raw_ab =
                        @view rho_raw[:,:,state_a,state_b]

                    rho_raw_ab .= rho_ab

                    if block_a_index == block_b_index
                        nothing
                    else
                        rho_raw_ba =
                            @view rho_raw[:,:,state_b,state_a]

                        rho_raw_ba .= rho_ab'
                    end
                end
            end
        end
    end

    return rho_raw, first_minor_matrix
end


function do_Gq(args::Vector{Float64}
              )


     flux1=args[1]
     flux2=args[2]
     q1=args[3]
     q2=args[4]
     N1=Int(args[5])
     N2=Int(args[6])
     Npa_LL=Int(args[7])
     Nvec=Int(args[8])
     file_pos=Int(args[9])
     grid_cutoff=args[10]
     NL=Int(args[11])
     moiream=args[12]
     gaus_damp=args[13]




               




    θ = π/3;
    l1 = moiream;
    l2 = moiream;

    a1 = Complex(l1)
    a2 = l2*(cos(θ) + im*sin(θ))
    lb=(1/(2π)*1/2*abs(a1*conj(a2)-a2*conj(a1)))^(1/2)

    b1=-im*a2/lb^2
    b2=im*a1/lb^2


    Nphi=N1*N2
    Nelectron=Nphi-Npa_LL

    

    Lb=(N1*N2)^(1/2)*lb
    L1=N1*a1
    L2=N2*a2
    T1=b1/N1
    T2=b2/N2

 
    Deltatheta=(T1*flux1+T2*flux2)
    qq=q1*b1+q2*b2


    T1_vec=[real(T1),imag(T1)]
    T2_vec=[real(T2),imag(T2)]
    b1_vec=[real(b1),imag(b1)]
    b2_vec=[real(b2),imag(b2)]

    b1T=Int.(round.(inv([T1_vec T2_vec])*b1_vec))
    b2T=Int.(round.(inv([T1_vec T2_vec])*b2_vec))


    a1_vec=[real(a1),imag(a1)]
    a2_vec=[real(a2),imag(a2)]
    L1_vec=[real(L1),imag(L1)]
    L2_vec=[real(L2),imag(L2)]


    eta1 = wzeta(a1 / (2 * lb), omega=(a1 / (2 * lb), a2 / (2 * lb)))
    eta2 = wzeta(a2 / (2 * lb), omega=(a1 / (2 * lb), a2 / (2 * lb)))
    α = eta1 * lb / a1 - conj(a1) / (4 * a1)



    Tgrid=Vector{Int}[]
  
    cutoffstandard=norm(b1_vec)*grid_cutoff
    cutoff=5*Int(round.(max(cutoffstandard/norm(T1_vec),cutoffstandard/norm(T2_vec))))
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
    gg=ja*[1,0]+jb*[0,1]
    gg_vec=[T1_vec T2_vec]*gg

    if norm(gg_vec)<cutoffstandard
        push!(Tgrid,gg) 
    end

    end

    Tgrid_dict=Dict{Vector{Int64},Int64}()
    for ja in eachindex(Tgrid)
        Tgrid_dict[Tgrid[ja]]=ja
    end


    ggrid=Vector{Int}[]
    cutoffstandard=norm(b1_vec)*(grid_cutoff+2)
    cutoff=5*Int(round.(max(cutoffstandard/norm(b1_vec),cutoffstandard/norm(b2_vec))))
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        gg=ja*b1T+jb*b2T
        gg_vec=[T1_vec T2_vec]*gg

        if norm(gg_vec)<cutoffstandard
            push!(ggrid,gg) 
        end

    end

    ggrid_dict=Dict{Vector{Int64},Int64}()
    for ja in eachindex(ggrid)
        ggrid_dict[ggrid[ja]]=ja
    end





    scratch_dir = ENV["SCRATCH"]
    LL_filepath=joinpath(scratch_dir, "constrcut_twothirds_Lauhglin_ED_v2/data_output$(Int(args[9]))/$(args[1])f1$(args[2])f2$(args[3])q1$(args[4])q2$(args[5])N1$(args[6])N2$(args[7])Npa$(args[8])Nvec.jld2")


    LL_data=load(LL_filepath)
    LL_allowedq=LL_data["allowedq"]
    LL_T1=LL_data["T1"]
    LL_T2=LL_data["T2"]

    @assert norm(LL_T1/moiream-T1_vec)<10^(-9)
    @assert norm(LL_T2/moiream-T2_vec)<10^(-9)
    

    PH_eigvector=LL_data["PH_eigvector"]
    PH_state_can=LL_data["PH_state_can"]
    @assert length(LL_data["PH_eigvector"])==Nvec

    LL_k_set=[complex_vec([LL_T1 LL_T2]*LL_allowedq[ja])/moiream for ja in eachindex(LL_allowedq)]
    k_set=vec([wrapparallel(n1*T1+n2*T2,b1,b2) for n1 in 1:N1, n2 in 1:N2])

    LL_k_map=zeros(Int,length(k_set))
    for ja in eachindex(LL_k_map)
        ccount=0
        for jb in eachindex(k_set)
        if abs(k_set[jb]-wrapparallel(LL_k_set[ja],b1,b2))<10^(-9)
            ccount+=1
            LL_k_map[ja]=jb
            
        end
        
        end
        @assert ccount==1
    end

    @assert sort(LL_k_map)==collect(1:1:length(k_set))


    k_set_n=Vector{Int}[]
    for ja in eachindex(k_set)
    nn=round.(inv([T1_vec T2_vec])*[real(k_set[ja]), imag(k_set[ja])])
    push!(k_set_n,nn)
    end


    

    Bloch_states=zeros(ComplexF64,length(Tgrid),N1*N2)

    for k1i in eachindex(k_set)
        k1_n=k_set_n[k1i]
        k1=k_set[k1i]


            for gi in eachindex(ggrid)
                g_n=ggrid[gi]
                gg=g_n[1]*T1+g_n[2]*T2
                g_n_in_b=inv([b1_vec b2_vec])*[real(gg),imag(gg)]

                Rk2=-im*lb^2*(-qq)
                pp=get(Tgrid_dict,k1_n+g_n,0)
            
                ct_2=exp(im*π* g_n_in_b[1]*g_n_in_b[2])*exp(-lb^2/4*(conj(gg)*gg+2*(k1-Deltatheta)*conj(gg)))
                ct_3=exp(-im*c_dot(k1+gg-Deltatheta,Rk2))
                ct_4=1/get_spinor_norm(NL,[real(k1-Deltatheta+gg),imag(k1-Deltatheta+gg)])

                if pp≠0

                Bloch_states[pp,k1i]+=ct_2*ct_3*ct_4
                
                end
            
            
            end

    end

    # ------------------------------------------------------------
    # Future Gaussian spreading:
    #
    # for ja in eachindex(Tgrid)
    #     kk = Tgrid[ja][1]*T1 + Tgrid[ja][2]*T2 - Deltatheta
    #     Bloch_states[ja,:] .*= exp(-xi * lb^2 * abs2(kk))
    # end
    #
    # The columns must then be renormalized.
    # Nk is calculated afterward, so the many-body coefficients
    # will be modified consistently through the products of Nk.
    # ------------------------------------------------------------

    for ja in 1:N1*N2
        Bloch_states[:,ja] /= norm(Bloch_states[:,ja])
    end




    Nk=zeros(ComplexF64,N1*N2)
    for whichk in 1:N1*N2
        whichk_map=LL_k_map[whichk]
        success=0
        last_result=0.0
        for attemp in 1:10
            z_test=rand(ComplexF64)*abs(a1)
            t1=landaulevel(z_test,LL_k_set[whichk]-Deltatheta-qq,a1,a2,lb,α,b1,b2)
            t2=landaulevel(z_test,-qq,a1,a2,lb,α,b1,b2)
            t3=0.0
            for ja in eachindex(Tgrid)
            kk=Tgrid[ja][1]*T1+Tgrid[ja][2]*T2-Deltatheta
            t3+=Bloch_states[ja,whichk_map]*get_spinor_norm(NL,[real(kk),imag(kk)])*exp(im*c_dot(z_test,kk))
            end

            if abs(t1.pv*conj(t2.pv)*exp(t1.ev+conj(t2.ev)))>0.01 && abs(t3)>0.01
                Nk[whichk_map]=t1.pv*conj(t2.pv)*exp(t1.ev+conj(t2.ev))/t3
                if success==0
                    last_result=t1.pv*conj(t2.pv)*exp(t1.ev+conj(t2.ev))/t3
                    success+=1
                else
                    @assert abs(t1.pv*conj(t2.pv)*exp(t1.ev+conj(t2.ev))/t3-last_result)<10^(-8)
                    success+=1
                end

            end
        end
        @assert success>=2
    end

    
    for ja in eachindex(Tgrid)
            kk=T1*Tgrid[ja][1]+T2*Tgrid[ja][2]-Deltatheta
            ff=exp(-gaus_damp*c_dot(kk,kk)*lb^2*3/2)# I added this 3/2 just to have the same convention
        Bloch_states[ja,:]= Bloch_states[ja,:]*ff
    
    end

   renorm_record=zeros(Float64,N1*N2)
    for ja in 1:N1*N2
        renorm_record[ja]=norm(Bloch_states[:,ja])
        Bloch_states[:,ja]=Bloch_states[:,ja]/norm(Bloch_states[:,ja])
    end










    blocks = Vector{Vector{Int}}()
    used = falses(Nvec)

    for i in 1:Nvec
        used[i] && continue

        block = [
            j for j in 1:Nvec
            if PH_state_can[j] == PH_state_can[i]
        ]

        push!(blocks, block)
        used[block] .= true
    end

            



        Nblock = length(blocks)

        block_configurations =
            Vector{Vector{Vector{Int}}}(undef, Nblock)

        block_config_ranges =
            Vector{UnitRange{Int}}(undef, Nblock)

        possible_config = Vector{Vector{Int}}()

        config_start = 1
     

        

        for block_index in 1:Nblock
            block = blocks[block_index]
            reference_state = block[1]
            configurations = PH_state_can[reference_state]

            for state_index in block
                if PH_state_can[state_index] == configurations
                    nothing
                else
                    error("States in block $block do not share the same PH_state_can.")
                end
            end

            block_configurations[block_index] = configurations

            config_stop = config_start + length(configurations) - 1
            block_config_ranges[block_index] = config_start:config_stop

            append!(possible_config, copy.(configurations))

            config_start = config_stop + 1
        end

        Nconfig = length(possible_config)

        @assert config_start == Nconfig + 1
           
            
            orbital_basis_orthogonal_list =
            Vector{Matrix{ComplexF64}}(undef, Nconfig)

        for configuration_index in 1:Nconfig
            orbital_basis_orthogonal_list[configuration_index] =
                construct_Slatermatrix(
                    possible_config[configuration_index],
                    LL_k_map,
                    Bloch_states,
                )
        end

        Nk_factor_list = Vector{ComplexF64}(undef, Nconfig)

        for configuration_index in 1:Nconfig
            Nk_factor_list[configuration_index] =
                configuration_Nk_factor(
                    possible_config[configuration_index],
                    LL_k_map,
                    Nk,
                )
        end

        renorm_factor_list = Vector{Float64}(undef, Nconfig)

        for configuration_index in 1:Nconfig
            configuration = possible_config[configuration_index]
            renorm_factor = 1.0

            for LL_orbital_index in configuration
                orbital_index = LL_k_map[LL_orbital_index]
                renorm_factor *= renorm_record[orbital_index]
            end

            renorm_factor_list[configuration_index] = renorm_factor
        end



        overall_mag_matrix =
            zeros(ComplexF64, Nconfig, Nvec)

        for block_index in 1:Nblock
            block = blocks[block_index]
            configurations = block_configurations[block_index]
            config_range = block_config_ranges[block_index]

            for state_index in block
                if length(PH_eigvector[state_index]) == length(configurations)
                    nothing
                else
                    error(
                        "State $state_index has the wrong number of coefficients."
                    )
                end

                for local_config in eachindex(configurations)
                    global_config = config_range[local_config]
                    overall_mag_matrix[global_config,state_index] =PH_eigvector[state_index][local_config]*Nk_factor_list[global_config]*renorm_factor_list[global_config]
                end
            end
        end

        spinor_set = Vector{Vector{ComplexF64}}(undef, length(Tgrid))

        for ja in eachindex(Tgrid)
            kk = Tgrid[ja][1]*T1 + Tgrid[ja][2]*T2 - Deltatheta
            spinor_set[ja] = get_spinor(NL, [real(kk), imag(kk)])
        end

        precomp = GqPrecomp(Tgrid, spinor_set, Nelectron)

        Nq = length(precomp.qkeys)
        M_basis = length(Tgrid)

        deno_raw = zeros(ComplexF64, Nvec, Nvec)
        Gq_raw_array = zeros(ComplexF64, Nvec, Nvec, Nq)
        Norbital_pair = length(precomp.Slaterorb_pairs)
        second_minor_matrix = zeros(
            ComplexF64,
            Norbital_pair,
            Norbital_pair,
            Nconfig,
            Nconfig,
        )
      

        for block_index in 1:Nblock
            block = blocks[block_index]
            config_range = block_config_ranges[block_index]

            slater_overlap_block,
            Gq_slater_block,
            second_minor_block =
                calculate_block_overlap_Gq(
                    orbital_basis_orthogonal_list,
                    config_range,
                    precomp,
                )

            @views second_minor_matrix[
                :,
                :,
                config_range,
                config_range,
            ] .= second_minor_block

            coefficient_block =
                overall_mag_matrix[config_range,block]

            deno_raw[block,block] .=
                coefficient_block' *
                slater_overlap_block *
                coefficient_block

            for iq in 1:Nq
                Gq_slater_q = @view Gq_slater_block[:,:,iq]

                Gq_raw_array[block,block,iq] .=
                    coefficient_block' *
                    Gq_slater_q *
                    coefficient_block
            end
        end

        deno_raw .= (deno_raw + deno_raw') / 2

        if isapprox(
            deno_raw,
            deno_raw';
            atol=1e-9,
            rtol=1e-9)
            nothing
        else
            error("deno_raw is not Hermitian.")
        end

        for iq in 1:Nq
            minus_iq = precomp.minus_q_index[iq]

            if isapprox(
                @view(Gq_raw_array[:,:,iq])',
                @view(Gq_raw_array[:,:,minus_iq]);
                atol=1e-7,
                rtol=1e-7,
            )
                nothing
            else
                error("Gq_raw_array failed Hermiticity at q index $iq.")
            end
        end

        iq0 = precomp.qindex[(0,0)]

        if isapprox(
            @view(Gq_raw_array[:,:,iq0]),
            Nelectron * (Nelectron - 1) * deno_raw;
            atol=1e-7,
            rtol=1e-7,
        )
            nothing
        else
            error("The q = 0 G(q) check failed.")
        end


        rho_raw, first_minor_matrix =
        construct_manybody_rho_blockwise(
            orbital_basis_orthogonal_list,
            overall_mag_matrix,
            blocks,
            block_config_ranges,
            precomp,
        )


        for state_a in 1:Nvec
            for state_b in 1:Nvec
                rho_ab = @view rho_raw[:,:,state_a,state_b]
                rho_ba = @view rho_raw[:,:,state_b,state_a]

                if isapprox(
                    rho_ab,
                    rho_ba';
                    atol=1e-7,
                    rtol=1e-7,
                )
                    nothing
                else
                    error(
                        "rho_raw Hermiticity failed for states " *
                        "$state_a and $state_b."
                    )
                end
            end
        end

        for state_a in 1:Nvec
            for state_b in 1:Nvec
                calculated_trace =
                    tr(@view rho_raw[:,:,state_a,state_b])

                expected_trace =
                    Nelectron * deno_raw[state_a,state_b]

                if isapprox(
                    calculated_trace,
                    expected_trace;
                    atol=1e-7,
                    rtol=1e-7,
                )
                    nothing
                else
                    error(
                        "rho_raw trace failed for states " *
                        "$state_a and $state_b.\n" *
                        "calculated = $calculated_trace\n" *
                        "expected   = $expected_trace"
                    )
                end
            end
        end

        reconstruct_matrix = zeros(ComplexF64, Nvec, Nvec)

        for block_index in 1:Nblock
            block = blocks[block_index]
            overlap_block = Matrix(Hermitian(deno_raw[block,block]))

            reconstruct_block = metric_gram_schmidt(overlap_block)
            reconstruct_matrix[block,block] .= reconstruct_block
        end

        overall_mag_matrix_orth = overall_mag_matrix * reconstruct_matrix

        deno_orth =
            reconstruct_matrix' *
            deno_raw *
            reconstruct_matrix

        Gq_orth_array =
            orthogonalize_Gq(
                Gq_raw_array,
                reconstruct_matrix,
            )

        rho_orth =
            orthogonalize_rho(
                rho_raw,
                reconstruct_matrix,
            )

       identity_state = Matrix{ComplexF64}(I, Nvec, Nvec)

        if isapprox(
            deno_orth,
            identity_state;
            atol=1e-8,
            rtol=1e-8,
        )
            nothing
        else
            error(
                "Orthogonalization failed.\n" *
                "maximum error = $(maximum(abs.(deno_orth - identity_state)))"
            )
        end


        for iq in 1:Nq
            minus_iq = precomp.minus_q_index[iq]

            if isapprox(
                @view(Gq_orth_array[:,:,iq])',
                @view(Gq_orth_array[:,:,minus_iq]);
                atol=1e-7,
                rtol=1e-7,
            )
                nothing
            else
                error(
                    "Gq_orth_array failed Hermiticity at q index $iq."
                )
            end
        end


        if isapprox(
            @view(Gq_orth_array[:,:,iq0]),
            Nelectron * (Nelectron - 1) * deno_orth;
            atol=1e-7,
            rtol=1e-7,
        )
            nothing
        else
            error("The orthogonalized q = 0 G(q) check failed.")
        end

        params = (
            flux1=flux1,
            flux2=flux2,
            q1=q1,
            q2=q2,
            N1=N1,
            N2=N2,
            Npa_LL=Npa_LL,
            Nvec=Nvec,
            file_pos=file_pos,
            grid_cutoff=grid_cutoff,
            NL=NL,
            moiream=moiream,
            Nphi=Nphi,
            Nelectron=Nelectron,
            lb=lb,
            Lb=Lb,
            L1=L1,
            L2=L2,
            T1=T1,
            T2=T2,
            Deltatheta=Deltatheta,
            qq=qq,
            Tgrid=Tgrid,
        )


    return (
    params = params,

    blocks = blocks,
    block_config_ranges = block_config_ranges,
    possible_config = possible_config,

    orbital_basis_orthogonal_list =
        orbital_basis_orthogonal_list,

    overall_mag_matrix =
        overall_mag_matrix,

    reconstruct_matrix =
        reconstruct_matrix,

    first_minor_matrix =
        first_minor_matrix,

    deno_raw =
        deno_raw,

    Gq_raw_array =
        Gq_raw_array,

    rho_raw =
        rho_raw,

    deno_orth =
        deno_orth,

    Gq_orth_array =
        Gq_orth_array,

    rho_orth =
        rho_orth,

    qkeys =
        precomp.qkeys,

    qindex =
        precomp.qindex,

    minus_q_index =
        precomp.minus_q_index,

    spinor_set =
        spinor_set,

    second_minor_matrix =
    second_minor_matrix,
)




end

function configuration_Nk_factor(
    configuration,
    LL_k_map,
    Nk,
)
    Nk_factor = 1.0 + 0.0im

    for LL_orbital_index in configuration
        orbital_index = LL_k_map[LL_orbital_index]
        Nk_factor *= Nk[orbital_index]
    end

    return Nk_factor
end

function construct_Slatermatrix(
    configuration,
    LL_k_map,
    Bloch_states,
)
    Nelectron = length(configuration)
    M_basis = size(Bloch_states, 1)

    Slatermatrix = zeros(ComplexF64, M_basis, Nelectron)

    for electron_index in 1:Nelectron
        LL_orbital_index = configuration[electron_index]
        orbital_index = LL_k_map[LL_orbital_index]

        Slatermatrix[:,electron_index] .=
            Bloch_states[:,orbital_index]
    end

    return Slatermatrix
end






