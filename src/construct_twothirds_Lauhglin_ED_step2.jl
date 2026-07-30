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
    LL_filepath=joinpath(scratch_dir, "constrcut_twothirds_Lauhglin_ED/data_output$(Int(args[9]))/$(args[1])f1$(args[2])f2$(args[3])q1$(args[4])q2$(args[5])N1$(args[6])N2$(args[7])Npa$(args[8])Nvec.jld2")


    LL_data=load(LL_filepath)
    LL_allowedq=LL_data["allowedq"]
    LL_T1=LL_data["T1"]
    LL_T2=LL_data["T2"]

    @assert norm(LL_T1/moiream-T1_vec)<10^(-9)
    @assert norm(LL_T2/moiream-T2_vec)<10^(-9)
    

    PH_eigvector=LL_data["PH_eigvector"]
    PH_state_can=LL_data["PH_state_can"]
    @assert length(LL_data["PH_eigvector"])==Nvec

    LL_k_set=[wrapparallel(complex_vec([LL_T1 LL_T2]*LL_allowedq[ja])/moiream,b1,b2) for ja in eachindex(LL_allowedq)]
    k_set=vec([wrapparallel(n1*T1+n2*T2,b1,b2) for n1 in 1:N1, n2 in 1:N2])

    LL_k_map=zeros(Int,length(k_set))
    for ja in eachindex(LL_k_map)
        ccount=0
        for jb in eachindex(k_set)
        if abs(k_set[jb]-LL_k_set[ja])<10^(-9)
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

    for ja in 1:N1*N2
    Bloch_states[:,ja]=Bloch_states[:,ja]/norm(Bloch_states[:,ja])
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


    
    orbital_basis_record=[Vector{Matrix{ComplexF64}}(undef,length(PH_eigvector[s])) for s in 1:Nvec] #This is just to construct it for the orthogonalized ones.
    coeff_record=[Vector{ComplexF64}(undef,length(PH_eigvector[s])) for s in 1:Nvec]


    for jveci in 1:Nvec 
        for ja in eachindex(PH_eigvector[jveci])
            pdd=1
            Slatermatrix=zeros(ComplexF64,length(Tgrid),Nelectron)
            for jb in eachindex(PH_state_can[jveci][ja])
                pdd*=Nk[LL_k_map[PH_state_can[jveci][ja][jb]]]
                Slatermatrix[:,jb]=Bloch_states[:,LL_k_map[PH_state_can[jveci][ja][jb]]]
                
            end
            coeff_record[jveci][ja]=PH_eigvector[jveci][ja]*pdd
            orbital_basis_record[jveci][ja]=Slatermatrix
        end
    end




    # The next step will be to orthogonalize it.


        # 1. Divide states into blocks
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


        # 2. Orthonormalize each block
        #
        # |Psi_orth[s]> = sum_a |Psi_raw[a]> * reconstruct_matrix[a,s]
                reconstruct_matrix = zeros(ComplexF64, Nvec, Nvec)
            coeff_orthogonal = Vector{Vector{ComplexF64}}(undef, Nvec)
            orbital_basis_orthogonal =
                Vector{Vector{Matrix{ComplexF64}}}(undef, Nvec)

            overlap_manybody = zeros(ComplexF64, Nvec, Nvec)

            for block in blocks
                Cblock = hcat(coeff_record[block]...)
                overlap_manybody[block, block] .= Cblock' * Cblock
            end

            overlap_manybody =
                (overlap_manybody + overlap_manybody') / 2

            for block in blocks
                Cblock = hcat(coeff_record[block]...)

                F = qr(Cblock)

                nblock = length(block)
                Qblock = Matrix(F.Q)[:, 1:nblock]
                Rblock = Matrix(F.R)[1:nblock, 1:nblock]

                Ablock =
                    Rblock \ Matrix{ComplexF64}(I, nblock, nblock)

                reconstruct_matrix[block, block] .= Ablock

                for (local_index, global_index) in enumerate(block)
                    coeff_orthogonal[global_index] =
                        copy(Qblock[:, local_index])

                    # Same PH_state_can means the same ordered Slater basis.
                    orbital_basis_orthogonal[global_index] =
                        orbital_basis_record[block[1]]
                end
            end

                @assert isapprox(
            reconstruct_matrix' *
            overlap_manybody *
            reconstruct_matrix,
            Matrix{ComplexF64}(I, Nvec, Nvec);
            atol = 1e-10,
            rtol = 1e-10,
        )



   
















        spinor_set=Vector{ComplexF64}[]
        for ja in eachindex(Tgrid)
            tvec=Tgrid[ja][1]*T1+Tgrid[ja][2]*T2-(Deltatheta)
            vv=get_spinor(NL,[real(tvec),imag(tvec)])
            push!(spinor_set,vv)

        end

     





        precomp = QResolvedPrecomp(Tgrid, spinor_set)

        rho_matrix_all =
            Vector{Matrix{ComplexF64}}(undef, Nvec)

        final_G_q_all =
            Vector{Dict{Tuple{Int,Int},ComplexF64}}(undef, Nvec)

        deno_all = Vector{ComplexF64}(undef, Nvec)

        for which_state in 1:Nvec
            coeff_list =
                coeff_orthogonal[which_state]

            orbital_basis_orthogonal_list =
                orbital_basis_orthogonal[which_state]

            # Your existing deno, G(q), and rho calculation goes here.
                                                deno=0.0+0.0*im


                        for ja in eachindex(coeff_list)
                            for  jb in eachindex(coeff_list)
                                coeff=coeff_list[ja]'*coeff_list[jb]
                                Dmatrix=orbital_basis_orthogonal_list[ja]
                                Fmatrix=orbital_basis_orthogonal_list[jb]
                                deno+=det(Dmatrix'*Fmatrix)*coeff


                            end
                        end

                    


                            


                        final_G_q=Dict{Tuple{Int, Int}, ComplexF64}()

                        local_Gq_list=Matrix{Dict{Tuple{Int,Int},ComplexF64}}(undef,length(coeff_list),length(coeff_list))
                        print_lock = ReentrantLock()
                        Threads.@threads :greedy for ja in eachindex(coeff_list)

                                lock(print_lock) do
                                            println("starting $(ja)/$(length(coeff_list)) in Threads$(Threads.threadid()),")
                                            flush(stdout)
                                        end
                            for  jb in eachindex(coeff_list)
                        
                                Dmatrix=orbital_basis_orthogonal_list[ja]
                                Fmatrix=orbital_basis_orthogonal_list[jb]
                                local_Gq_list[ja,jb],_=q_resolved_pair_weight(
                                                                        Dmatrix,
                                                                        Fmatrix,
                                                                        precomp;
                                                                        skip_q0 = false,
                                                                    )


                            end
                        end

                    for ja in eachindex(coeff_list), jb in eachindex(coeff_list)
                        coeff=coeff_list[ja]'*coeff_list[jb]
                        for (qkey,vals) in local_Gq_list[ja,jb]
                                final_G_q[qkey] = (get(final_G_q, qkey, 0.0 + 0.0im) + vals*coeff/deno)
                        end

                    end


                    # ------------------------------------------------------------
                    # One-body density matrix:
                    #
                    # rho_matrix[a,b] = <Psi| c_b^dagger c_a |Psi> / <Psi|Psi>
                    # ------------------------------------------------------------

                    M_basis = size(orbital_basis_orthogonal_list[1], 1)
                    rho_num = zeros(ComplexF64, M_basis, M_basis)
                
                    for ja in eachindex(coeff_list), jb in eachindex(coeff_list)
                        coeff = coeff_list[ja]' * coeff_list[jb]

                        Dmatrix = orbital_basis_orthogonal_list[ja]
                        Fmatrix = orbital_basis_orthogonal_list[jb]

                        rho_ij, _ = onebody_transition_rho(Dmatrix, Fmatrix)

                        rho_num .+= coeff .* rho_ij
                    
                    end

                    rho_matrix = rho_num ./ deno


            rho_matrix_all[which_state] = rho_matrix
            final_G_q_all[which_state] = final_G_q
            deno_all[which_state] = deno
        end


       



    params=(a1=a1,
                a2=a2,
                lb=lb,
                b1=b1,
                b2=b2,
                α=α,

                L1=L1,
                L2=L2,
                Lb=Lb,
                T1=T1,
                T2=T2,
            

     

                N1=N1,
                N2=N2,

                Nphi=Nphi,
                Nelectron=Nelectron,

                Deltatheta=Deltatheta,
                qq=qq,
               
                T1_vec=T1_vec,
                T2_vec=T2_vec,
                b1_vec=b1_vec,
                b2_vec=b2_vec,
             
                b1T=b1T,
                b2T=b2T,
          

                a1_vec=a1_vec,
                a2_vec=a2_vec,
                L1_vec=L1_vec,
                L2_vec=L2_vec,

                Tgrid=Tgrid,
                Tgrid_dict=Tgrid_dict,
                ggrid=ggrid,
                ggrid_dict=ggrid_dict,
              
                grid_cutoff=grid_cutoff,

                k_set=k_set,
                k_set_n=k_set_n,

          
                NL=NL,
                spinor_set=spinor_set,
                LL_k_set=LL_k_set,
                LL_k_map=LL_k_map,
                PH_eigvector= PH_eigvector,
                PH_state_can=PH_state_can
               
        )







    return (
        rho_matrix_all,
        final_G_q_all,
        deno_all,
        reconstruct_matrix,
        overlap_manybody,
        coeff_orthogonal,
        orbital_basis_orthogonal,
        params,
    )




end







