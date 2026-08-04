using LinearAlgebra
using JLD2

BLAS.set_num_threads(1)

struct sigma_related
    pv::ComplexF64
    ev::ComplexF64
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

    return overlap,C1
end

function calculate_all_slater_pair_Gq(
    orbital_basis_orthogonal_list,
    precomp::GqPrecomp,
)
    Nconfig = length(orbital_basis_orthogonal_list)
    Nelectron = size(orbital_basis_orthogonal_list[1], 2)
    Nq = length(precomp.qkeys)

    # --------------------------------------------------------
    # Full pair data that we want to KEEP.
    # --------------------------------------------------------

    slater_overlap_matrix =
        zeros(ComplexF64, Nconfig, Nconfig)

    Gq_slater_pair =
        zeros(ComplexF64, Nconfig, Nconfig, Nq)
    
    first_minor_matrix =
    zeros(
        ComplexF64,
        Nelectron,
        Nelectron,
        Nconfig,
        Nconfig,
    )

    # --------------------------------------------------------
    # Memory estimate.
    # --------------------------------------------------------

    memory_Gq =
        sizeof(ComplexF64) * Nconfig^2 * Nq

    println(
        "Gq_slater_pair memory = ",
        round(memory_Gq / 1024^3, digits=3),
        " GiB",
    )

    # --------------------------------------------------------
    # Only explicitly calculate i <= j.
    #
    # Flattening this list gives much better load balance than
    # threading directly over i, because every pair has roughly
    # the same computational cost.
    # --------------------------------------------------------

    configuration_pairs = Tuple{Int,Int}[]
    sizehint!(
        configuration_pairs,
        Nconfig * (Nconfig + 1) ÷ 2,
    )

    for i in 1:Nconfig
        for j in i:Nconfig
            push!(configuration_pairs, (i,j))
        end
    end
        Npair = length(configuration_pairs)
        Nthread = Threads.nthreads()
        Nworkspace = Threads.maxthreadid()

        println("Nconfig       = ", Nconfig)
        println("Npair         = ", Npair)
        println("Nq            = ", Nq)
        println("Threads       = ", Nthread)
        println("Max thread ID = ", Nworkspace)

        workspaces = [
            SlaterPairWorkspace(Nelectron, Nq)
            for _ in 1:Nworkspace
        ]

    print_lock = ReentrantLock()

    Threads.@threads :static for ipair in eachindex(configuration_pairs)
        tid = Threads.threadid()
        workspace = workspaces[tid]

        i, j = configuration_pairs[ipair]

        Dmatrix = orbital_basis_orthogonal_list[i]
        Fmatrix = orbital_basis_orthogonal_list[j]

        overlap,C1 =
            slater_pair_Gq!(
                workspace,
                Dmatrix,
                Fmatrix,
                precomp,
            )

        # ----------------------------------------------------
        # Store i,j directly.
        # ----------------------------------------------------

        slater_overlap_matrix[i,j] = overlap
        @views first_minor_matrix[:,:,i,j] .= C1

        for iq in 1:Nq
            Gq_slater_pair[i,j,iq] = workspace.Gq[iq]
        end

        # ----------------------------------------------------
        # The reverse pair follows from Hermiticity:
        #
        # <D_j|G(q)|D_i>
        #     = <D_i|G(-q)|D_j>*
        # ----------------------------------------------------

        if i != j
            slater_overlap_matrix[j,i] = conj(overlap)
            @views first_minor_matrix[:,:,j,i] .= C1'

            for iq in 1:Nq
                iminus = precomp.minus_q_index[iq]

                Gq_slater_pair[j,i,iq] =
                    conj(workspace.Gq[iminus])
            end
        end

        if ipair % 100 == 0 || ipair == Npair
            lock(print_lock) do
                println(
                    "Slater pairs: $ipair / $Npair, " *
                    "thread $tid"
                )
                flush(stdout)
            end
        end
    end

    return (
    slater_overlap_matrix,
    Gq_slater_pair,
    first_minor_matrix)
end






function construct_manybody_Gq(
    Gq_slater_pair,
    overall_mag_matrix,
)
    Nconfig, Nconfig2, Nq =
        size(Gq_slater_pair)

    @assert Nconfig == Nconfig2
    @assert size(overall_mag_matrix, 1) == Nconfig

    Nvec = size(overall_mag_matrix, 2)

    Gq_raw_array =
        zeros(ComplexF64, Nvec, Nvec, Nq)

    Threads.@threads :static for iq in 1:Nq
        GDD =
            @view Gq_slater_pair[:,:,iq]

        Gq_raw_array[:,:,iq] .=
            overall_mag_matrix' *
            GDD *
            overall_mag_matrix
    end

    return Gq_raw_array
end

function metric_gram_schmidt(overlap)
    N = size(overlap, 1)

    @assert size(overlap, 2) == N
    @assert abs(det(overlap)) > 1e-8

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



function construct_manybody_rho(
    orbital_basis_orthogonal_list,
    first_minor_matrix,
    overall_mag_matrix,
)
    Nconfig = length(orbital_basis_orthogonal_list)
    M_basis, Nelectron = size(orbital_basis_orthogonal_list[1])
    Nvec = size(overall_mag_matrix, 2)

    U = hcat(orbital_basis_orthogonal_list...)
    Nlarge = Nconfig * Nelectron

    @assert size(U) == (M_basis, Nlarge)

    rho_raw = zeros(
        ComplexF64,
        M_basis,
        M_basis,
        Nvec,
        Nvec,
    )

    Threads.@threads :static for iab in 1:Nvec^2
        a = (iab - 1) ÷ Nvec + 1
        b = (iab - 1) % Nvec + 1

        B = zeros(ComplexF64, Nlarge, Nlarge)

        for i in 1:Nconfig
            col_range = (i-1)*Nelectron+1:i*Nelectron

            for j in 1:Nconfig
                row_range = (j-1)*Nelectron+1:j*Nelectron

                weight = conj(overall_mag_matrix[i,a]) *
                         overall_mag_matrix[j,b]

                C1 = @view first_minor_matrix[:,:,i,j]

                @views B[row_range,col_range] .=
                    weight .* transpose(C1)
            end
        end

        rho_ab = U * B * U'

        @views rho_raw[:,:,a,b] .= rho_ab
    end

    return rho_raw
end



function orthogonalize_Gq(
    Gq_raw_array,
    orthogonalize_matrix,
)
    Nvec, Nvec2, Nq = size(Gq_raw_array)

    @assert Nvec == Nvec2

    Gq_orth_array =
        zeros(ComplexF64, Nvec, Nvec, Nq)

    Threads.@threads :static for iq in 1:Nq
        Graw = @view Gq_raw_array[:,:,iq]

        Gq_orth_array[:,:,iq] .=
            orthogonalize_matrix' *
            Graw *
            orthogonalize_matrix
    end

    return Gq_orth_array
end


function orthogonalize_rho(
    rho_raw,
    orthogonalize_matrix,
)
    M_basis, M_basis2, Nvec, Nvec2 = size(rho_raw)

    @assert M_basis == M_basis2
    @assert Nvec == Nvec2

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
                orthogonalize_matrix' *
                rho_state *
                orthogonalize_matrix
        end
    end

    return rho_orth
end


function Gq_array_to_dict(Gq_array, qkeys)
    Nq = length(qkeys)

    @assert size(Gq_array, 3) == Nq

    Gq_dict =
        Dict{Tuple{Int,Int},Matrix{ComplexF64}}()

    for iq in 1:Nq
        Gq_dict[qkeys[iq]] =
            Matrix(@view Gq_array[:,:,iq])
    end

    return Gq_dict
end

function get_input_path(args)
    scratch_dir = ENV["SCRATCH"]
    file_pos = Int(args[14])

    filename =
        "$(args[1])f1$(args[2])f2$(args[3])NL$(args[4])am" *
        "$(args[5])N1$(args[6])N2$(args[7])N1f$(args[8])N2f" *
        "$(args[9])xir$(args[10])xii$(args[11])grid" *
        "$(args[12])type$(args[13])gaud.jld2"

    return joinpath(
        scratch_dir,
        "constrcut_twothirds_skv_v2/data_output$(file_pos)",
        filename,
    )
end


function main_func(args::Vector{Float64})
    filepath = get_input_path(args)
    isfile(filepath) || error("No file exists: $filepath")

    st = load(filepath)

    params = st["params"]
    topo_sectors = st["topo_sectors"]
    overall_mag_matrix = st["overall_mag_matrix"]
    spinor_set = st["spinor_set"]
    possible_config = st["possible_config"]
    orbital_basis_orthogonal_list =
        st["orbital_basis_orthogonal_list"]

    Nconfig = length(orbital_basis_orthogonal_list)
    Nvec = size(overall_mag_matrix, 2)

    @assert size(overall_mag_matrix, 1) == Nconfig
    @assert Nvec == 3

    M_basis, Nelectron =
        size(orbital_basis_orthogonal_list[1])

    for Dmatrix in orbital_basis_orthogonal_list
        @assert size(Dmatrix) == (M_basis, Nelectron)
    end

    # --------------------------------------------------------
    # Precompute all q-dependent single-particle information.
    # --------------------------------------------------------

    precomp =
        GqPrecomp(
            params.Tgrid,
            spinor_set,
            Nelectron,
        )

    Nq = length(precomp.qkeys)

    println("M_basis   = ", M_basis)
    println("Nelectron = ", Nelectron)
    println("Nconfig   = ", Nconfig)
    println("Nq        = ", Nq)
    println("Threads   = ", Threads.nthreads())

    # --------------------------------------------------------
    # Expensive step:
    #
    # slater_overlap_matrix[i,j] = <D_i|D_j>
    #
    # Gq_slater_pair[i,j,iq]
    #     = <D_i|G(q_iq)|D_j>
    #
    # These are retained and saved.
    # --------------------------------------------------------

    slater_overlap_matrix, Gq_slater_pair, first_minor_matrix =
        calculate_all_slater_pair_Gq(
            orbital_basis_orthogonal_list,
            precomp,
        )
    for i in 1:Nconfig
    @assert isapprox(
        @view(first_minor_matrix[:,:,i,i]),
        Matrix{ComplexF64}(I, Nelectron, Nelectron);
        atol = 1e-8,
        rtol = 1e-8,
    )
end

    # --------------------------------------------------------
    # Checks directly in the Slater-determinant basis.
    # --------------------------------------------------------

    @assert isapprox(
        slater_overlap_matrix,
        slater_overlap_matrix';
        atol = 1e-9,
        rtol = 1e-9,
    )

    iq0 = precomp.qindex[(0,0)]

    @assert isapprox(
        @view(Gq_slater_pair[:,:,iq0]),
        Nelectron * (Nelectron - 1) *
        slater_overlap_matrix;
        atol = 1e-7,
        rtol = 1e-7,
    )

    for iq in 1:Nq
        iminus = precomp.minus_q_index[iq]

        @assert isapprox(
            @view(Gq_slater_pair[:,:,iq])',
            @view(Gq_slater_pair[:,:,iminus]);
            atol = 1e-7,
            rtol = 1e-7,
        )
    end

    # --------------------------------------------------------
    # Raw three-state many-body overlap.
    #
    # |Psi_s> = sum_i C[i,s] |D_i>
    #
    # Therefore:
    #
    # deno_raw = C' * D_slater * C
    # --------------------------------------------------------

    deno_raw =
        overall_mag_matrix' *
        slater_overlap_matrix *
        overall_mag_matrix

    @assert isapprox(
        deno_raw,
        deno_raw';
        atol = 1e-9,
        rtol = 1e-9,
    )

    # Remove only floating-point non-Hermiticity.
    deno_raw .= (deno_raw + deno_raw') / 2

    println("Raw many-body overlap matrix:")
    display(deno_raw)

    # --------------------------------------------------------
    # Raw three-state G(q).
    #
    # G_raw(q) = C' * G_DD(q) * C
    # --------------------------------------------------------

    Gq_raw_array =
        construct_manybody_Gq(
            Gq_slater_pair,
            overall_mag_matrix,
        )
    println("Constructing rho_raw...")

    rho_raw = construct_manybody_rho(
        orbital_basis_orthogonal_list,
        first_minor_matrix,
        overall_mag_matrix,
    )
    # q = 0 must satisfy
    #
    # G_raw(0) = Ne(Ne-1) deno_raw.
    @assert isapprox(
        @view(Gq_raw_array[:,:,iq0]),
        Nelectron * (Nelectron - 1) * deno_raw;
        atol = 1e-7,
        rtol = 1e-7,
    )

    for iq in 1:Nq
        iminus = precomp.minus_q_index[iq]

        @assert isapprox(
            @view(Gq_raw_array[:,:,iq])',
            @view(Gq_raw_array[:,:,iminus]);
            atol = 1e-7,
            rtol = 1e-7,
        )
    end

    for a in 1:Nvec, b in 1:Nvec
    @assert isapprox(
        tr(@view rho_raw[:,:,a,b]),
        Nelectron * deno_raw[a,b];
        atol = 1e-7,
        rtol = 1e-7,
    )

    @assert isapprox(
        @view(rho_raw[:,:,a,b])',
        @view(rho_raw[:,:,b,a]);
        atol = 1e-8,
        rtol = 1e-8,
    )
    end

    # --------------------------------------------------------
    # Orthonormalize the three many-body states.
    # --------------------------------------------------------

    orthogonalize_matrix =
        metric_gram_schmidt(deno_raw)

    deno_orth =
        orthogonalize_matrix' *
        deno_raw *
        orthogonalize_matrix

    rho_orth = orthogonalize_rho(
    rho_raw,
    orthogonalize_matrix,
    )

    @assert isapprox(
        deno_orth,
        Matrix{ComplexF64}(I, Nvec, Nvec);
        atol = 1e-10,
        rtol = 1e-10,
    )

    for a in 1:Nvec, b in 1:Nvec
        expected = a == b ? Nelectron : 0.0

        @assert isapprox(
            tr(@view rho_orth[:,:,a,b]),
            expected;
            atol = 1e-7,
            rtol = 1e-7,
        )
    end

    println("Orthogonalization matrix:")
    display(orthogonalize_matrix)

    println("Overlap after orthogonalization:")
    display(deno_orth)

    # --------------------------------------------------------
    # Transform G(q) into the orthonormal manifold.
    # --------------------------------------------------------

    Gq_orth_array =
        orthogonalize_Gq(
            Gq_raw_array,
            orthogonalize_matrix,
        )

    # At q=0:
    #
    # G_orth(0) = Ne(Ne-1) I.
    @assert isapprox(
        @view(Gq_orth_array[:,:,iq0]),
        Nelectron * (Nelectron - 1) *
        Matrix{ComplexF64}(I, Nvec, Nvec);
        atol = 1e-7,
        rtol = 1e-7,
    )

    for iq in 1:Nq
        iminus = precomp.minus_q_index[iq]

        @assert isapprox(
            @view(Gq_orth_array[:,:,iq])',
            @view(Gq_orth_array[:,:,iminus]);
            atol = 1e-7,
            rtol = 1e-7,
        )
    end

    # --------------------------------------------------------
    # Dict versions for compatibility with old post-processing.
    # --------------------------------------------------------

    Gq_raw =
        Gq_array_to_dict(
            Gq_raw_array,
            precomp.qkeys,
        )

    Gq_orth =
        Gq_array_to_dict(
            Gq_orth_array,
            precomp.qkeys,
        )

    # --------------------------------------------------------
    # Return a NAMED TUPLE.
    #
    # I prefer this over returning 15 positional variables:
    # it prevents accidentally swapping objects in main_test.
    # --------------------------------------------------------

    return (
        params = params,
        topo_sectors = topo_sectors,
        possible_config = possible_config,
        overall_mag_matrix = overall_mag_matrix,
        orbital_basis_orthogonal_list =orbital_basis_orthogonal_list,
        spinor_set = spinor_set,

        qkeys = precomp.qkeys,
        qindex = precomp.qindex,
        minus_q_index = precomp.minus_q_index,

        slater_overlap_matrix = slater_overlap_matrix,
        Gq_slater_pair = Gq_slater_pair,
        first_minor_matrix = first_minor_matrix,
        rho_raw = rho_raw,

        deno_raw = deno_raw,
        Gq_raw_array = Gq_raw_array,
        Gq_raw = Gq_raw,

        orthogonalize_matrix = orthogonalize_matrix,

        deno_orth = deno_orth,
        rho_orth = rho_orth,
        Gq_orth_array = Gq_orth_array,
        Gq_orth = Gq_orth,
    )
end