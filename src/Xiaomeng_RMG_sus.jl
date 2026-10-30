using LinearAlgebra


function get_Ham(k::Vector{Float64}, uD::Float64,
                 valley::Int64, stacking::Int, NL::Int)
    @assert valley in (-1, 1)
    @assert stacking in (-1, 1)

    a = 0.246  # nm
    t0, t1, t2, t3, t4 = 3100.0, 380.0, -15.0, 290.0, 141.0

    # k is relative to the valley center, in nm^-1.
    f = -sqrt(3)*a/2 * (valley*k[1] - im*stacking*k[2])

    Ham = zeros(ComplexF64, 2*NL, 2*NL)

    for l in 1:(NL-1)
        Ham[2*l-1:2*l, 2*l+1:2*l+2] =
            [t4*f t3*conj(f); t1 t4*f]
    end

    for l in 1:(NL-2)
        Ham[2*l-1, 2*l+4] = t2/2
    end

    Ham = Ham + Ham'

    for l in 1:NL
        U = uD*(l-(NL+1)/2)
        Ham[2*l-1:2*l, 2*l-1:2*l] =
            [U -t0*f; -t0*conj(f) U]
    end

    return Ham
end


function get_HamQ(kpoints, Q, D, U, UQ, FQ, w, εr, dgate)
    N = length(kpoints)
    H = Matrix{ComplexF64}(undef, N, N)
    VQ = Vscreened(norm(Q), εr, dgate)

    @views for b in 1:N, a in 1:N
        q = hypot(kpoints[a][1] - kpoints[b][1],
                  kpoints[a][2] - kpoints[b][2])

        hartree = VQ * FQ[a] * conj(FQ[b])

        exchange = Vscreened(q, εr, dgate) *
                   dot(U[:, b], U[:, a]) *
                   dot(UQ[:, a], UQ[:, b])

        H[a, b] = w * (hartree - exchange)
    end

    for a in 1:N
        H[a, a] += D[a]
    end

    return H
end


# q in nm^-1; returned interaction in meV·nm²
function Vscreened(q, εr, dgate)
    return (9047.5636 / εr) * (iszero(q) ? dgate : tanh(q*dgate)/q)
end



function get_mu(;
    uD, n_h,
    NL=6, valley=1, stacking=1, band=NL,
    kmax=1.0, Nk=501
)
    @assert kmax > 0 && Nk >= 2
    @assert n_h > 0
    @assert 1 <= band <= 2NL

    # Uniform Cartesian grid for determining μ
    dk = 2kmax / Nk
    kgrid = range(-kmax + dk/2, kmax - dk/2; length=Nk)
    w = dk^2 / (2π)^2

    # Rows correspond to ky; columns correspond to kx.
    εgrid = [
        eigvals(Hermitian(get_Ham([kx, ky], uD, valley, stacking, NL)))[band]
        for ky in kgrid, kx in kgrid
    ]

    # Holes occupy the highest-energy states of this band.
    εsorted = sort(vec(εgrid); rev=true)
    Nholes = round(Int, n_h / w)

    @assert 1 <= Nholes < length(εsorted) "Density is outside the grid's range."

    # Place μ between the last hole state and first occupied electron state.
    μ = (εsorted[Nholes] + εsorted[Nholes + 1]) / 2

    fgrid = εgrid .< μ
    n_h_grid = w * count(.!fgrid)

    # Check that μ produces the selected discrete hole count.
    Nholes = count(.!fgrid)  # Record the actual count.

    # Check that the sampled hole region stays inside the box.
    @assert all(fgrid[1, :]) && all(fgrid[end, :]) &&
            all(fgrid[:, 1]) && all(fgrid[:, end]) "Increase kmax: holes reach the grid boundary."

    return (; μ, n_h_grid, kgrid, εgrid, fgrid, dk, w, Nholes)
end



function get_pair_grid(
    Q_index::AbstractVector{<:Integer},
    μ, uD, valley, stacking, NL, band,
    Rpair, dk_pair, gridshape::Int
)
    @assert length(Q_index) == 2
    @assert Rpair > 0 && dk_pair > 0
    @assert 1 <= band <= 2NL

    radius2 = (Rpair / dk_pair)^2

    if gridshape == 1
        b1 = [1.0, 0.0]
        b2 = [0.5, sqrt(3)/2]

        N = 2*ceil(Int, 2Rpair / (sqrt(3)*dk_pair))

        # |i*b1 + j*b2|² = i² + i*j + j²
        ij = [
            (i, j)
            for j in -N:N for i in -N:N
            if i^2 + i*j + j^2 < radius2
        ]

    elseif gridshape == 2
        b1 = [1.0, 0.0]
        b2 = [0.0, 1.0]

        N = 2*ceil(Int, Rpair / dk_pair)

        # |i*b1 + j*b2|² = i² + j²
        ij = [
            (i, j)
            for j in -N:N for i in -N:N
            if i^2 + j^2 < radius2
        ]

    else
        error("Unknown gridshape: $gridshape")
    end

    # Actual primitive vectors and physical Q, in nm^-1.
    B = dk_pair * hcat(b1, b2)
    Q = B * Q_index

    @assert norm(Q) > 0 "This calculation requires a nonzero Q."

    w_pair = abs(det(B)) / (2π)^2

    # Origin-centered lattice with circular cutoff |k| < Rpair.
    kmesh = [B * [i, j] for (i, j) in ij]

    # No additional cutoff on k + Q.
    kQmesh = [k + Q for k in kmesh]

    # Band energies at both endpoints
    εpair = [
        eigvals(Hermitian(get_Ham(k, uD, valley, stacking, NL)))[band]
        for k in kmesh
    ]

    εpairQ = [
        eigvals(Hermitian(get_Ham(kQ, uD, valley, stacking, NL)))[band]
        for kQ in kQmesh
    ]

    # Electron occupations and active pairs
    fpair = εpair .< μ
    fpairQ = εpairQ .< μ
    active = fpair .!= fpairQ

    Nactive = count(active)
    @assert Nactive > 0 "No active pairs; check μ or refine dk_pair."

    # Check the outermost shell, one grid spacing thick.
    boundary = [norm(k) >= Rpair - dk_pair for k in kmesh]
    @assert !any(active .& boundary) "Increase Rpair: active pairs approach the cutoff."

    # Collect active pairs.
    k_active = kmesh[active]
    kQ_active = kQmesh[active]
    D = abs.(εpairQ[active] .- εpair[active])

    # Active-pair index -> full-grid index
    active_index = findall(active)

    # Full-grid index -> active-pair index; zero means inactive.
    index_grid = zeros(Int, length(kmesh))
    index_grid[active] = 1:Nactive

    return (;
        Q_index, Q, gridshape, dk_pair, b1, b2,
        ij, w_pair,
        kmesh, kQmesh,
        εpair, εpairQ, fpair, fpairQ,
        active, active_index, index_grid,
        k_active, kQ_active, D, Nactive,
    )
end


function get_pair_spinors(
    k_active, kQ_active,
    uD, valley, stacking, NL, band, modeltype
)
    @assert modeltype in (:spinor, :no_spinor)
    @assert 1 <= band <= 2NL

    Nactive = length(k_active)
    @assert length(kQ_active) == Nactive
    @assert Nactive > 0

    if modeltype == :spinor
        U = Matrix{ComplexF64}(undef, 2NL, Nactive)
        UQ = Matrix{ComplexF64}(undef, 2NL, Nactive)

        for a in 1:Nactive
            eig_k = eigen(Hermitian(
                get_Ham(k_active[a], uD, valley, stacking, NL)
            ))

            eig_kQ = eigen(Hermitian(
                get_Ham(kQ_active[a], uD, valley, stacking, NL)
            ))

            U[:, a] = eig_k.vectors[:, band]
            UQ[:, a] = eig_kQ.vectors[:, band]
        end

    else
        # No-spinor model: replace each spinor by the scalar 1.
        U = ones(ComplexF64, 1, Nactive)
        UQ = ones(ComplexF64, 1, Nactive)
    end

    # FQ[a] = ⟨u_(k_a+Q) | u_(k_a)⟩
    FQ = [
        dot(UQ[:, a], U[:, a])
        for a in 1:Nactive
    ]

    return (; U, UQ, FQ)
end


function get_lowest_modes(
    HQ, FQ, w_pair, nmodes;
    hermiticity_tol=1e-10,
    residual_tol=1e-8
)
    N = size(HQ, 1)

    @assert N > 0 && size(HQ, 2) == N
    @assert length(FQ) == N
    @assert nmodes isa Integer && nmodes > 0
    @assert w_pair > 0

    # Check Hermiticity without allocating another full matrix.
    hermiticity_error = maximum(
        abs(HQ[a, b] - conj(HQ[b, a]))
        for b in 1:N for a in 1:b
    )

    @assert hermiticity_error < hermiticity_tol "HQ is not Hermitian within tolerance."

    # Lowest eigenvalues in ascending order
    nsave = min(nmodes, N)
    eigQ = eigen(Hermitian(HQ), 1:nsave)

    eigenvalues = eigQ.values
    modes = eigQ.vectors / sqrt(w_pair)

    normalizations = zeros(Float64, nsave)
    residuals = zeros(Float64, nsave)
    dQ = zeros(ComplexF64, nsave)

    for n in 1:nsave
        X = @view modes[:, n]
        κ = eigenvalues[n]

        normalizations[n] = w_pair * sum(abs2, X)

        residuals[n] = norm(HQ*X - κ*X) / norm(X)

        # FQ[k] = ⟨u_(k+Q) | u_k⟩
        dQ[n] = w_pair * dot(FQ, X)
    end

    @assert maximum(residuals) < residual_tol "An eigenmode residual exceeds the tolerance."

    return (;
        eigenvalues, modes,
        normalizations, residuals, dQ,
        hermiticity_error,
    )
end