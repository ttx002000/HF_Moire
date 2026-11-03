using LinearAlgebra
using FastGaussQuadrature

using Random

using KrylovKit

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



function pair_energies(r, θ, Q, uD, valley, stacking, NL, band)
    # Midpoint on the ray at angle θ.
    p = r * [cos(θ), sin(θ)]

    # The two physical momenta belonging to this pair.
    kminus = p - Q/2
    kplus  = p + Q/2

    εminus = eigvals(Hermitian(get_Ham(kminus, uD, valley, stacking, NL)))[band]
    εplus  = eigvals(Hermitian(get_Ham(kplus, uD, valley, stacking, NL)))[band]

    return εminus, εplus
end



function ray_brackets(θ, Q, μ, uD, valley, stacking, NL, band;
                      rmax=0.8, nscan=800)
    @assert rmax > 0 && nscan >= 1

    rgrid = range(0.0, rmax; length=nscan+1)
    brackets = Tuple{Float64,Float64,Int}[]

    εleft = pair_energies(rgrid[1], θ, Q, uD, valley, stacking, NL, band)

    for i in 1:nscan
        εright = pair_energies(rgrid[i+1], θ, Q, uD, valley, stacking, NL, band)

        for s in 1:2
            if (εleft[s] - μ) * (εright[s] - μ) <= 0
                push!(brackets, (rgrid[i], rgrid[i+1], s))
            end
        end

        εleft = εright
    end

    return brackets
end


function refine_crossing(bracket, θ, Q, μ, uD, valley, stacking, NL, band;
                         r_tol=1e-10)
    @assert r_tol > 0
    a, b, s = bracket

    fa = pair_energies(a, θ, Q, uD, valley, stacking, NL, band)[s] - μ
    fb = pair_energies(b, θ, Q, uD, valley, stacking, NL, band)[s] - μ

    fa == 0 && return a
    fb == 0 && return b
    @assert fa * fb < 0

    while b - a > r_tol
        c = (a + b)/2
        (c == a || c == b) && return c
        fc = pair_energies(c, θ, Q, uD, valley, stacking, NL, band)[s] - μ

        fc == 0 && return c

        if fa * fc < 0
            b = c
        else
            a = c
            fa = fc
        end
    end

    return (a + b)/2
end


function ray_intervals(θ, Q, μ, uD, valley, stacking, NL, band;
                       rmax=0.8, nscan=400, r_tol=1e-10)
    brackets = ray_brackets(θ, Q, μ, uD, valley, stacking, NL, band;
                            rmax=rmax, nscan=nscan)

    edges = [0.0, rmax]
    for bracket in brackets
        r = refine_crossing(bracket, θ, Q, μ, uD, valley, stacking, NL, band;
                            r_tol=r_tol)
        push!(edges, r)
    end
    sort!(unique!(edges))

    intervals = Tuple{Float64,Float64}[]
    for i in 1:(length(edges)-1)
        a, b = edges[i], edges[i+1]
        r = (a + b)/2
        εminus, εplus = pair_energies(r, θ, Q, uD, valley, stacking, NL, band)

        if (εminus < μ) != (εplus < μ)
            push!(intervals, (a, b))
        end
    end

    return intervals
end


function radial_quadrature(intervals, Nr)
    @assert Nr >= 1
    x, w = gausslegendre(Nr)

    radii = Float64[]
    weights = Float64[]

    for (a, b) in intervals
        # Map the standard nodes from [-1, 1] to [a, b].
        r = (a + b)/2 .+ (b - a)/2 .* x

        # Radial integration uses r dr in polar coordinates.
        wr = (b - a)/2 .* w .* r

        append!(radii, r)
        append!(weights, wr)
    end

    return radii, weights
end


function polar_grid(Q, μ, uD, valley, stacking, NL, band;
                    Nθ=96, Nr=16, rmax=0.8, nscan=400, r_tol=1e-10)

    @assert length(Q) == 2 && norm(Q) > 0
    @assert Nθ >= 6 && Nθ % 6 == 0 "Nθ must be a positive multiple of 6."

    # Look for a physical mirror that preserves or reverses Q.
    M = nothing
    α = 0.0

    for axis in 0:2
        αtest = axis * π/3
        Mtest = [cos(2αtest)  sin(2αtest);
                 sin(2αtest) -cos(2αtest)]

        if isapprox(Mtest * Q, -Q; atol=0, rtol=1e-10) ||
           isapprox(Mtest * Q,  Q; atol=0, rtol=1e-10)
            M = Mtest
            α = αtest
            break
        end
    end

    mirror_found = M !== nothing

    p_points = Vector{Float64}[]
    weights = Float64[]
    Δθ = 2π / Nθ

    last_ray = mirror_found ? Nθ÷2 : Nθ-1

    for j in 0:last_ray
        θ = α + j * Δθ
        direction = [cos(θ), sin(θ)]

        # Put boundary rays along the selected mirror axis.
        if mirror_found && j == 0
            direction = [cos(α), sin(α)]
        elseif mirror_found && j == Nθ÷2
            direction = -[cos(α), sin(α)]
        end

        intervals = ray_intervals(
            θ, Q, μ, uD, valley, stacking, NL, band;
            rmax=rmax, nscan=nscan, r_tol=r_tol
        )

        radii, radial_weights = radial_quadrature(intervals, Nr)

        for (r, wr) in zip(radii, radial_weights)
            p = r * direction
            W = Δθ * wr / (2π)^2

            push!(p_points, p)
            push!(weights, W)

            # Reflect interior rays; boundary rays appear only once.
            if mirror_found && 0 < j < Nθ÷2
                push!(p_points, M * p)
                push!(weights, W)
            end
        end
    end

    return p_points, weights
end


function pair_data(p_points, Q, uD, valley, stacking, NL, band)
    N = length(p_points)

    D = zeros(Float64, N)
    Uminus = zeros(ComplexF64, 2*NL, N)
    Uplus  = zeros(ComplexF64, 2*NL, N)
    FQ = zeros(ComplexF64, N)

    for i in 1:N
        kminus = p_points[i] - Q/2
        kplus  = p_points[i] + Q/2

        eminus = eigen(Hermitian(get_Ham(kminus, uD, valley, stacking, NL)))
        eplus  = eigen(Hermitian(get_Ham(kplus, uD, valley, stacking, NL)))

        D[i] = abs(eplus.values[band] - eminus.values[band])

        Uminus[:, i] = eminus.vectors[:, band]
        Uplus[:, i]  = eplus.vectors[:, band]

        FQ[i] = dot(Uplus[:, i], Uminus[:, i])
    end

    return D, Uminus, Uplus, FQ
end


function Vscreened(q, εr, interaction_type::Integer)
    if interaction_type == 1
        # Gate screening: fixed gate distance in nm.
        dgate = 30.0

        return (9047.5636 / εr) *
               (iszero(q) ? dgate : tanh(q*dgate)/q)

    elseif interaction_type == 2
        # Thomas–Fermi screening: provisional qTF in nm^-1.
        qTF = 0.1

        return (9047.5636 / εr) / (q + qTF)

    else
        error("Unknown interaction_type: $interaction_type")
    end
end



function dense_hessian(p_points, Q, weights, D,
                       Uminus, Uplus, FQ, εr, interaction_type)
    N = length(p_points)
    H = zeros(ComplexF64, N, N)

    sqrtw = sqrt.(weights)
    VQ = Vscreened(norm(Q), εr, interaction_type)

    for b in 1:N, a in 1:N
        q = norm(p_points[a] - p_points[b])

        hartree = VQ * FQ[a] * conj(FQ[b])

        exchange = Vscreened(q, εr, interaction_type) *
                   dot(Uminus[:, b], Uminus[:, a]) *
                   dot(Uplus[:, a], Uplus[:, b])

        H[a, b] = sqrtw[a] * sqrtw[b] * (hartree - exchange)
    end

    for a in 1:N
        H[a, a] += D[a]
    end

    return H
end

function apply_hessian(x, p_points, Q, weights, D,
                       Uminus, Uplus, FQ, εr, interaction_type)
    N = length(p_points)
    sqrtw = sqrt.(weights)
    g = sqrtw .* FQ

    # Kinetic term and rank-one Hartree term
    y = D .* x + Vscreened(norm(Q), εr, interaction_type) * g * dot(g, x)

    # Exchange term
    for b in 1:N, a in 1:N
        q = norm(p_points[a] - p_points[b])

        exchange = Vscreened(q, εr, interaction_type) *
                   dot(Uminus[:, b], Uminus[:, a]) *
                   dot(Uplus[:, a], Uplus[:, b])

        y[a] -= sqrtw[a] * sqrtw[b] * exchange * x[b]
    end

    return y
end


struct MomentumCluster
    positions::UnitRange{Int}
    lower::Vector{Float64}
    upper::Vector{Float64}
    children::Vector{MomentumCluster}
end

#=
function momentum_tree(p_points; leafsize=64)
    @assert leafsize >= 1 && !isempty(p_points)

    order = collect(eachindex(p_points))

    function split_cluster!(positions)
        lower = [
            minimum(p_points[order[i]][d] for i in positions)
            for d in 1:2
        ]

        upper = [
            maximum(p_points[order[i]][d] for i in positions)
            for d in 1:2
        ]

        children = MomentumCluster[]

        if length(positions) > leafsize
            # Split along the direction with the larger spread.
            axis = argmax(upper - lower)

            sort!(
                view(order, positions);
                by=i -> p_points[i][axis]
            )

            middle = first(positions) + length(positions)÷2 - 1

            push!(children, split_cluster!(first(positions):middle))
            push!(children, split_cluster!((middle+1):last(positions)))
        end

        return MomentumCluster(positions, lower, upper, children)
    end

    tree = split_cluster!(1:length(order))
    return order, tree
end
=#

function momentum_tree(p_points; leafsize=64)
    @assert leafsize >= 1 && !isempty(p_points)

    function build_cluster(indices, first_position)
        n = length(indices)

        # Coordinates of the points in this cluster.
        x = [p_points[i][1] for i in indices]
        y = [p_points[i][2] for i in indices]

        lower = [minimum(x), minimum(y)]
        upper = [maximum(x), maximum(y)]

        # Where this cluster will appear in the final ordering.
        last_position = first_position + n - 1
        positions = first_position:last_position

        # Small enough: keep these points as one leaf.
        if n <= leafsize
            cluster = MomentumCluster(
                positions, lower, upper, MomentumCluster[]
            )
            return indices, cluster
        end

        # Sort along the direction with the larger spread.
        width_x = upper[1] - lower[1]
        width_y = upper[2] - lower[2]

        if width_x >= width_y
            permutation = sortperm(x)
        else
            permutation = sortperm(y)
        end

        sorted_indices = indices[permutation]

        # Divide into two groups of nearly equal size.
        middle = div(n, 2)
        left_indices = sorted_indices[1:middle]
        right_indices = sorted_indices[middle+1:end]

        # Build each child independently.
        left_order, left_cluster = build_cluster(
            left_indices, first_position
        )

        right_order, right_cluster = build_cluster(
            right_indices, first_position + middle
        )

        # Combine the children's results.
        order = vcat(left_order, right_order)
        children = [left_cluster, right_cluster]

        cluster = MomentumCluster(
            positions, lower, upper, children
        )

        return order, cluster
    end

    indices = collect(1:length(p_points))
    return build_cluster(indices, 1)
end











function well_separated(A::MomentumCluster, B::MomentumCluster; η=1.0)
    @assert η > 0

    gap_x = max(
        0.0,
        A.lower[1] - B.upper[1],
        B.lower[1] - A.upper[1]
    )

    gap_y = max(
        0.0,
        A.lower[2] - B.upper[2],
        B.lower[2] - A.upper[2]
    )

    distance = hypot(gap_x, gap_y)

    size_A = norm(A.upper - A.lower)
    size_B = norm(B.upper - B.lower)

    return distance > 0 && max(size_A, size_B) <= η * distance
end


function exchange_block(I, J, p_points, weights,
                        Uminus, Uplus, εr, interaction_type)
    B = zeros(ComplexF64, length(I), length(J))

    sqrtw_I = sqrt.(weights[I])
    sqrtw_J = sqrt.(weights[J])

    @views for b in eachindex(J), a in eachindex(I)
        i = I[a]
        j = J[b]

        q = hypot(
            p_points[i][1] - p_points[j][1],
            p_points[i][2] - p_points[j][2]
        )

        B[a, b] = sqrtw_I[a] * sqrtw_J[b] *
                  Vscreened(q, εr, interaction_type) *
                  dot(Uminus[:, j], Uminus[:, i]) *
                  dot(Uplus[:, i], Uplus[:, j])
    end

    return B
end


function aca_block(I, J, p_points, weights,
                   Uminus, Uplus, εr, interaction_type; maxrank=20)
    @assert maxrank >= 1

    m, n = length(I), length(J)
    rank_limit = min(maxrank, m, n)

    L = zeros(ComplexF64, m, rank_limit)
    R = zeros(ComplexF64, rank_limit, n)

    used_rows = falses(m)
    pivot_columns = Int[]

    i = 1
    rank = 0

    for step in 1:rank_limit
        row = vec(exchange_block(
            I[i:i], J, p_points, weights,
            Uminus, Uplus, εr, interaction_type
        ))
        row_scale = maximum(abs.(row))

        if rank > 0
            row -= vec(transpose(L[i, 1:rank]) * R[1:rank, :])
        end

        j = argmax(abs.(row))
        pivot = row[j]

        if abs(pivot) <= 100 * eps(Float64) * row_scale
            break
        end

        column = vec(exchange_block(
            I, J[j:j], p_points, weights,
            Uminus, Uplus, εr, interaction_type
        ))

        if rank > 0
            column -= L[:, 1:rank] * R[1:rank, j]
        end

        rank += 1
        L[:, rank] = column / pivot
        R[rank, :] = row

        # Record the pivots used to construct the approximation.
        used_rows[i] = true
        push!(pivot_columns, j)

        scores = abs.(column)
        scores[used_rows] .= -Inf
        i = argmax(scores)
    end

    pivot_rows = findall(used_rows)

    return L[:, 1:rank], R[1:rank, :], pivot_rows, pivot_columns
end

function check_aca(I, J, L, R, pivot_rows, pivot_columns,
                   p_points, weights, Uminus, Uplus, εr, interaction_type;
                   tol=1e-8, ncheck=8)
    @assert tol > 0 && ncheck >= 1

    m, n = length(I), length(J)

    available_rows = setdiff(1:m, pivot_rows)
    columns = setdiff(1:n, pivot_columns)

    # No independent entries remain: do not accept compression
    # on the strength of an empty validation sample.
    if isempty(available_rows) || isempty(columns)
        return false, Inf
    end

    nsample = min(ncheck, length(available_rows))
    selection = randperm(length(available_rows))[1:nsample]
    rows = sort(available_rows[selection])

    exact = exchange_block(
        I[rows], J[columns], p_points, weights,
        Uminus, Uplus, εr, interaction_type
    )

    approximation = L[rows, :] * R[:, columns]

    error = norm(exact - approximation)
    scale = norm(exact)

    if scale == 0
        relative_error = error == 0 ? 0.0 : Inf
    else
        relative_error = error / scale
    end

    return relative_error <= tol, relative_error
end

struct ExchangeBlock
    I::UnitRange{Int}
    J::UnitRange{Int}
    L::Matrix{ComplexF64}
    R::Union{Nothing, Matrix{ComplexF64}}
    diagonal::Bool
end

function build_exchange_blocks(
    tree, p_points, weights, Uminus, Uplus, εr, interaction_type;
    tol=1e-8, maxrank=20, η=1.0, ncheck=8
)
    blocks = ExchangeBlock[]

    function visit(A, B)
        I = A.positions
        J = B.positions
        diagonal = I == J

        # Try compression for sufficiently separated clusters.
        if !diagonal && well_separated(A, B; η=η)
            L, R, pivot_rows, pivot_columns = aca_block(
                I, J, p_points, weights,
                Uminus, Uplus, εr, interaction_type;
                maxrank=maxrank
            )

            # Only keep factors if they use less storage.
            if length(L) + length(R) < length(I) * length(J)
                passed, _ = check_aca(
                    I, J, L, R, pivot_rows, pivot_columns,
                    p_points, weights, Uminus, Uplus, εr, interaction_type;
                    tol=tol, ncheck=ncheck
                )

                if passed
                    push!(blocks, ExchangeBlock(I, J, L, R, false))
                    return
                end
            end
        end

        # Both clusters are leaves: store the block densely.
        if isempty(A.children) && isempty(B.children)
            dense = exchange_block(
                I, J, p_points, weights,
                Uminus, Uplus, εr, interaction_type
            )

            push!(blocks, ExchangeBlock(I, J, dense, nothing, diagonal))

        elseif diagonal
            # The other off-diagonal block is the adjoint.
            A1, A2 = A.children

            visit(A1, A1)
            visit(A1, A2)
            visit(A2, A2)

        elseif !isempty(A.children) &&
               (isempty(B.children) || length(I) >= length(J))
            # Split the larger available cluster.
            for child in A.children
                visit(child, B)
            end

        else
            for child in B.children
                visit(A, child)
            end
        end
    end

    visit(tree, tree)
    return blocks
end


function apply_compressed_hessian(
    x, blocks, Q, weights, D, FQ, εr, interaction_type
)
    g = sqrt.(weights) .* FQ

    # Exact kinetic and Hartree terms.
    y = D .* x + Vscreened(norm(Q), εr, interaction_type) * g * dot(g, x)

    @views for b in blocks
        if b.R === nothing
            # Dense exchange block.
            y[b.I] .-= b.L * x[b.J]

            if !b.diagonal
                y[b.J] .-= b.L' * x[b.I]
            end
        else
            # Compressed exchange block: B ≈ L * R.
            y[b.I] .-= b.L * (b.R * x[b.J])

            if !b.diagonal
                y[b.J] .-= b.R' * (b.L' * x[b.I])
            end
        end
    end

    return y
end




function get_compressed_problem(
    Q, μ, uD, valley, stacking, NL, band, εr, interaction_type;
    Nθ=96, Nr=8, rmax=0.8, nscan=800, r_tol=1e-10,
    modeltype=:spinor,
    leafsize=64, compression_tol=1e-8,
    maxrank=20, η=1.0, ncheck=8
)
    @assert modeltype in (:spinor, :no_spinor)

    # Construct the midpoint grid and integration weights.
    p_points, weights = polar_grid(
        Q, μ, uD, valley, stacking, NL, band;
        Nθ=Nθ, Nr=Nr, rmax=rmax,
        nscan=nscan, r_tol=r_tol
    )

    @assert !isempty(p_points) "No active pairs were found."

    # Compute energy differences and endpoint spinors.
    D, Uminus, Uplus, FQ = pair_data(
        p_points, Q, uD, valley, stacking, NL, band
    )

    # For the no-spinor model, retain the same dispersion
    # but replace all spinor overlaps by one.
    if modeltype == :no_spinor
        N = length(p_points)
        Uminus = ones(ComplexF64, 1, N)
        Uplus = ones(ComplexF64, 1, N)
        FQ = ones(ComplexF64, N)
    end

    # Group nearby momenta.
    order, tree = momentum_tree(p_points; leafsize=leafsize)

    # Put every array into the same ordering.
    p_points = p_points[order]
    weights = weights[order]
    D = D[order]
    Uminus = Uminus[:, order]
    Uplus = Uplus[:, order]
    FQ = FQ[order]

    # Construct the compressed exchange representation.
    blocks = build_exchange_blocks(
        tree, p_points, weights, Uminus, Uplus, εr, interaction_type;
        tol=compression_tol,
        maxrank=maxrank,
        η=η,
        ncheck=ncheck
    )

    return (;
        Q, εr, interaction_type,
        p_points, weights, D, Uminus, Uplus, FQ,
        order, blocks
    )
end


function get_compressed_modes(
    problem, nmodes;
    eig_tol=1e-10, krylovdim=60, maxiter=400
)
    N = length(problem.D)

    @assert 1 <= nmodes < N
    @assert eig_tol > 0 && maxiter >= 1

    krylovdim = min(krylovdim, N)
    @assert krylovdim > nmodes

    # The operator acting on weighted solver vectors.
    H_action = x -> apply_compressed_hessian(
        x, problem.blocks, problem.Q,
        problem.weights, problem.D, problem.FQ,
        problem.εr, problem.interaction_type
    )

    x0 = randn(ComplexF64, N)
    x0 /= norm(x0)

    values, vectors, info = eigsolve(
        H_action, x0, nmodes, :SR;
        ishermitian=true,
        tol=eig_tol,
        krylovdim=krylovdim,
        maxiter=maxiter
    )

    @assert info.converged >= nmodes (
        "Not all requested eigenpairs converged."
    )

    # Keep the requested lowest eigenvalues in ascending order.
    selected = sortperm(real.(values))[1:nmodes]
    eigenvalues = real.(values[selected])

    solver_modes = zeros(ComplexF64, N, nmodes)
    modes = zeros(ComplexF64, N, nmodes)

    residuals = zeros(Float64, nmodes)
    normalizations = zeros(Float64, nmodes)
    dQ = zeros(ComplexF64, nmodes)

    sqrtw = sqrt.(problem.weights)

    for a in 1:nmodes
        λ = eigenvalues[a]

        # Ordinary normalized eigenvector of the weighted Hessian.
        v_tilde = vectors[selected[a]]
        v_tilde = v_tilde / norm(v_tilde)

        solver_modes[:, a] = v_tilde

        # Convert to the physical mode shape.
        modes[:, a] = v_tilde ./ sqrtw

        # Physical normalization: sum_i W_i |X_i|² = 1.
        normalizations[a] = sum(
            problem.weights .* abs2.(modes[:, a])
        )

        # Charge overlap: sum_i W_i conj(FQ_i) X_i.
        dQ[a] = dot(
            problem.FQ,
            problem.weights .* modes[:, a]
        )

        # The weighted Hessian acts on v_tilde, not directly on X.
        residuals[a] = norm(
            H_action(v_tilde) - λ * v_tilde
        )
    end

    @assert all(isfinite, eigenvalues) && all(isfinite, residuals)
    @assert maximum(residuals) <= 10 * eig_tol (
        "An eigenpair residual exceeds 10 * eig_tol."
    )

    return (;
        eigenvalues, solver_modes, modes,
        normalizations, dQ, residuals
    )
end


function get_mu(;
    uD, n_h,
    NL=6, valley=1, stacking=1, band=NL,
    kmax=1.0, Nk=501
)
    @assert kmax > 0 && Nk >= 2
    @assert n_h > 0
    @assert 1 <= band <= 2NL

    # Uniform Cartesian grid for determining μ.
    dk = 2kmax / Nk
    kgrid = range(
        -kmax + dk/2, kmax - dk/2;
        length=Nk
    )
    w = dk^2 / (2π)^2

    # Rows correspond to ky; columns correspond to kx.
    εgrid = [
        eigvals(Hermitian(
            get_Ham([kx, ky], uD, valley, stacking, NL)
        ))[band]
        for ky in kgrid, kx in kgrid
    ]

    # Holes occupy the highest-energy states of this band.
    εsorted = sort(vec(εgrid); rev=true)
    Nholes = round(Int, n_h / w)

    @assert 1 <= Nholes < length(εsorted) (
        "Density is outside the grid's range."
    )

    # Place μ between the last hole state and
    # the first occupied electron state.
    μ = (εsorted[Nholes] + εsorted[Nholes + 1]) / 2

    fgrid = εgrid .< μ
    n_h_grid = w * count(.!fgrid)
    Nholes = count(.!fgrid)

    # The sampled hole region should stay inside the box.
    @assert all(fgrid[1, :]) && all(fgrid[end, :]) &&
            all(fgrid[:, 1]) && all(fgrid[:, end]) (
        "Increase kmax: holes reach the grid boundary."
    )

    return (; μ, n_h_grid, kgrid, εgrid, fgrid, dk, w, Nholes)
end