using LinearAlgebra
using Arpack
using Combinatorics
using Random
using StaticArrays
using JLD2
BLAS.set_num_threads(1)




function get_Ham(k::Vector{Float64}, λ::Float64, Δ::Float64, rs::Float64)
    kx = k[1]
    ky = k[2]

    return -Δ * [
        λ^2 * norm(k)^2      -λ * (kx - im * ky)
        -λ * (kx + im * ky)   1
    ] + norm(k)^2 / rs^2 * Matrix{Float64}(I, 2, 2)
end


# ============================================================
# State indexing
#
# g = plane-wave/momentum index
# a = local active-band index, a=1:length(kept_bands)
# kept_bands[a] is the physical band label, either 1 or 2.
# ============================================================

@inline function state_index(g::Int, a::Int, nband::Int)::Int
    return a + (g - 1) * nband
end

@inline function state_g(I::Int, nband::Int)::Int
    return (I - 1) ÷ nband + 1
end

@inline function state_band(I::Int, nband::Int)::Int
    return (I - 1) % nband + 1
end



# ============================================================
# kept_bands helpers
#
# Static saved args convention:
# args = [
#   gcutoff, λ, whether_kept[1], whether_kept[2],
#   enlarge_factor, V2_scalar, ϕ, pin_coeff,
#   rs, dedis, defec_pos, constq, filling, Δ
# ]
# ============================================================

function kept_bands_from_flags(keep1::Real, keep2::Real)::Vector{Int64}
    kept = Int64[]

    if Int(round(keep1)) == 1
        push!(kept, 1)
    end

    if Int(round(keep2)) == 1
        push!(kept, 2)
    end

    @assert !isempty(kept) "No active band selected."
    return kept
end


function kept_bands_from_saved_args(seed_args)::Vector{Int64}
    return kept_bands_from_flags(seed_args[3], seed_args[4])
end

function kept_flags_from_bands(kept_bands::Vector{Int64})
    keep1 = 1 in kept_bands ? 1.0 : 0.0
    keep2 = 2 in kept_bands ? 1.0 : 0.0
    return keep1, keep2
end


# ============================================================
# Reference density
#
# This is always the simple diagonal 1/0 matrix in the current
# active basis.
#
# kept_bands=[1]   -> Pref = I
# kept_bands=[2]   -> Pref = 0
# kept_bands=[1,2] -> Pref = 1 on physical band 1 block, 0 on band 2
# ============================================================

function build_Pref(
    wave::Vector{Vector{Int64}},
    kept_bands::Vector{Int64},
)::Matrix{ComplexF64}

    Ng = length(wave)
    nband = length(kept_bands)
    Pref = zeros(ComplexF64, Ng * nband, Ng * nband)

    local_band_1 = findfirst(==(1), kept_bands)

    if local_band_1 !== nothing
        for g in 1:Ng
            I = state_index(g, local_band_1, nband)
            Pref[I, I] = 1.0 + 0.0im
        end
    end

    return Pref
end


# ============================================================
# Hermitian cleanup
# ============================================================

@inline function hermitize!(A::Matrix{ComplexF64})
    A .= (A .+ A') ./ 2
    return A
end

@inline function symmetrize_from_lower!(A::Matrix{ComplexF64})
    n = size(A, 1)

    @inbounds for i in 1:n
        A[i, i] = complex(real(A[i, i]), 0.0)

        for j in (i + 1):n
            A[i, j] = conj(A[j, i])
        end
    end

    return A
end


# ============================================================
# Stable potentials
# ============================================================

function Coulomb(k::Vector{Int64}, T1::Vector{Float64}, T2::Vector{Float64})::Float64
    q = norm([T1 T2] * k)
    return q < 1e-14 ? 0.0 : 4π / q
end

function get_Fourier_potential(
    k::Vector{Int64},
    T1::Vector{Float64},
    T2::Vector{Float64},
    d::Float64,
)
    kvec = [T1 T2] * k
    q = norm(kvec)
    D = 25.0

    if q < 1e-14
        return (D - d) * 4π
    else
        # stable version of sinh(q*(D-d)) / cosh(q*D)
        ratio = (exp(-q * d) - exp(-q * (2D - d))) / (1.0 + exp(-2q * D))
        return ratio / q * 4π
    end
end



struct WaveLookup
    pos::Matrix{Int64}   # 0 means “missing”
    n1min::Int64
    n2min::Int64
end

@inline function lookup(wl::WaveLookup, n1::Int, n2::Int)::Int64
    i1 = n1 - wl.n1min + 1
    i2 = n2 - wl.n2min + 1
     #bounds check + read
    return (1 <= i1 <= size(wl.pos,1) && 1 <= i2 <= size(wl.pos,2)) ? wl.pos[i1,i2] : Int64(0)
end



function build_wave_lookup(wave::Vector{Vector{Int64}})
    Ng = length(wave)
    n1min = minimum(w -> w[1], wave)
    n1max = maximum(w -> w[1], wave)
    n2min = minimum(w -> w[2], wave)
    n2max = maximum(w -> w[2], wave)

    
 

    pos = fill(Int64(0), n1max - n1min + 1, n2max - n2min + 1)
    @inbounds for p in 1:Ng
        pos[wave[p][1] - n1min + 1, wave[p][2] - n2min + 1] = Int64(p)
    end
    return WaveLookup(pos, n1min, n2min)
end



struct ShiftIndexer
    dn1min::Int
    dn1max::Int
    dn2min::Int
    dn2max::Int
    S1::Int
    NSHIFT::Int
end

function ShiftIndexer(wave_n1::Vector{Int64}, wave_n2::Vector{Int64})
    n1min = Int(minimum(wave_n1));  n1max = Int(maximum(wave_n1))
    n2min = Int(minimum(wave_n2));  n2max = Int(maximum(wave_n2))

    dn1min = n1min - n1max
    dn1max = n1max - n1min
    dn2min = n2min - n2max
    dn2max = n2max - n2min

    S1 = dn1max - dn1min + 1
    S2 = dn2max - dn2min + 1
    NSHIFT = S1 * S2

    return ShiftIndexer(dn1min, dn1max, dn2min, dn2max, S1, NSHIFT)
end

@inline function shift_id(ix::ShiftIndexer, dn1::Int, dn2::Int)::Int
    return (dn1 - ix.dn1min + 1) + (dn2 - ix.dn2min) * ix.S1
end

struct ShiftCSR
    ix::ShiftIndexer
    offsets::Vector{Int64}   # length NSHIFT+1
    g2_list::Vector{Int64}
    g3_list::Vector{Int64}
end


function build_shiftcsr(wl::WaveLookup, wave_n1::Vector{Int64}, wave_n2::Vector{Int64})
    Ng = length(wave_n1)
    ix = ShiftIndexer(wave_n1, wave_n2)
    NSHIFT = ix.NSHIFT

    # -------- PASS 1: counts per shift --------
    counts = fill(Int64(0), NSHIFT)

    for dn2 in ix.dn2min:ix.dn2max
        for dn1 in ix.dn1min:ix.dn1max
            sid = shift_id(ix, dn1, dn2)
            c = Int64(0)
            @inbounds for g2 in 1:Ng
                g3 = lookup(wl, wave_n1[g2] + dn1, wave_n2[g2] + dn2)
                c += (g3 != 0)
            end
            counts[sid] = c
        end
    end

    # -------- offsets (prefix sum) --------
    offsets = Vector{Int}(undef, NSHIFT + 1)
    offsets[1] = 1
    @inbounds for sid in 1:NSHIFT
        offsets[sid+1] = offsets[sid] + counts[sid]
    end

    total = Int(offsets[end] - 1)
    g2_list = Vector{Int64}(undef, total)
    g3_list = Vector{Int64}(undef, total)

    # -------- PASS 2: fill packed arrays --------
    for dn2 in ix.dn2min:ix.dn2max
        for dn1 in ix.dn1min:ix.dn1max
            sid = shift_id(ix, dn1, dn2)
            p = Int(offsets[sid])
            @inbounds for g2 in 1:Ng
                g3 = lookup(wl, wave_n1[g2] + dn1, wave_n2[g2] + dn2)
                if g3 != 0
                    g2_list[p] = Int(g2)
                    g3_list[p] = g3
                    p += 1
                end
            end
        end
    end

    return ShiftCSR(ix, offsets, g2_list, g3_list)
end

@inline function csr_range(csr::ShiftCSR, dn1::Int, dn2::Int)
    sid = shift_id(csr.ix, dn1, dn2)
    lo = Int(csr.offsets[sid])
    hi = Int(csr.offsets[sid+1]) - 1
    return lo, hi
end

function csr_stats(csr::ShiftCSR, Ng::Int)
    total_pairs = Int(csr.offsets[end] - 1)
    mean_pairs = total_pairs / csr.ix.NSHIFT
    println("CSR stats: NSHIFT=$(csr.ix.NSHIFT), total_pairs=$total_pairs, mean_pairs/shift=$mean_pairs, Ng=$Ng")
end



function make_wave_A(
    T1::Vector{Float64},
    T2::Vector{Float64},
    b1T::Vector{Int64},
    b2T::Vector{Int64},
    gcutoff_work::Float64,
    Ashift::Vector{Float64},
)
    wave = Vector{Int64}[]

    cutoffstandard = gcutoff_work * norm(T1)
    cutoff = Int(ceil(gcutoff_work)) * 6

    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        w = ja * b1T + jb * b2T
        p = [T1 T2] * w + Ashift

        if dot(p, p) < cutoffstandard^2
            push!(wave, Int64.(w))
        end
    end

    return wave
end


# ============================================================
# Static lattice geometry from lambda-jellium conventions
# ============================================================

function lambda_geometry(enlarge_factor::Int)
    LLL = sqrt(2π / sqrt(3))

    b1 = 4π / (sqrt(3) * LLL) * [1.0, 0.0] / enlarge_factor
    b2 = 4π / (sqrt(3) * LLL) * [-0.5, sqrt(3) / 2] / enlarge_factor

    T1 = b1
    T2 = b2

    a1m = inv([b1'; b2']) * [2π, 0.0]
    a2m = inv([b1'; b2']) * [0.0, 2π]

    am = norm(a1m)
    Area = sqrt(3) / 2 * am^2

    b1T = Int64.(round.(inv([T1 T2]) * b1))
    b2T = Int64.(round.(inv([T1 T2]) * b2))

    return T1, T2, a1m, a2m, b1, b2, b1T, b2T, Area
end



# ============================================================
# Basis bookkeeping
# ============================================================

function build_basis_bookkeeping(
    wave::Vector{Vector{Int64}},
    T1::Vector{Float64},
    T2::Vector{Float64},
    rs::Float64,
    constq::Float64,
)
    Ng = length(wave)

    wave_n1 = Vector{Int64}(undef, Ng)
    wave_n2 = Vector{Int64}(undef, Ng)

    @inbounds for g in 1:Ng
        wave_n1[g] = wave[g][1]
        wave_n2[g] = wave[g][2]
    end

    wave_dict = Dict{Tuple{Int64, Int64}, Int64}()

    @inbounds for g in eachindex(wave)
        wave_dict[(wave[g][1], wave[g][2])] = Int64(g)
    end

    wave_diff = Vector{Vector{Int64}}()
    sizehint!(wave_diff, Ng^2)

    @inbounds for g1 in eachindex(wave), g2 in eachindex(wave)
        push!(wave_diff, wave[g1] - wave[g2])
    end

    wave_diff = unique(wave_diff)

    wl = build_wave_lookup(wave)
    csr = build_shiftcsr(wl, wave_n1, wave_n2)

    Coulomb_transfer = Matrix{Float64}(undef, Ng, Ng)

    @inbounds for g1 in 1:Ng, g3 in 1:Ng
        Coulomb_transfer[g1, g3] =
            Coulomb(wave[g1] - wave[g3], T1, T2) / rs + constq
    end

    return (
        wave_n1 = wave_n1,
        wave_n2 = wave_n2,
        wave_dict = wave_dict,
        wave_diff = wave_diff,
        csr = csr,
        Coulomb_transfer = Coulomb_transfer,
    )
end



# ============================================================
# Seed-density embedding
#
# Input P_seed is relative density from static HF.
# Output P_new is relative density in TDHF active basis.
#
# This allows:
#   [1]   -> [1]
#   [2]   -> [2]
#   [1,2] -> [1,2]
#   [1]   -> [1,2]
#   [2]   -> [1,2]
#
# It forbids [1,2] -> [1] or [2] unless explicitly allowed.
# ============================================================

function remap_seed_density_to_tdhf_basis(
    P_seed::Matrix{ComplexF64},
    wave_seed::Vector{Vector{Int64}},
    kept_seed::Vector{Int64},
    wave_new::Vector{Vector{Int64}},
    kept_new::Vector{Int64};
    allow_seed_truncation::Bool = false,
)
    nband_seed = length(kept_seed)
    nband_new = length(kept_new)

    @assert size(P_seed) == (length(wave_seed) * nband_seed,
                             length(wave_seed) * nband_seed)

    discarded = setdiff(kept_seed, kept_new)

    if !isempty(discarded) && !allow_seed_truncation
        error(
            "Seed kept_bands=$(kept_seed), TDHF kept_bands=$(kept_new). " *
            "This would discard band(s) $(discarded)."
        )
    end

    P_new = zeros(ComplexF64, length(wave_new) * nband_new,
                              length(wave_new) * nband_new)

    old_g_index = Dict{Tuple{Int64, Int64}, Int64}()

    for g in eachindex(wave_seed)
        old_g_index[(wave_seed[g][1], wave_seed[g][2])] = Int64(g)
    end

    old_local = Dict{Int64, Int64}()
    for a in eachindex(kept_seed)
        old_local[kept_seed[a]] = Int64(a)
    end

    new_local = Dict{Int64, Int64}()
    for a in eachindex(kept_new)
        new_local[kept_new[a]] = Int64(a)
    end

    for gnew1 in eachindex(wave_new), gnew2 in eachindex(wave_new)
        gold1 = get(old_g_index, (wave_new[gnew1][1], wave_new[gnew1][2]), 0)
        gold2 = get(old_g_index, (wave_new[gnew2][1], wave_new[gnew2][2]), 0)

        if gold1 == 0 || gold2 == 0
            continue
        end

        for phys_b1 in kept_seed
            if !haskey(new_local, phys_b1)
                continue
            end

            for phys_b2 in kept_seed
                if !haskey(new_local, phys_b2)
                    continue
                end

                aold1 = old_local[phys_b1]
                aold2 = old_local[phys_b2]

                anew1 = new_local[phys_b1]
                anew2 = new_local[phys_b2]

                Iold = state_index(gold1, aold1, nband_seed)
                Jold = state_index(gold2, aold2, nband_seed)

                Inew = state_index(gnew1, anew1, nband_new)
                Jnew = state_index(gnew2, anew2, nband_new)

                P_new[Inew, Jnew] = P_seed[Iold, Jold]
            end
        end
    end

    hermitize!(P_new)
    return P_new
end


# ============================================================
# Relative-density transport
#
# This is the convention we agreed on:
# transport Q=P+Pref, then subtract the destination Pref.
# ============================================================

function transport_relative_density(
    P_old::Matrix{ComplexF64},
    Pref_old::Matrix{ComplexF64},
    Pref_new::Matrix{ComplexF64},
    S::Matrix{ComplexF64},
)
    Q_old = P_old + Pref_old
    Q_new = S * Q_old * S'
    P_new = Q_new - Pref_new

    hermitize!(P_new)
    return P_new
end



# ============================================================
# TDHF parameters
# ============================================================

struct TDHFParams
    # lambda-jellium/static parameters
    gcutoff_seed::Float64
    λ::Float64
    enlarge_factor::Int
    V2_scalar::Float64
    ϕ::Float64
    pin_coeff::Float64
    rs::Float64
    dedis::Float64
    defec_pos::Int
    constq::Float64
    filling::Int
    Δ::Float64

    # TDHF active basis
    kept_bands::Vector{Int64}
    gcutoff_work::Float64

    # time evolution
    gamma::Float64
    temp::Float64
    dt::Float64
    dAshift_dt::Vector{Float64}

    # geometry
    T1::Vector{Float64}
    T2::Vector{Float64}
    a1m::Vector{Float64}
    a2m::Vector{Float64}
    b1::Vector{Float64}
    b2::Vector{Float64}
    b1T::Vector{Int64}
    b2T::Vector{Int64}
    Area::Float64
end


function make_tdhf_params(
    static_args::Vector{Float64},
    kept_bands_tdhf::Vector{Int64};
    gcutoff_work::Float64 = static_args[1],
    gamma::Float64,
    temp::Float64,
    dt::Float64,
    dAshift_dt::Vector{Float64},
)
    gcutoff_seed  = static_args[1]
    λ             = static_args[2]
    enlarge_factor = Int(round(static_args[5]))
    V2_scalar     = static_args[6]
    ϕ             = static_args[7]
    pin_coeff     = static_args[8]
    rs            = static_args[9]
    dedis         = static_args[10]
    defec_pos     = Int(round(static_args[11]))
    constq        = static_args[12]
    filling       = Int(round(static_args[13]))
    Δ             = static_args[14]

    T1, T2, a1m, a2m, b1, b2, b1T, b2T, Area =
        lambda_geometry(enlarge_factor)

    return TDHFParams(
        gcutoff_seed,
        λ,
        enlarge_factor,
        V2_scalar,
        ϕ,
        pin_coeff,
        rs,
        dedis,
        defec_pos,
        constq,
        filling,
        Δ,
        kept_bands_tdhf,
        gcutoff_work,
        gamma,
        temp,
        dt,
        dAshift_dt,
        T1,
        T2,
        a1m,
        a2m,
        b1,
        b2,
        b1T,
        b2T,
        Area,
    )
end



# ============================================================
# TDHF work buffers
# ============================================================

mutable struct TDHFWork
    HartreeMatrix::Matrix{ComplexF64}
    FockMatrix::Matrix{ComplexF64}
    H_phys::Matrix{ComplexF64}
    HartreeAccShift::Vector{ComplexF64}
end

function TDHFWork(dimension::Int, NSHIFT::Int)
    Z = zeros(ComplexF64, dimension, dimension)

    return TDHFWork(
        copy(Z),
        copy(Z),
        copy(Z),
        zeros(ComplexF64, NSHIFT),
    )
end


# ============================================================
# Defect position
# Same convention as static lambda code.
# ============================================================

function defect_position_vec(
    a1m::Vector{Float64},
    a2m::Vector{Float64},
    enlarge_factor::Int,
    defec_pos::Int,
)
    if defec_pos == 1
        return (2 * a2m - a1m) / enlarge_factor * 2 / 3
    elseif defec_pos == 2
        return (2 * a2m - a1m) / enlarge_factor * 1 / 3
    elseif defec_pos == 3
        return (2 * a2m - a1m) / enlarge_factor * 0
    elseif defec_pos == 4
        return (2 * a2m - a1m) / enlarge_factor * (-2 / 3)
    else
        error("Unknown defec_pos = $(defec_pos)")
    end
end





# ============================================================
# Rebuild instantaneous projected data at a given Ashift
#
# Returns:
#   spinor_set[g,a]
#   overlapmatrix[I,J]
#   single_Ham
#   single_MoirePo
#   pinning_po
#   Pref
# ============================================================

function refresh_projected_data(
    wave::Vector{Vector{Int64}},
    book,
    Ashift::Vector{Float64},
    params::TDHFParams,
)
    kept_bands = params.kept_bands
    nband = length(kept_bands)
    Ng = length(wave)
    dimension = Ng * nband

    T1 = params.T1
    T2 = params.T2

    λ = params.λ
    Δ = params.Δ
    rs = params.rs

    enlarge_factor = params.enlarge_factor
    V2_scalar = params.V2_scalar
    ϕ = params.ϕ

    pin_coeff = params.pin_coeff
    dedis = params.dedis
    defec_pos = params.defec_pos

    a1m = params.a1m
    a2m = params.a2m
    b1T = params.b1T
    b2T = params.b2T
    Area = params.Area

    # ------------------------------------------------------------
    # Spinors and kinetic Hamiltonian
    # ------------------------------------------------------------

    spinor_set = Matrix{Vector{ComplexF64}}(undef, Ng, nband)
    single_Ham = zeros(ComplexF64, dimension, dimension)

    for g in 1:Ng
        kvec = [T1 T2] * wave[g] + Ashift
        F = eigen(Hermitian(get_Ham(kvec, λ, Δ, rs)))

        for a in 1:nband
            physical_band = kept_bands[a]

            v = F.vectors[:, physical_band]

            # deterministic phase convention, same style as static lambda code
            imax = argmax(abs.(v))
            v = v * exp(-im * angle(v[imax])) / norm(v)

            spinor_set[g, a] = v

            I = state_index(g, a, nband)
            single_Ham[I, I] = real(F.values[physical_band])
        end
    end

    # ------------------------------------------------------------
    # Form-factor matrix
    # ------------------------------------------------------------

    overlapmatrix = zeros(ComplexF64, dimension, dimension)

    for g1 in 1:Ng, a1 in 1:nband
        I1 = state_index(g1, a1, nband)

        for g2 in 1:Ng, a2 in 1:nband
            I2 = state_index(g2, a2, nband)
            overlapmatrix[I1, I2] =
                spinor_set[g1, a1]' * spinor_set[g2, a2]
        end
    end

    # ------------------------------------------------------------
    # Scalar moire potential
    # ------------------------------------------------------------

    single_MoirePo = zeros(ComplexF64, dimension, dimension)
    op_5 = Matrix{ComplexF64}(I, 2, 2)

    for jc in 1:Ng
        # -b1
        pos = get(book.wave_dict, (wave[jc][1] - enlarge_factor * b1T[1],
                                   wave[jc][2] - enlarge_factor * b1T[2]), 0)

        if pos != 0
            for a in 1:nband, b in 1:nband
                I = state_index(jc, a, nband)
                J = state_index(pos, b, nband)

                single_MoirePo[I, J] +=
                    V2_scalar * exp(-im * ϕ) *
                    (spinor_set[jc, a]' * op_5 * spinor_set[pos, b])
            end
        end

        # -b2
        pos = get(book.wave_dict, (wave[jc][1] - enlarge_factor * b2T[1],
                                   wave[jc][2] - enlarge_factor * b2T[2]), 0)

        if pos != 0
            for a in 1:nband, b in 1:nband
                I = state_index(jc, a, nband)
                J = state_index(pos, b, nband)

                single_MoirePo[I, J] +=
                    V2_scalar * exp(-im * ϕ) *
                    (spinor_set[jc, a]' * op_5 * spinor_set[pos, b])
            end
        end

        # +(b1+b2)
        pos = get(book.wave_dict, (wave[jc][1] + enlarge_factor * (b1T[1] + b2T[1]),
                                   wave[jc][2] + enlarge_factor * (b1T[2] + b2T[2])), 0)

        if pos != 0
            for a in 1:nband, b in 1:nband
                I = state_index(jc, a, nband)
                J = state_index(pos, b, nband)

                single_MoirePo[I, J] +=
                    V2_scalar * exp(-im * ϕ) *
                    (spinor_set[jc, a]' * op_5 * spinor_set[pos, b])
            end
        end
    end

    single_MoirePo = single_MoirePo + single_MoirePo'

    # ------------------------------------------------------------
    # Defect/pinning potential
    # ------------------------------------------------------------

    pinning_po = zeros(ComplexF64, dimension, dimension)

    if abs(pin_coeff) > 0
        defec_pos_vec = defect_position_vec(a1m, a2m, enlarge_factor, defec_pos)

        for jc in 1:Ng
            for diff in book.wave_diff
                pos = get(book.wave_dict, (wave[jc][1] + diff[1],
                                           wave[jc][2] + diff[2]), 0)

                if pos != 0
                    qvec = [T1 T2] * diff
                    phase = exp(im * dot(qvec, defec_pos_vec))

                    prefactor =
                        pin_coeff / (Area * rs) *
                        get_Fourier_potential(diff, T1, T2, dedis)

                    for a in 1:nband, b in 1:nband
                        I = state_index(pos, a, nband)
                        J = state_index(jc, b, nband)

                        pinning_po[I, J] += prefactor *
                            (spinor_set[pos, a]' * op_5 * spinor_set[jc, b]) *
                            phase
                    end
                end
            end
        end
    end

    pinning_po = (pinning_po + pinning_po') / 2

    # ------------------------------------------------------------
    # Reference matrix
    # ------------------------------------------------------------

    Pref = build_Pref(wave, kept_bands)

    return (
        spinor_set = spinor_set,
        overlapmatrix = overlapmatrix,
        single_Ham = single_Ham,
        single_MoirePo = single_MoirePo,
        pinning_po = pinning_po,
        Pref = Pref,
    )
end


# ============================================================
# Build instantaneous HF Hamiltonian for lambda-jellium TDHF
#
# Input:
#   P = relative density matrix
#
# Output stored in work:
#   work.HartreeMatrix
#   work.FockMatrix
#   work.H_phys = single_Ham + single_MoirePo + pinning_po
#                 + HartreeMatrix - FockMatrix
#
# This follows the static lambda-jellium Construct_DensityMatrix
# Hartree/Fock convention, but does not diagonalize and does not
# update P.
# ============================================================

function Build_HHF_lambda!(
    work::TDHFWork,
    csr::ShiftCSR,
    wave::Vector{Vector{Int64}},
    wave_n1::Vector{Int64},
    wave_n2::Vector{Int64},
    P::Matrix{ComplexF64},
    single_Ham::Matrix{ComplexF64},
    single_MoirePo::Matrix{ComplexF64},
    pinning_po::Matrix{ComplexF64},
    overlapmatrix::Matrix{ComplexF64},
    Area::Float64,
    Coulomb_transfer::Matrix{Float64},
    nband::Int,
)
    Ng = length(wave)
    dimension = Ng * nband

    @assert size(P) == (dimension, dimension)
    @assert size(single_Ham) == (dimension, dimension)
    @assert size(single_MoirePo) == (dimension, dimension)
    @assert size(pinning_po) == (dimension, dimension)
    @assert size(overlapmatrix) == (dimension, dimension)
    @assert size(Coulomb_transfer) == (Ng, Ng)

    fill!(work.HartreeMatrix, 0)
    fill!(work.FockMatrix, 0)

    # ------------------------------------------------------------
    # Fock term
    #
    # Static lambda convention:
    #
    # Fock[I1,I4] =
    #   sum_{g2,g3,a2,a3}
    #       P[I3,I2] *
    #       V(g1-g3) *
    #       F[I1,I3] *
    #       F[I2,I4]
    #
    # where momentum conservation is handled by CSR:
    #   g1 + g2 = g3 + g4
    # ------------------------------------------------------------

    Threads.@threads :greedy for I1 in 1:dimension
        @inbounds begin
            g1 = state_g(I1, nband)

            w1n1 = wave_n1[g1]
            w1n2 = wave_n2[g1]

            for I4 in 1:I1
                g4 = state_g(I4, nband)

                dg_n1 = w1n1 - wave_n1[g4]
                dg_n2 = w1n2 - wave_n2[g4]

                lo, hi = csr_range(csr, dg_n1, dg_n2)

                acc = 0.0 + 0.0im

                for p in lo:hi
                    g2 = Int(csr.g2_list[p])
                    g3 = Int(csr.g3_list[p])

                    V13 = Coulomb_transfer[g1, g3]

                    for a2 in 1:nband, a3 in 1:nband
                        I2 = state_index(g2, a2, nband)
                        I3 = state_index(g3, a3, nband)

                        acc += P[I3, I2] *
                            V13 *
                            overlapmatrix[I1, I3] *
                            overlapmatrix[I2, I4]
                    end
                end

                work.FockMatrix[I1, I4] = acc
            end
        end
    end

    # ------------------------------------------------------------
    # Hartree term
    #
    # First build density Fourier component for each shift:
    #
    #   accShift[shift] =
    #       sum_{g2,g4,a2,a4 with g4-g2=shift}
    #           P[I4,I2] * F[I2,I4]
    #
    # Then:
    #
    #   Hartree[I1,I3] =
    #       accShift[g1-g3] *
    #       V(g1-g3) *
    #       F[I1,I3]
    # ------------------------------------------------------------

    ix = csr.ix
    accShift = work.HartreeAccShift

    @assert length(accShift) == ix.NSHIFT

    offsets = csr.offsets
    g2_list = csr.g2_list
    g4_list = csr.g3_list

    Threads.@threads :greedy for sid in 1:ix.NSHIFT
        acc = 0.0 + 0.0im

        lo = Int(offsets[sid])
        hi = Int(offsets[sid + 1]) - 1

        @inbounds for p in lo:hi
            g2 = Int(g2_list[p])
            g4 = Int(g4_list[p])

            for a2 in 1:nband, a4 in 1:nband
                I2 = state_index(g2, a2, nband)
                I4 = state_index(g4, a4, nband)

                acc += P[I4, I2] * overlapmatrix[I2, I4]
            end
        end

        accShift[sid] = acc
    end

    Threads.@threads :greedy for I1 in 1:dimension
        @inbounds begin
            g1 = state_g(I1, nband)

            w1n1 = wave_n1[g1]
            w1n2 = wave_n2[g1]

            for I3 in 1:I1
                g3 = state_g(I3, nband)

                sid = shift_id(
                    ix,
                    w1n1 - wave_n1[g3],
                    w1n2 - wave_n2[g3],
                )

                V13 = Coulomb_transfer[g1, g3]

                work.HartreeMatrix[I1, I3] =
                    accShift[sid] *
                    V13 *
                    overlapmatrix[I1, I3]
            end
        end
    end

    symmetrize_from_lower!(work.HartreeMatrix)
    symmetrize_from_lower!(work.FockMatrix)

    work.HartreeMatrix ./= Area
    work.FockMatrix ./= Area

    copy!(work.H_phys, single_Ham)
    work.H_phys .+= single_MoirePo
    work.H_phys .+= pinning_po
    work.H_phys .+= work.HartreeMatrix
    work.H_phys .-= work.FockMatrix

    hermitize!(work.H_phys)

    return nothing
end






# ============================================================
# Basis transport
#
# S maps old active basis to new active basis:
#
#   Q_new = S * Q_old * S'
#
# where Q = P + Pref is the physical density.
#
# If nband == 1:
#   use phase-only scalar overlap to avoid trace leakage.
#
# If kept_bands == [1,2]:
#   use the full 2x2 overlap block. This is unitary at fixed g
#   because we keep the complete two-band Hilbert space.
#
# Cutoff changes still cause rectangular S and physical loss/gain
# at the boundary, just like in the RMG TDHF code.
# ============================================================

function build_S_matrix(
    wave_new::Vector{Vector{Int64}},
    spinor_new::Matrix{Vector{ComplexF64}},
    wave_old::Vector{Vector{Int64}},
    spinor_old::Matrix{Vector{ComplexF64}},
    kept_bands::Vector{Int64},
)
    nband = length(kept_bands)

    dim_new = length(wave_new) * nband
    dim_old = length(wave_old) * nband

    S = zeros(ComplexF64, dim_new, dim_old)

    old_index = Dict{Tuple{Int64, Int64}, Int64}()

    for gold in eachindex(wave_old)
        old_index[(wave_old[gold][1], wave_old[gold][2])] = Int64(gold)
    end

    for gnew in eachindex(wave_new)
        gold = get(old_index, (wave_new[gnew][1], wave_new[gnew][2]), 0)

        if gold == 0
            continue
        end

        if nband == 1
            amp = spinor_new[gnew, 1]' * spinor_old[gold, 1]

            phase = abs(amp) > 1e-12 ? amp / abs(amp) : 1.0 + 0.0im

            Inew = state_index(gnew, 1, nband)
            Iold = state_index(gold, 1, nband)

            S[Inew, Iold] = phase

        elseif nband == 2
            # Full two-band subspace. Use the actual 2x2 overlap block.
            # No polar decomposition by default.
            for anew in 1:nband, aold in 1:nband
                Inew = state_index(gnew, anew, nband)
                Iold = state_index(gold, aold, nband)

                S[Inew, Iold] = spinor_new[gnew, anew]' * spinor_old[gold, aold]
            end

        else
            error("Unexpected nband=$(nband). For lambda-jellium this should be 1 or 2.")
        end
    end

    return S
end


function transport_relative_density(
    P_old::Matrix{ComplexF64},
    Pref_old::Matrix{ComplexF64},
    Pref_new::Matrix{ComplexF64},
    S::Matrix{ComplexF64},
)
    Q_old = P_old + Pref_old
    Q_new = S * Q_old * S'
    P_new = Q_new - Pref_new

    hermitize!(P_new)
    return P_new
end


function transport_diagnostics(S::Matrix{ComplexF64})
    dim_new, dim_old = size(S)

    right_error = norm(S' * S - Matrix{ComplexF64}(I, dim_old, dim_old))
    left_error  = norm(S * S' - Matrix{ComplexF64}(I, dim_new, dim_new))

    return (
        dim_new = dim_new,
        dim_old = dim_old,
        right_unitarity_error = right_error,
        left_unitarity_error = left_error,
    )
end




# ============================================================
# Fermi factors for bath relaxation
# ============================================================

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
    eps::Vector{Float64},
    target_particle_number::Float64,
    temp::Float64;
    maxiter::Int = 300,
)
    if temp <= 0.0
        fill!(fermifactor, 0.0)

        nfill = Int(round(target_particle_number))
        @assert abs(target_particle_number - nfill) < 1e-8
        @assert 0 <= nfill <= length(eps)

        for n in 1:nfill
            fermifactor[n] = 1.0
        end

        mu = nfill == 0 ? minimum(eps) - 1.0 :
             nfill == length(eps) ? maximum(eps) + 1.0 :
             0.5 * (eps[nfill] + eps[nfill + 1])

        return mu, sum(fermifactor)
    end

    lo = minimum(eps) - 100.0 * temp - 10.0
    hi = maximum(eps) + 100.0 * temp + 10.0

    mu = 0.5 * (lo + hi)
    particle_number = 0.0

    tol = target_particle_number == 0.0 ? 1e-9 : abs(1e-8 * target_particle_number)

    for _ in 1:maxiter
        mu = 0.5 * (lo + hi)
        particle_number = 0.0

        @inbounds for i in eachindex(eps)
            occ = fermi_occ((eps[i] - mu) / temp)
            fermifactor[i] = occ
            particle_number += occ
        end

        if abs(particle_number - target_particle_number) < tol
            return mu, particle_number
        elseif particle_number > target_particle_number
            hi = mu
        else
            lo = mu
        end
    end

    return mu, particle_number
end


# ============================================================
# Old-basis Hamiltonian + bath evolution
#
# Input:
#   Q = physical density matrix in the current basis
#
# Output:
#   Q_new = physical density after dt evolution under frozen H_HF
#
# In HF eigenbasis:
#   Q_mn(t+dt) = exp[-i(e_m-e_n)dt] exp[-gamma dt] Q_mn
# for offdiagonal entries.
#
# Diagonal entries relax to Fermi factors if gamma > 0.
# ============================================================

function evolve_physical_density_old_basis(
    Q::Matrix{ComplexF64},
    H_HF::Matrix{ComplexF64},
    dt::Float64;
    gamma::Float64,
    temp::Float64,
    target_particle_number::Float64,
)
    dim = size(Q, 1)

    @assert size(Q) == (dim, dim)
    @assert size(H_HF) == (dim, dim)

    F = eigen(Hermitian(H_HF))
    eps = real(F.values)
    U = F.vectors

    Qtilde = U' * Q * U

    fermifactor = zeros(Float64, dim)

    mu, particle_number_check = find_FL_iterative!(
        fermifactor,
        eps,
        target_particle_number,
        temp,
    )

    decay = exp(-gamma * dt)

    for m in 1:dim
        for n in 1:dim
            phase = exp(-im * (eps[m] - eps[n]) * dt)

            if m == n
                if gamma > 0.0
                    Qtilde[m, m] = fermifactor[m] + decay * (Qtilde[m, m] - fermifactor[m])
                else
                    Qtilde[m, m] = Qtilde[m, m]
                end
            else
                if gamma > 0.0
                    Qtilde[m, n] = decay * phase * Qtilde[m, n]
                else
                    Qtilde[m, n] = phase * Qtilde[m, n]
                end
            end
        end
    end

    Q_new = U * Qtilde * U'
    hermitize!(Q_new)

    return (
        Q = Q_new,
        eps = eps,
        HF_eigenvector = U,
        mu = mu,
        particle_number_check = particle_number_check,
    )
end




# ============================================================
# One Strang TDHF step
#
# Input:
#   P_old      relative density at Ashift_old
#   wave_old   active wave list at Ashift_old
#   book_old   bookkeeping for wave_old
#   data_old   projected data for wave_old, Ashift_old
#
# Output:
#   P_new      relative density at Ashift_new
#   wave_new, book_new, data_new
#
# Strang structure:
#   old basis -> half basis transport
#   frozen HF+bath evolution at half basis for dt
#   half basis -> new basis transport
# ============================================================

# ============================================================
# One Strang TDHF step, save-friendly ordering
#
# Input state is assumed to be a saved/checkpoint-consistent state:
#
#   P_old, Ashift_old, wave_old, book_old, data_old,
#   H_old = H_HF[P_old, Ashift_old]
#
# Step:
#   1. half evolve Q_old = P_old + Pref_old using H_old
#   2. transport full physical density to A_new basis
#   3. build temporary H_HF[P_trans, A_new]
#   4. half evolve using temporary new HF
#   5. rebuild final H_HF[P_new, A_new] for saving
#
# Output:
#   P_new, wave_new, book_new, data_new, Ashift_new,
#   H_HF[P_new,A_new], Hartree[P_new,A_new], Fock[P_new,A_new]
# ============================================================

function tdhf_one_step_strang!(
    work::TDHFWork,
    Prj::Matrix{ComplexF64},
    wave_work::Vector{Vector{Int64}},
    book,
    proj,
    Ashift::Vector{Float64},
    params::TDHFParams,
)
    kept_bands = params.kept_bands
    nband = length(kept_bands)

    dt = params.dt
    half_dt = 0.5 * dt

    T1 = params.T1
    T2 = params.T2
    b1T = params.b1T
    b2T = params.b2T
    gcutoff_work = params.gcutoff_work

    Ashift_new = Ashift + dt * params.dAshift_dt

    dim_old = length(wave_work) * nband

    @assert size(Prj) == (dim_old, dim_old)
    @assert size(work.H_phys) == (dim_old, dim_old)
    @assert size(proj.Pref) == (dim_old, dim_old)

    # ------------------------------------------------------------
    # 1. First half-step under already-built H_HF[P_n,A_n].
    # work.H_phys is assumed to have been built in the outer loop.
    # ------------------------------------------------------------

    Q_old = Prj + proj.Pref

    target_particle_number_old =
        Float64(params.filling) + real(tr(proj.Pref))

    evolved_old_half = evolve_physical_density_old_basis(
        Q_old,
        work.H_phys,
        half_dt;
        gamma = params.gamma,
        temp = params.temp,
        target_particle_number = target_particle_number_old,
    )

    Q_mid_oldbasis = evolved_old_half.Q
    P_mid_oldbasis = Q_mid_oldbasis - proj.Pref
    hermitize!(P_mid_oldbasis)

    # ------------------------------------------------------------
    # 2. Build basis/data at A_{n+1}.
    # ------------------------------------------------------------

    wave_new = wave_work
    book_new = book

    proj_new = refresh_projected_data(
        wave_new,
        book_new,
        Ashift_new,
        params,
    )

    dim_new = length(wave_new) * nband

    # ------------------------------------------------------------
    # 3. Transport physical density from old basis to new basis.
    #
    # P_trans = S*(P_mid_oldbasis + Pref_old)*S' - Pref_new
    # ------------------------------------------------------------

    S_old_to_new = build_S_matrix(
        wave_new,
        proj_new.spinor_set,
        wave_work,
        proj.spinor_set,
        kept_bands,
    )

    P_trans = transport_relative_density(
        P_mid_oldbasis,
        proj.Pref,
        proj_new.Pref,
        S_old_to_new,
    )

    @assert size(P_trans) == (dim_new, dim_new)

    # ------------------------------------------------------------
    # 4. Build temporary new HF from transported density.
    #
    # This overwrites work if the dimension is unchanged. If the cutoff
    # changes and dim_new != dim_old, create a fresh work object.
    # ------------------------------------------------------------

    if dim_new == dim_old && book_new.csr.ix.NSHIFT == length(work.HartreeAccShift)
        work_temp = work
    else
        work_temp = TDHFWork(dim_new, book_new.csr.ix.NSHIFT)
    end

    Build_HHF_lambda!(
        work_temp,
        book_new.csr,
        wave_new,
        book_new.wave_n1,
        book_new.wave_n2,
        P_trans,
        proj_new.single_Ham,
        proj_new.single_MoirePo,
        proj_new.pinning_po,
        proj_new.overlapmatrix,
        params.Area,
        book_new.Coulomb_transfer,
        nband,
    )

    # ------------------------------------------------------------
    # 5. Second half-step under temporary H_HF[P_trans,A_{n+1}].
    # ------------------------------------------------------------

    Q_trans = P_trans + proj_new.Pref

    target_particle_number_new =
        Float64(params.filling) + real(tr(proj_new.Pref))

    evolved_new_half = evolve_physical_density_old_basis(
        Q_trans,
        work_temp.H_phys,
        half_dt;
        gamma = params.gamma,
        temp = params.temp,
        target_particle_number = target_particle_number_new,
    )

    Q_new = evolved_new_half.Q
    Prj_new = Q_new - proj_new.Pref
    hermitize!(Prj_new)

    # ------------------------------------------------------------
    # 6. Diagnostics.
    #
    # Important: work_temp.H_phys is not necessarily H[Prj_new,A_new].
    # The outer loop will rebuild H_HF before saving the next checkpoint.
    # ------------------------------------------------------------

    step_transport_diagnostics = transport_diagnostics(S_old_to_new)

    diagnostics = (
        trace_P_old = real(tr(Prj)),
        trace_Q_old = real(tr(Q_old)),

        trace_P_mid_oldbasis = real(tr(P_mid_oldbasis)),
        trace_Q_mid_oldbasis = real(tr(Q_mid_oldbasis)),

        trace_P_trans = real(tr(P_trans)),
        trace_Q_trans = real(tr(Q_trans)),

        trace_P_new = real(tr(Prj_new)),
        trace_Q_new = real(tr(Q_new)),

        target_particle_number_old = target_particle_number_old,
        target_particle_number_new = target_particle_number_new,

        particle_number_check_old_half = evolved_old_half.particle_number_check,
        particle_number_check_new_half = evolved_new_half.particle_number_check,

        mu_old_half = evolved_old_half.mu,
        mu_new_half = evolved_new_half.mu,

        Ashift_old = copy(Ashift),
        Ashift_new = copy(Ashift_new),

        dim_old = dim_old,
        dim_new = dim_new,

        transport = step_transport_diagnostics,
    )

    return (
        Prj = Prj_new,
        Ashift = Ashift_new,
        wave_work = wave_new,
        book = book_new,
        proj = proj_new,
        step_transport_diagnostics = diagnostics,
    )
end

# ============================================================
# Save TDHF checkpoint
# Mirrors your RMG save_tdhf_file! convention.
# Assumes work.H_phys has already been built for the current Prj/Ashift.
# ============================================================

function save_tdhf_file!(
    checkpoint_dir::String,
    Prj::Matrix{ComplexF64},
    Ashift::Vector{Float64},
    wave_work::Vector{Vector{Int64}},
    step_index::Int,
    time_now::Float64,
    diagnostics_record,
    seed_file_path::String,
    work::TDHFWork,
    proj,
    params::TDHFParams,
    args::Vector{Float64},
)
    mkpath(checkpoint_dir)

    save_file_path = checkpoint_path_from_args(
        checkpoint_dir,
        args,
        step_index,
    )

    JLD2.jldsave(
        save_file_path;
        Prj = Prj,
        Ashift = Ashift,
        wave_work = wave_work,
        step_index = step_index,
        time_now = time_now,
        diagnostics_record = diagnostics_record,
        seed_file_path = seed_file_path,

        args = args,
        params = params,

        H_HF = copy(work.H_phys),
        HartreeMatrix = copy(work.HartreeMatrix),
        FockMatrix = copy(work.FockMatrix),

        spinor_set = proj.spinor_set,
        overlapmatrix = proj.overlapmatrix,
        single_Ham = proj.single_Ham,
        single_MoirePo = proj.single_MoirePo,
        pinning_po = proj.pinning_po,
        Pref = proj.Pref,
    )

    println("Saved TDHF file: ", save_file_path)
    return save_file_path
end


# ============================================================
# TDHF run-state loader from static seed
# Keep this simple first: start_step = 0 only.
# Resume can be added later in the same pattern as RMG.
# ============================================================

function load_tdhf_starting_point(
    seed_file_path::String,
    output_dir::String,
    checkpoint_dir::String,
    kept_bands_tdhf::Vector{Int64},
    args::Vector{Float64};
    start_step::Int,
    total_steps::Int,
    save_every::Int,
    refresh_every::Int,
    gcutoff_work::Float64,
    gamma::Float64,
    temp::Float64,
    dt::Float64,
    dAshift_dt::Vector{Float64},
    allow_seed_truncation::Bool = false,
)
    if !isfile(seed_file_path)
        error("Seed file does not exist: $(seed_file_path)")
    end

    seed = load(seed_file_path)
    static_args = Vector{Float64}(seed["arguments"])
    kept_bands_seed = kept_bands_from_saved_args(static_args)

    params = make_tdhf_params(
        static_args,
        kept_bands_tdhf;
        gcutoff_work = gcutoff_work,
        gamma = gamma,
        temp = temp,
        dt = dt,
        dAshift_dt = dAshift_dt,
    )

    checkpoint_file_path = checkpoint_path_from_args(
        checkpoint_dir,
        args,
        start_step,
    )

    if isfile(checkpoint_file_path)
        println("Loading checkpoint: ", checkpoint_file_path)

        ck = load(checkpoint_file_path)

        Prj = Matrix{ComplexF64}(ck["Prj"])
        Ashift = Vector{Float64}(ck["Ashift"])
        wave_work = Vector{Vector{Int64}}(ck["wave_work"])
        time_now = Float64(ck["time_now"])

        diagnostics_record =
            haskey(ck, "diagnostics_record") ? ck["diagnostics_record"] : Any[]

        book = build_basis_bookkeeping(
            wave_work,
            params.T1,
            params.T2,
            params.rs,
            params.constq,
        )

        proj = refresh_projected_data(
            wave_work,
            book,
            Ashift,
            params,
        )

        nband = length(params.kept_bands)
        dimension = length(wave_work) * nband
        work = TDHFWork(dimension, book.csr.ix.NSHIFT)

        println("Resumed TDHF state:")
        println("  output_dir       = ", output_dir)
        println("  checkpoint_dir   = ", checkpoint_dir)
        println("  start_step       = ", start_step)
        println("  total_steps      = ", total_steps)
        println("  save_every       = ", save_every)
        println("  refresh_every    = ", refresh_every)
        println("  time_now         = ", time_now)
        println("  dimension        = ", dimension)
        println("  Tr(P_rel)        = ", real(tr(Prj)))
        println("  Tr(Pref)         = ", real(tr(proj.Pref)))
        println("  Tr(Q)            = ", real(tr(Prj + proj.Pref)))

        return (
            Prj = Prj,
            Ashift = Ashift,
            wave_work = wave_work,
            book = book,
            proj = proj,
            work = work,
            params = params,

            current_step = start_step,
            total_steps = total_steps,
            save_every = save_every,
            refresh_every = refresh_every,
            time_now = time_now,
            diagnostics_record = diagnostics_record,

            seed_file_path = seed_file_path,
            output_dir = output_dir,
            checkpoint_dir = checkpoint_dir,
        )
    end

    if start_step != 0
        error(
            "Requested start_step=$(start_step), but checkpoint does not exist:\n" *
            checkpoint_file_path
        )
    end

    println("No checkpoint_step_0 found. Initializing from seed: ", seed_file_path)

    Ashift = [0.0, 0.0]
    time_now = 0.0
    diagnostics_record = Any[]

    wave_work = make_wave_A(
        params.T1,
        params.T2,
        params.b1T,
        params.b2T,
        params.gcutoff_work,
        Ashift,
    )

    book = build_basis_bookkeeping(
        wave_work,
        params.T1,
        params.T2,
        params.rs,
        params.constq,
    )

    proj = refresh_projected_data(
        wave_work,
        book,
        Ashift,
        params,
    )

    Prj = remap_seed_density_to_tdhf_basis(
        Matrix{ComplexF64}(seed["densitymatrix"]),
        Vector{Vector{Int64}}(seed["wave"]),
        kept_bands_seed,
        wave_work,
        kept_bands_tdhf;
        allow_seed_truncation = allow_seed_truncation,
    )

    nband = length(params.kept_bands)
    dimension = length(wave_work) * nband
    work = TDHFWork(dimension, book.csr.ix.NSHIFT)

    println("Initialized TDHF state from seed:")
    println("  output_dir       = ", output_dir)
    println("  checkpoint_dir   = ", checkpoint_dir)
    println("  kept_bands_seed  = ", kept_bands_seed)
    println("  kept_bands_tdhf  = ", kept_bands_tdhf)
    println("  start_step       = 0")
    println("  total_steps      = ", total_steps)
    println("  save_every       = ", save_every)
    println("  refresh_every    = ", refresh_every)
    println("  dt               = ", dt)
    println("  dAshift_dt       = ", dAshift_dt)
    println("  gamma            = ", gamma)
    println("  temp             = ", temp)
    println("  gcutoff_work     = ", gcutoff_work)
    println("  dimension        = ", dimension)
    println("  Tr(P_rel)        = ", real(tr(Prj)))
    println("  Tr(Pref)         = ", real(tr(proj.Pref)))
    println("  Tr(Q)            = ", real(tr(Prj + proj.Pref)))

    return (
        Prj = Prj,
        Ashift = Ashift,
        wave_work = wave_work,
        book = book,
        proj = proj,
        work = work,
        params = params,

        current_step = 0,
        total_steps = total_steps,
        save_every = save_every,
        refresh_every = refresh_every,
        time_now = time_now,
        diagnostics_record = diagnostics_record,

        seed_file_path = seed_file_path,
        output_dir = output_dir,
        checkpoint_dir = checkpoint_dir,
    )
end

# ============================================================
# Main run loop
# Directly mirrors your RMG convention.
# ============================================================

function run_lambda_tdhf_from_seed!(
    seed_file_path::String,
    kept_bands_tdhf::Vector{Int64};
    output_dir::String,
    checkpoint_dir::String,
    args::Vector{Float64},
    start_step::Int,
    total_fluxes::Real,
    steps_per_flux::Int,
    save_every::Int,
    refresh_every::Int,
    direction_angle_degree::Float64 = 0.0,
    gcutoff_work::Float64,
    gamma::Float64 = 0.0,
    temp::Float64 = 1e-3,
    dt::Float64,
    allow_seed_truncation::Bool = false,
)

 seed = load(seed_file_path)
    static_args = Vector{Float64}(seed["arguments"])

    enlarge_factor = Int(round(static_args[5]))
    T1, T2, _, _, _, _, _, _, _ = lambda_geometry(enlarge_factor)

    dAshift_dt, deltaA_step, one_flux_vector =
        dAshift_dt_from_steps_per_flux(
            T1,
            direction_angle_degree,
            steps_per_flux,
            dt,
        )

    total_steps = Int(round(total_fluxes * steps_per_flux))

    println("Flux insertion setup:")
    println("  direction_angle_degree = ", direction_angle_degree)
    println("  steps_per_flux         = ", steps_per_flux)
    println("  total_fluxes           = ", total_fluxes)
    println("  total_steps            = ", total_steps)
    println("  one_flux_vector        = ", one_flux_vector)
    println("  deltaA_step            = ", deltaA_step)
    println("  dAshift_dt             = ", dAshift_dt)



    state = load_tdhf_starting_point(
        seed_file_path,
        output_dir,
        checkpoint_dir,
        kept_bands_tdhf,
        args;
        start_step = start_step,
        total_steps = total_steps,
        save_every = save_every,
        refresh_every = refresh_every,
        gcutoff_work = gcutoff_work,
        gamma = gamma,
        temp = temp,
        dt = dt,
        dAshift_dt = dAshift_dt,
        allow_seed_truncation = allow_seed_truncation,
    )

    Prj = state.Prj
    Ashift = state.Ashift
    wave_work = state.wave_work
    book = state.book
    proj = state.proj
    work = state.work
    params = state.params

    time_now = state.time_now
    diagnostics_record = state.diagnostics_record

    current_step = state.current_step

    while current_step <= state.total_steps
        # Build H_HF[P_n,A_n] before saving/evolving.
        Build_HHF_lambda!(
            work,
            book.csr,
            wave_work,
            book.wave_n1,
            book.wave_n2,
            Prj,
            proj.single_Ham,
            proj.single_MoirePo,
            proj.pinning_po,
            proj.overlapmatrix,
            params.Area,
            book.Coulomb_transfer,
            length(params.kept_bands),
        )

        if current_step == 0 ||
           current_step % state.save_every == 0 ||
           current_step == state.total_steps

            save_tdhf_file!(
            state.checkpoint_dir,
            Prj,
            Ashift,
            wave_work,
            current_step,
            time_now,
            diagnostics_record,
            state.seed_file_path,
            work,
            proj,
            params,
            args,
        )

            save_charge_density!(
            state.checkpoint_dir,
            args,
            Prj,
            Ashift,
            wave_work,
            proj,
            params,
            current_step,
        )
        end

        if current_step == state.total_steps
            break
        end

        out = tdhf_one_step_strang!(
            work,
            Prj,
            wave_work,
            book,
            proj,
            Ashift,
            params,
        )

            

        Prj = out.Prj
        Ashift = out.Ashift
        wave_work = out.wave_work
        book = out.book
        proj = out.proj

        time_now += params.dt
        next_step = current_step + 1

        trace_change_refresh = 0.0
        refresh_diagnostics = nothing

        # ------------------------------------------------------------
        # Refresh cutoff basis after arriving at next_step.
        # This mirrors the RMG convention.
        # ------------------------------------------------------------
        if next_step % state.refresh_every == 0
            trace_before_refresh = real(tr(Prj))

            refreshed = refresh_cutoff_basis(
                Prj,
                wave_work,
                book,
                proj,
                Ashift,
                params,
            )

            Prj = refreshed.Prj
            wave_work = refreshed.wave_work
            book = refreshed.book
            proj = refreshed.proj
            work = refreshed.work

            trace_after_refresh = real(tr(Prj))
            trace_change_refresh = trace_after_refresh - trace_before_refresh
            refresh_diagnostics = refreshed.refresh_diagnostics
        else
            dimension = size(Prj, 1)

            if size(work.H_phys, 1) != dimension || length(work.HartreeAccShift) != book.csr.ix.NSHIFT
                work = TDHFWork(dimension, book.csr.ix.NSHIFT)
            end
        end

        current_step = next_step

        push!(
            diagnostics_record,
            (
                step = current_step,
                step_diagnostics = out.step_transport_diagnostics,
                trace_change_refresh = trace_change_refresh,
                refresh_diagnostics = refresh_diagnostics,
            ),
        )

            Q = Prj + proj.Pref

        eigs_Prj = eigvals(Hermitian(Prj))
        eigs_Q = eigvals(Hermitian(Q))

        println(
            "step = ", current_step,
            " time = ", time_now,
            " Tr(P_rel) = ", real(tr(Prj)),
            " Tr(Q) = ", real(tr(Q)),
            " eig(P_rel) min/max = ", minimum(eigs_Prj), " / ", maximum(eigs_Prj),
            " eig(Q) min/max = ", minimum(eigs_Q), " / ", maximum(eigs_Q),
            " refresh ΔTr(P) = ", trace_change_refresh,
            " Ashift = ", Ashift,
            " dim = ", size(Prj, 1),
        )
    end

    return nothing
end



function rotate_vector(v::Vector{Float64}, angle_degree::Float64)
    c = cosd(angle_degree)
    s = sind(angle_degree)

    return [
        c * v[1] - s * v[2],
        s * v[1] + c * v[2],
    ]
end


function dAshift_dt_from_steps_per_flux(
    T1::Vector{Float64},
    direction_angle_degree::Float64,
    steps_per_flux::Int,
    dt::Float64,
)
    @assert steps_per_flux > 0
    @assert dt > 0

    one_flux_vector = rotate_vector(T1, direction_angle_degree)

    deltaA_step = one_flux_vector / steps_per_flux
    dAshift_dt = deltaA_step / dt

    return dAshift_dt, deltaA_step, one_flux_vector
end


function complete_reference_for_new_momenta!(
    Q_new::Matrix{ComplexF64},
    Pref_new::Matrix{ComplexF64},
    wave_new::Vector{Vector{Int64}},
    wave_old::Vector{Vector{Int64}},
    kept_bands::Vector{Int64},
)
    nband = length(kept_bands)

    old_index = Dict{Tuple{Int64, Int64}, Bool}()

    for gold in eachindex(wave_old)
        old_index[(wave_old[gold][1], wave_old[gold][2])] = true
    end

    for gnew in eachindex(wave_new)
        key = (wave_new[gnew][1], wave_new[gnew][2])

        if !haskey(old_index, key)
            # This momentum point newly entered the working shell.
            # Initialize it to the reference density, so P_rel = 0 there.
            for a in 1:nband, b in 1:nband
                I = state_index(gnew, a, nband)
                J = state_index(gnew, b, nband)

                Q_new[I, J] = Pref_new[I, J]
            end
        end
    end

    return Q_new
end




# ============================================================
# Refresh cutoff shell at fixed Ashift
#
# This is the analogue of the RMG refresh_cutoff_basis.
#
# Important:
#   - Ashift is fixed.
#   - wave_work may change.
#   - common states are transported by S.
#   - newly entered states are initialized to Pref_new,
#     so their relative density is zero.
# ============================================================

function refresh_cutoff_basis(
    Prj::Matrix{ComplexF64},
    wave_old::Vector{Vector{Int64}},
    book_old,
    proj_old,
    Ashift::Vector{Float64},
    params::TDHFParams,
)
    kept_bands = params.kept_bands
    nband = length(kept_bands)

    wave_new = make_wave_A(
        params.T1,
        params.T2,
        params.b1T,
        params.b2T,
        params.gcutoff_work,
        Ashift,
    )

    book_new = build_basis_bookkeeping(
        wave_new,
        params.T1,
        params.T2,
        params.rs,
        params.constq,
    )

    proj_new = refresh_projected_data(
        wave_new,
        book_new,
        Ashift,
        params,
    )

    S = build_S_matrix(
        wave_new,
        proj_new.spinor_set,
        wave_old,
        proj_old.spinor_set,
        kept_bands,
    )

    Q_old = Prj + proj_old.Pref
    Q_new = S * Q_old * S'

    complete_reference_for_new_momenta!(
        Q_new,
        proj_new.Pref,
        wave_new,
        wave_old,
        kept_bands,
    )

    Prj_new = Q_new - proj_new.Pref
    hermitize!(Prj_new)

    dim_new = length(wave_new) * nband
    work_new = TDHFWork(dim_new, book_new.csr.ix.NSHIFT)

    refresh_diag = (
        dim_old = size(Prj, 1),
        dim_new = dim_new,

        Ng_old = length(wave_old),
        Ng_new = length(wave_new),

        trace_P_before = real(tr(Prj)),
        trace_Q_before = real(tr(Q_old)),

        trace_P_after = real(tr(Prj_new)),
        trace_Q_after = real(tr(Prj_new + proj_new.Pref)),

        trace_P_change = real(tr(Prj_new)) - real(tr(Prj)),
        trace_Q_change = real(tr(Prj_new + proj_new.Pref)) - real(tr(Q_old)),

        transport = transport_diagnostics(S),
    )

    return (
        Prj = Prj_new,
        wave_work = wave_new,
        book = book_new,
        proj = proj_new,
        work = work_new,
        refresh_diagnostics = refresh_diag,
    )
end



# ============================================================
# Production path helpers
# ============================================================
function output_dir_from_args(args::Vector{Float64})
    output_index = Int(round(args[15]))

    scratch_root = ENV["SCRATCH"]

    return joinpath(
        scratch_root,
        "jellium_TDHF",
        "data_output$(output_index)",
    )
end

function seed_dir_from_output(output_dir::String)
    return joinpath(output_dir, "seed")
end

function checkpoint_dir_from_output(output_dir::String)
    return joinpath(output_dir, "checkpoint")
end


function seed_filename_from_args(args::Vector{Float64})
    return string(
        args[1], "gcut",
        args[2], "lda",
        args[3], "ban1",
        args[4], "band2",
        args[5], "elg",
        args[6], "V2",
        args[7], "phi",
        args[8], "pco",
        args[9], "rs",
        args[10], "ddis",
        args[11], "dpos",
        args[12], "Cq",
        args[13], "nu",
        args[14], "Del",
        ".jld2",
    )
end


function find_seed_file_from_args(seed_dir::String, args::Vector{Float64})
    seed_name = seed_filename_from_args(args)
    seed_path = joinpath(seed_dir, seed_name)

    if !isfile(seed_path)
        error(
            "Seed file not found.\n" *
            "Expected: $(seed_path)\n" *
            "Available files in seed_dir:\n" *
            join(readdir(seed_dir), "\n")
        )
    end

    return seed_path
end


function checkpoint_filename_from_args(args::Vector{Float64})
    return string(
        args[1], "gcut",
        args[2], "lda",
        args[3], "ban1",
        args[4], "band2",
        args[5], "elg",
        args[6], "V2",
        args[7], "phi",
        args[8], "pco",
        args[9], "rs",
        args[10], "ddis",
        args[11], "dpos",
        args[12], "Cq",
        args[13], "nu",
        args[14], "Del",
        args[18], "SPF",
        args[20], "rfr",
        args[21], "gma",
        args[22], "tmp",
        args[23], "dt",
        args[24], "angle",
        args[25], "cutwork",
        ".jld2",
    )
end


function checkpoint_path_from_args(checkpoint_dir::String, args::Vector{Float64}, step_index::Int)
    base = checkpoint_filename_from_args(args)
    stem = replace(base, ".jld2" => "")
    return joinpath(checkpoint_dir, "$(stem)_step$(step_index).jld2")
end


# ============================================================
# ARGS convention
# ============================================================
#
# args[1]  = gcutoff
# args[2]  = λ
# args[3]  = keep_band_1
# args[4]  = keep_band_2
# args[5]  = enlarge_factor
# args[6]  = V2_scalar
# args[7]  = ϕ
# args[8]  = pin_coeff
# args[9]  = rs
# args[10] = dedis
# args[11] = defec_pos
# args[12] = constq
# args[13] = filling
# args[14] = Δ
#
# args[15] = output_index          # data_output$(output_index)
# args[16] = start_step
# args[17] = total_fluxes
# args[18] = steps_per_flux
# args[19] = save_every
# args[20] = refresh_every
# args[21] = gamma
# args[22] = temp
# args[23] = dt
# args[24] = direction_angle_degree
# args[25] = gcutoff_work



function run_tdhf_from_args!(args::Vector{Float64})
    if length(args) < 25
        error("Expected at least 25 Float64 args. See ARGS convention comment.")
    end

    output_dir = output_dir_from_args(args)
    seed_dir = seed_dir_from_output(output_dir)
    checkpoint_dir = checkpoint_dir_from_output(output_dir)

    mkpath(output_dir)
    mkpath(seed_dir)
    mkpath(checkpoint_dir)

    seed_file_path = find_seed_file_from_args(seed_dir, args)

    kept_bands_tdhf = kept_bands_from_flags(args[3], args[4])

    start_step = Int(round(args[16]))
    total_fluxes = args[17]
    steps_per_flux = Int(round(args[18]))
    save_every = Int(round(args[19]))
    refresh_every = Int(round(args[20]))

    gamma = args[21]
    temp = args[22]
    dt = args[23]
    direction_angle_degree = args[24]
    gcutoff_work = args[25]

    println("Production jellium TDHF:")
    println("  pwd            = ", pwd())
    println("  output_dir     = ", output_dir)
    println("  seed_dir       = ", seed_dir)
    println("  checkpoint_dir = ", checkpoint_dir)
    println("  seed_file_path = ", seed_file_path)
    println("  start_step     = ", start_step)
    println("  checkpoint stem = ", checkpoint_filename_from_args(args))

    run_lambda_tdhf_from_seed!(
        seed_file_path,
        kept_bands_tdhf;
        output_dir = output_dir,
        checkpoint_dir = checkpoint_dir,
        args = args,
        start_step = start_step,
        total_fluxes = total_fluxes,
        steps_per_flux = steps_per_flux,
        save_every = save_every,
        refresh_every = refresh_every,
        direction_angle_degree = direction_angle_degree,
        gcutoff_work = gcutoff_work,
        gamma = gamma,
        temp = temp,
        dt = dt,
    )

    return nothing
end

function charge_density_dir_from_checkpoint_dir(checkpoint_dir::String)
    return joinpath(checkpoint_dir, "CD_data")
end


function charge_density_path_from_args(
    checkpoint_dir::String,
    args::Vector{Float64},
    step_index::Int,
)
    cd_dir = charge_density_dir_from_checkpoint_dir(checkpoint_dir)
    mkpath(cd_dir)

    checkpoint_name = basename(
        checkpoint_path_from_args(checkpoint_dir, args, step_index)
    )

    return joinpath(cd_dir, "CD" * checkpoint_name)
end


function Densitymap_customize(
    xrange::Vector{Float64},
    yrange::Vector{Float64},
    DM::Matrix{ComplexF64},
    wave::Vector{Vector{Int64}},
    T1::Vector{Float64},
    T2::Vector{Float64},
)
    N = length(wave)

    coeff = Dict{Tuple{Int64, Int64}, ComplexF64}()

    @inbounds for j in 1:N, i in 1:N
        d = (
            wave[i][1] - wave[j][1],
            wave[i][2] - wave[j][2],
        )

        coeff[d] = get(coeff, d, 0.0 + 0.0im) + DM[i, j]
    end

    deltas = collect(keys(coeff))
    c = ComplexF64[coeff[d] for d in deltas]

    kvecs = Vector{SVector{2, Float64}}(undef, length(deltas))

    T1s = SVector{2, Float64}(T1[1], T1[2])
    T2s = SVector{2, Float64}(T2[1], T2[2])

    @inbounds for n in eachindex(deltas)
        d1, d2 = deltas[n]
        kvecs[n] = d1 * T1s + d2 * T2s
    end

    zvec = zeros(Float64, length(xrange), length(yrange))

    Threads.@threads for jb in eachindex(yrange)
        y = yrange[jb]

        @inbounds for ja in eachindex(xrange)
            r = SVector{2, Float64}(xrange[ja], y)
            s = 0.0 + 0.0im

            @simd for n in eachindex(c)
                s += c[n] * cis(dot(kvecs[n], r))
            end

            zvec[ja, jb] = real(s)
        end
    end

    return zvec
end


function save_charge_density!(
    checkpoint_dir::String,
    args::Vector{Float64},
    Prj::Matrix{ComplexF64},
    Ashift::Vector{Float64},
    wave_work::Vector{Vector{Int64}},
    proj,
    params::TDHFParams,
    step_index::Int,
)
    N = length(wave_work)
    nband = length(params.kept_bands)

    @assert size(Prj) == (N * nband, N * nband)

    densitymatrix_projected_rel = zeros(ComplexF64, N, N)

    norb = length(proj.spinor_set[1, 1])

    Threads.@threads for g1 in 1:N
        @inbounds for g2 in 1:N
            acc_rel = 0.0 + 0.0im

            for a1 in 1:nband, a2 in 1:nband
                I = state_index(g1, a1, nband)
                J = state_index(g2, a2, nband)

                ov = 0.0 + 0.0im

                for orb in 1:norb
                    ov += proj.spinor_set[g1, a1][orb] *
                          conj(proj.spinor_set[g2, a2][orb])
                end

                acc_rel += Prj[I, J] * ov
            end

            densitymatrix_projected_rel[g1, g2] = acc_rel
        end
    end

    xrange = collect(range(0.0, stop = 1.1 * norm(params.a1m), length = 160))
    yrange = collect(range(0.0, stop = 1.5 * norm(params.a1m), length = 160))

    zvec = Densitymap_customize(
        xrange,
        yrange,
        densitymatrix_projected_rel,
        wave_work,
        params.T1,
        params.T2,
    )

    save_file_path = charge_density_path_from_args(
        checkpoint_dir,
        args,
        step_index,
    )

    JLD2.jldsave(
        save_file_path;
        zvec = zvec,
        xrange = xrange,
        yrange = yrange,
        densitymatrix_projected_rel = densitymatrix_projected_rel,
        step_index = step_index,
        Ashift = Ashift,
        args = args,
        kept_bands = params.kept_bands,
    )

    println("Saved charge density: ", save_file_path)

    return save_file_path
end