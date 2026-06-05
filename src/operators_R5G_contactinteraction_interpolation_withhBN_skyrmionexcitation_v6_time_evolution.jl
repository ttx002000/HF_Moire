using LinearAlgebra
using Arpack
using Combinatorics
using Random
using JLD2
using StaticArrays
using Plots
BLAS.set_num_threads(1)


struct TDHFParams
    NL::Int
    uD::Float64
    λ::Float64
    enlarge_factor::Int

    V0_hBN::Float64
    V1_hBN::Float64
    ψ_hBN::Float64
    V2_scalar::Float64
    ϕ::Float64

    pin_coeff::Float64
    ϵr::Float64
    dedis::Float64
    defec_pos::Int
    constq::Float64

    filling::Int
    gamma::Float64
    temp::Float64
    dt::Float64
    gcutoff_work::Float64
    dAshift_dt::Vector{Float64}

    T1::Vector{Float64}
    T2::Vector{Float64}
    a1m::Vector{Float64}
    a2m::Vector{Float64}
    b1T::Vector{Int64}
    b2T::Vector{Int64}
    Area::Float64
end




mutable struct TDHFWork
    HartreeMatrix::Matrix{ComplexF64}
    FockMatrix::Matrix{ComplexF64}
    H_phys::Matrix{ComplexF64}
    HartreeAccShift::Vector{ComplexF64}
end

function TDHFWork(dimension::Int, num_shift_sectors::Int)
    zero_matrix = zeros(ComplexF64, dimension, dimension)

    return TDHFWork(
        copy(zero_matrix),
        copy(zero_matrix),
        copy(zero_matrix),
        zeros(ComplexF64, num_shift_sectors)
    )
end


# The goal of this piece of code is to do the calculation for specifically Nq=1 case(v4 is also for spinless).
function get_Ham(k::Vector{Float64},uD::Float64,valley::Int64,stacking::Int,NL::Int)
 
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




function get_Ham_Holomorphic(k::Vector{Float64},uD::Float64,valley::Int64,stacking::Int,NL::Int)
 
  Ham=zeros(ComplexF64,2*NL,2*NL)

  t0=3100
  t1=380
  fk=-√3/2*0.246*k[1]*valley+stacking*im*√3/2*0.246*k[2]

  for layer in 1:NL-1
     Ham[2*layer-1:2*layer,2*layer+1:2*layer+2]+=[0.0 0.0;t1 0.0]
  end



  Ham=Ham+Ham'

  for layer in 1:NL
      Ham[2*layer-1:2*layer,2*layer-1:2*layer]+=[0.0 -t0*fk;-t0*conj(fk) 0.0]
  end
  
   Ham[1:2,1:2]=[0.0 0.0; 0.0 0.0]

 for layer in 1:NL
   Ham[2*layer,2*layer]+=eigen(get_Ham(k,uD,valley,stacking,NL)).values[NL+1] 
 end

  for layer in 1:NL
   Ham[2*layer-1,2*layer-1]-=eigen(get_Ham(k,uD,valley,stacking,NL)).values[NL+1]  
 end
  
 return Ham
end




 

function get_f(k::Vector{Float64})
  delta1=1/√3*0.246*[0,1]
  delta2=1/√3*0.246*[√3/2,-1/2]
  delta3=1/√3*0.246*[-√3/2,-1/2]

  return exp(im*dot(k,delta1))+exp(im*dot(k,delta2))+exp(im*dot(k,delta3))
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



function get_Fourier_potential(k::Vector{Int},T1::Vector{Float64},T2::Vector{Float64},d::Float64)
   kvec=[T1 T2]*k
   D=25
    if k==[0,0]
        return (D-d)*9047.5646
    else
        return 1/cosh(norm(kvec)*D)*sinh(norm(kvec)*(D-d))/norm(kvec)*9047.5646
    end
end



@inline function symmetrize_from_lower!(A::Matrix{ComplexF64})
    n = size(A,1)
    @inbounds for i in 1:n
        A[i,i] = complex(real(A[i,i]), 0.0)
        for j in (i+1):n
            A[i,j] = conj(A[j,i])   # fill UPPER from LOWER
        end
    end
    return A
end



function Coulomb(k::Vector{Int},T1::Vector{Float64},T2::Vector{Float64})::Float64
   D=25
   return k==[0,0] ? D*9047.5636 : tanh(norm([T1 T2]*k*D))/norm(k[1]*T1+k[2]*T2)*9047.5636
end











function make_wave_A(T1, T2, b1T, b2T, gcutoff_work::Float64, Ashift::Vector{Float64})
    wave = Vector{Int64}[]
    cutoffstandard = gcutoff_work * norm(T1)
    cutoff = Int(ceil(gcutoff_work)) * 6

    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        w = ja*b1T + jb*b2T
        p = [T1 T2] * w + Ashift
        if dot(p, p) < cutoffstandard^2
            push!(wave, w)
        end
    end

    return wave
end



function build_basis_bookkeeping(
    wave::Vector{Vector{Int64}},
    T1::Vector{Float64},
    T2::Vector{Float64},
    ϵr::Float64,
    constq::Float64
)
    dimension = length(wave)

    wave_n1 = Vector{Int64}(undef, dimension)
    wave_n2 = Vector{Int64}(undef, dimension)

    @inbounds for g in 1:dimension
        wave_n1[g] = Int64(wave[g][1])
        wave_n2[g] = Int64(wave[g][2])
    end

    wave_dict = Dict{Vector{Int64}, Int64}()
    @inbounds for j in eachindex(wave)
        wave_dict[wave[j]] = Int64(j)
    end

    wave_diff = Vector{Int64}[]
    sizehint!(wave_diff, dimension^2)

    @inbounds for i in eachindex(wave), j in eachindex(wave)
        push!(wave_diff, wave[i] - wave[j])
    end
    wave_diff = unique(wave_diff)

    wl = build_wave_lookup(wave)
    csr = build_shiftcsr(wl, wave_n1, wave_n2)

    Coulomb_transfer = Matrix{Float64}(undef, dimension, dimension)

    @inbounds for i in eachindex(wave), j in eachindex(wave)
        Coulomb_transfer[i, j] = Coulomb(wave[i] - wave[j], T1, T2) / ϵr + constq
    end

    return (
        wave_n1 = wave_n1,
        wave_n2 = wave_n2,
        wave_dict = wave_dict,
        wave_diff = wave_diff,
        csr = csr,
        Coulomb_transfer = Coulomb_transfer
    )
end



function refresh_projected_data(
    wave::Vector{Vector{Int64}},
    book,
    Ashift::Vector{Float64},
    params::TDHFParams
)


    NL = params.NL
    uD = params.uD
    λ = params.λ
    enlarge_factor = params.enlarge_factor

    V0_hBN = params.V0_hBN
    V1_hBN = params.V1_hBN
    ψ_hBN = params.ψ_hBN
    V2_scalar = params.V2_scalar
    ϕ = params.ϕ

    pin_coeff = params.pin_coeff
    ϵr = params.ϵr
    dedis = params.dedis
    defec_pos = params.defec_pos

    T1 = params.T1
    T2 = params.T2
    a1m = params.a1m
    a2m = params.a2m
    b1T = params.b1T
    b2T = params.b2T
    Area = params.Area

    dimension = length(wave)

    spinor_set = Vector{Vector{ComplexF64}}(undef, dimension)
    single_Ham = zeros(ComplexF64, dimension, dimension)

    Threads.@threads :greedy for j in eachindex(wave)
        p = [T1 T2] * wave[j] + Ashift

        hh = (1.0 - λ) * get_Ham_Holomorphic(p, uD, 1, 1, NL) +
              λ        * get_Ham(p, uD, 1, 1, NL)

        F = eigen(Hermitian(hh))

        v1 = F.vectors[:, NL + 1]

        if uD > 0.0
            v1angle = angle(v1[2 * NL])
        else
            error("Current gauge fixing assumes uD > 0.")
        end

        spinor_set[j] = v1 * exp(-im * v1angle) / norm(v1)
        single_Ham[j, j] = real(F.values[NL + 1])
    end

    overlapmatrix = zeros(ComplexF64, dimension, dimension)

    Threads.@threads :greedy for i in eachindex(wave)
        @inbounds for j in eachindex(wave)
            overlapmatrix[i, j] = spinor_set[i]' * spinor_set[j]
        end
    end

    ω = exp(im * 2π / 3)

    op_1 = zeros(ComplexF64, 2 * NL, 2 * NL)
    op_1[1:2, 1:2] = [1 1; ω ω]

    op_2 = zeros(ComplexF64, 2 * NL, 2 * NL)
    op_2[1:2, 1:2] = [1 ω^2; ω^2 ω]

    op_3 = zeros(ComplexF64, 2 * NL, 2 * NL)
    op_3[1:2, 1:2] = [1 ω; 1 ω]

    op_4 = zeros(ComplexF64, 2 * NL, 2 * NL)
    op_4[1:2, 1:2] = [1 0; 0 1]

    op_5 = Matrix{ComplexF64}(I, 2 * NL, 2 * NL)

    single_MoirePo = zeros(ComplexF64, dimension, dimension)

   Threads.@threads :greedy for j in eachindex(wave)
        pos = get(book.wave_dict, wave[j] - enlarge_factor * b1T, 0)
        if pos != 0
            single_MoirePo[j, pos] += V1_hBN * exp(-im * ψ_hBN) *
                (spinor_set[j]' * op_1 * spinor_set[pos])

            single_MoirePo[j, pos] += V2_scalar * exp(-im * ϕ) *
                (spinor_set[j]' * op_5 * spinor_set[pos])
        end

        pos = get(book.wave_dict, wave[j] - enlarge_factor * b2T, 0)
        if pos != 0
            single_MoirePo[j, pos] += V1_hBN * exp(-im * ψ_hBN) *
                (spinor_set[j]' * op_2 * spinor_set[pos])

            single_MoirePo[j, pos] += V2_scalar * exp(-im * ϕ) *
                (spinor_set[j]' * op_5 * spinor_set[pos])
        end

        pos = get(book.wave_dict, wave[j] + enlarge_factor * (b1T + b2T), 0)
        if pos != 0
            single_MoirePo[j, pos] += V1_hBN * exp(-im * ψ_hBN) *
                (spinor_set[j]' * op_3 * spinor_set[pos])

            single_MoirePo[j, pos] += V2_scalar * exp(-im * ϕ) *
                (spinor_set[j]' * op_5 * spinor_set[pos])
        end

        single_MoirePo[j, j] += V0_hBN / 2 *
            (spinor_set[j]' * op_4 * spinor_set[j])
    end

    single_MoirePo = single_MoirePo + single_MoirePo'

    pinning_po = zeros(ComplexF64, dimension, dimension)

    if defec_pos == 1
        defec_pos_vec = (2 * a2m - a1m) / enlarge_factor * 2 / 3
    elseif defec_pos == 2
        defec_pos_vec = (2 * a2m - a1m) / enlarge_factor * 1 / 3
    elseif defec_pos == 3
        defec_pos_vec = (2 * a2m - a1m) / enlarge_factor * 0
    elseif defec_pos == 4
        defec_pos_vec = (2 * a2m - a1m) / enlarge_factor * (-2 / 3)
    else
        error("Unknown defec_pos = $defec_pos")
    end

    Threads.@threads :greedy for j in eachindex(wave)
        for d in eachindex(book.wave_diff)
            pos = get(book.wave_dict, wave[j] + book.wave_diff[d], 0)

            if pos != 0
                qvec = [T1 T2] * book.wave_diff[d]
                phase = exp(im * dot(qvec, defec_pos_vec))

                pinning_po[pos, j] += pin_coeff / (Area * ϵr) *
                    get_Fourier_potential(book.wave_diff[d], T1, T2, dedis) *
                    (spinor_set[pos]' * op_5 * spinor_set[j]) *
                    phase
            end
        end
    end

    pinning_po = (pinning_po + pinning_po') / 2

    Coulomb_matrix = book.Coulomb_transfer .* overlapmatrix

    return (
        spinor_set = spinor_set,
        overlapmatrix = overlapmatrix,
        single_Ham = single_Ham,
        single_MoirePo = single_MoirePo,
        pinning_po = pinning_po,
        Coulomb_matrix = Coulomb_matrix
    )
end




function Build_HHF!(work::TDHFWork,csr::ShiftCSR,wave::Vector{Vector{Int64}},wave_n1::Vector{Int64},wave_n2::Vector{Int64},
                               input_DensityMatrix::Matrix{ComplexF64},single_Ham::Matrix{ComplexF64},
                               single_MoirePo::Matrix{ComplexF64},pinning_po::Matrix{ComplexF64},overlapmatrix::Matrix{ComplexF64},
                              Area::Float64,Coulomb_matrix::Matrix{ComplexF64})
  
 
   dimension=length(wave)

 

    fill!(work.HartreeMatrix, 0)
    fill!(work.FockMatrix, 0)
  


  


  Threads.@threads  :greedy  for g1 in eachindex(wave)
    @inbounds begin
    
        w1n1=wave_n1[g1]; 
        w1n2=wave_n2[g1]
      for g4 in 1:g1
          
            w4n1=wave_n1[g4]; 
            w4n2=wave_n2[g4]
            acc = 0.0 + 0.0im
         
            dg_n1=w1n1- w4n1
            dg_n2=w1n2- w4n2
            lo, hi = csr_range(csr,dg_n1, dg_n2)
          
            
              @inbounds for k in lo:hi
                 g2 = Int(csr.g2_list[k])
                 g3 = Int(csr.g3_list[k])
              
              
                  tmp=Coulomb_matrix[g1,g3]*overlapmatrix[g2,g4]
                  acc += input_DensityMatrix[g3,g2] * tmp
                 
              end
                
          

              work.FockMatrix[g1,g4] = acc
        
              
          
  

      end
    end
  end


    
    ix = csr.ix
    accShift = work.HartreeAccShift

    offsets = csr.offsets
    g2_list  = csr.g2_list
    g3_list  = csr.g3_list   # in Hartree interpretation this is (g4)

    Threads.@threads :greedy for sid in 1:ix.NSHIFT
        acc = 0.0 + 0.0im
        lo = Int(offsets[sid])
        hi = Int(offsets[sid + 1]) - 1

        @inbounds for k in lo:hi
            g2 = Int(g2_list[k])
            g4 = Int(g3_list[k])
            acc += input_DensityMatrix[g4, g2] * overlapmatrix[g2, g4]
        end

        accShift[sid] = acc
    end
    Threads.@threads :greedy for g1 in eachindex(wave_n1)
    @inbounds begin
        w1n1 = Int(wave_n1[g1]); w1n2 = Int(wave_n2[g1])
        for g3 in 1:g1
            sid = shift_id(ix,
                w1n1 - Int(wave_n1[g3]),
                w1n2 - Int(wave_n2[g3])
            )
            work.HartreeMatrix[g1, g3] = accShift[sid] * Coulomb_matrix[g1, g3]
        end
    end
   end


   symmetrize_from_lower!(work.HartreeMatrix)
    symmetrize_from_lower!(work.FockMatrix)
    work.HartreeMatrix ./= Area
    work.FockMatrix    ./= Area

    copy!(work.H_phys, single_MoirePo)
    work.H_phys .+= single_Ham
    work.H_phys .+= pinning_po
    work.H_phys .+= work.HartreeMatrix
    work.H_phys .-= work.FockMatrix

     work.H_phys .= (work.H_phys .+ work.H_phys') ./ 2
   

     
 


 return nothing
end






function build_S_matrix(
    wave_new::Vector{Vector{Int64}},
    spinor_new::Vector{Vector{ComplexF64}},
    Ashift_new::Vector{Float64},
    wave_old::Vector{Vector{Int64}},
    spinor_old::Vector{Vector{ComplexF64}},
    Ashift_old::Vector{Float64},
    T1::Vector{Float64},
    T2::Vector{Float64},
    L1::Vector{Float64},
    L2::Vector{Float64}
)
    S = zeros(ComplexF64, length(wave_new), length(wave_old))

    old_index = Dict{Tuple{Int64, Int64}, Int64}()

    # Keep this serial: Dict writes are not thread-safe.
    @inbounds for old_momentum_index in eachindex(wave_old)
        old_label = (
            Int64(wave_old[old_momentum_index][1]),
            Int64(wave_old[old_momentum_index][2])
        )
        old_index[old_label] = Int64(old_momentum_index)
    end

    # Safe to thread: each thread writes to a different row of S,
    # and old_index is read-only here.
    Threads.@threads :greedy for new_momentum_index in eachindex(wave_new)
        @inbounds begin
            new_label = (
                Int64(wave_new[new_momentum_index][1]),
                Int64(wave_new[new_momentum_index][2])
            )

            old_momentum_index = get(old_index, new_label, Int64(0))

            if old_momentum_index != 0
                amp = spinor_new[new_momentum_index]' * spinor_old[old_momentum_index]

                if abs(amp) > 1e-12
                    S[new_momentum_index, old_momentum_index] = amp / abs(amp)
                else
                    S[new_momentum_index, old_momentum_index] = 1.0 + 0.0im
                end
            end
        end
    end

    return S
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
    target_particle_number::Float64,
    val_s::Float64,
    val_e::Float64,
    temp::Float64;
    maxiter::Int = 300,
)
    lo = val_s
    hi = val_e
    try_FL = 0.5 * (lo + hi)

    stan = target_particle_number == 0.0 ? 1e-9 : abs(1e-8 * target_particle_number)

    particle_number = 0.0

    for _ in 1:maxiter
        try_FL = 0.5 * (lo + hi)

        particle_number = 0.0

        @inbounds for i in eachindex(quasi_particle_energy)
            occ = fermi_occ((quasi_particle_energy[i] - try_FL) / temp)
            fermifactor[i] = occ
            particle_number += occ
        end

        if abs(particle_number - target_particle_number) < stan
            return try_FL, particle_number
        elseif particle_number > target_particle_number
            hi = try_FL
        else
            lo = try_FL
        end
    end

    return try_FL, particle_number
end




function evolve_old_basis_dissipative_step(
    Prj::Matrix{ComplexF64},
    H_HF::Matrix{ComplexF64},
    dt::Float64;
    gamma::Float64,
    temp::Float64,
    filling::Int
)
    bath_diagonalization = eigen(Hermitian(H_HF))
    eps = real(bath_diagonalization.values)
    bath_eigenvectors = bath_diagonalization.vectors

    # temporary change of coordinates to instantaneous HF eigenbasis
    Prj_tilde = bath_eigenvectors' * Prj * bath_eigenvectors

    # find bath chemical potential so that sum_n f_n = filling
    fermifactor = zeros(Float64, length(eps))

    val_s = minimum(eps) - 100.0 * temp - 10.0
    val_e = maximum(eps) + 100.0 * temp + 10.0

    mu, particle_number_check = find_FL_iterative!(
        fermifactor,
        eps,
        Float64(filling),
        val_s,
        val_e,
        temp
    )

    decay = exp(-gamma * dt)

        Threads.@threads :greedy for m in eachindex(eps)
            @inbounds for n in eachindex(eps)
                if m == n
                    Prj_tilde[m, n] = fermifactor[m] + decay * (Prj_tilde[m, n] - fermifactor[m])
                else
                    phase = exp(-im * (eps[m] - eps[n]) * dt)
                    Prj_tilde[m, n] = decay * phase * Prj_tilde[m, n]
                end
            end
        end

    # rotate back to the old projected-band basis
    Prj_star = bath_eigenvectors * Prj_tilde * bath_eigenvectors'
    Prj_star = (Prj_star + Prj_star') / 2

    return Prj_star, eps, bath_eigenvectors, mu, particle_number_check
end



function tdhf_one_step!(
    work::TDHFWork,
    Prj::Matrix{ComplexF64},
    wave_work::Vector{Vector{Int64}},
    book,
    proj,
    Ashift::Vector{Float64},
    params::TDHFParams
)
    T1 = params.T1
    T2 = params.T2
    a1m = params.a1m
    a2m = params.a2m
    Area = params.Area

    dt = params.dt
    gamma = params.gamma
    temp = params.temp
    filling = params.filling
    dAshift_dt = params.dAshift_dt




  
    


    # 2. Physical TDHF + reservoir step in the old A-basis
    Prj_star, eps, bath_eigenvectors, chemical_potential, particle_number_check  =
        evolve_old_basis_dissipative_step(
            Prj,
            work.H_phys,
            dt;
            gamma = gamma,
            temp = temp,
            filling = filling
        )

    # 3. Update vector potential shift
    Ashift_new = Ashift .+ dAshift_dt .* dt

    # 4. Rebuild A-dependent projected data at the new Ashift
    proj_new = refresh_projected_data(
        wave_work,
        book,
        Ashift_new,
        params
    )

    # 5. Build moving-basis overlap S_{new,old}
    Smatrix = build_S_matrix(
        wave_work,
        proj_new.spinor_set,
        Ashift_new,
        wave_work,
        proj.spinor_set,
        Ashift,
        T1,
        T2,
        a1m,
        a2m
    )

        # 6. Transport density matrix to new instantaneous basis.
        # We use the unitary/isometric part W of S for the actual projected-band dynamics.
        # The raw nonunitarity of S is kept as a diagnostic.

                trace_before_old_basis_step = real(tr(Prj))
        trace_after_old_basis_step = real(tr(Prj_star))

        old_blas_threads = BLAS.get_num_threads()
        BLAS.set_num_threads(min(Threads.nthreads(), 16))

        Prj_new = Smatrix * Prj_star * Smatrix'

        BLAS.set_num_threads(old_blas_threads)

        Prj_new = (Prj_new + Prj_new') / 2

        trace_after_transport = real(tr(Prj_new))

        step_transport_diagnostics = (
            trace_change_old_basis = trace_after_old_basis_step - trace_before_old_basis_step,
            transport_trace_change = trace_after_transport - trace_after_old_basis_step
        )

  return Prj_new, proj_new, Ashift_new, eps, bath_eigenvectors,
    chemical_potential, particle_number_check, step_transport_diagnostics
end



function refresh_cutoff_basis(
    Prj::Matrix{ComplexF64},
    wave_old::Vector{Vector{Int64}},
    book_old,
    proj_old,
    Ashift::Vector{Float64},
    params::TDHFParams
)


    
    T1 = params.T1
    T2 = params.T2
    a1m = params.a1m
    a2m = params.a2m
    b1T = params.b1T
    b2T = params.b2T

    ϵr = params.ϵr
    constq = params.constq
    gcutoff_work = params.gcutoff_work


    wave_new = make_wave_A(T1, T2, b1T, b2T, gcutoff_work, Ashift)

    book_new = build_basis_bookkeeping(
        wave_new,
        T1,
        T2,
        ϵr,
        constq
    )

    proj_new = refresh_projected_data(
        wave_new,
        book_new,
        Ashift,
        params
    )

    cutoff_overlap = build_S_matrix(
        wave_new,
        proj_new.spinor_set,
        Ashift,
        wave_old,
        proj_old.spinor_set,
        Ashift,
        T1,
        T2,
        a1m,
        a2m
    )

    Prj_new = cutoff_overlap * Prj * cutoff_overlap'
    Prj_new = (Prj_new + Prj_new') / 2

    work_new = TDHFWork(length(wave_new), book_new.csr.ix.NSHIFT)

    return Prj_new, wave_new, book_new, proj_new, work_new
end




function seed_filename_from_args(args)
    return "$(args[1])NL$(args[2])theta$(args[3])constq$(args[4])ϵr$(args[5])uD$(args[6])filling$(args[7])cutoff$(args[8])lambda$(args[9])trytime$(args[10])enlarge$(args[11])V0_hBN$(args[12])V1_hBN$(args[13])ψ_hBN$(args[14])V2_scalar$(args[15])ϕ$(args[16])pincof$(args[17])dedis$(args[18])depos.jld2"
end

function tdhf_filename_for_step(args, stepnum::Int)
    return "$(args[1])NL$(args[2])theta$(args[3])constq$(args[4])ϵr$(args[5])uD$(args[6])filling$(args[7])cutoff$(args[8])lambda$(args[9])trytime$(args[10])enlarge$(args[11])V0_hBN$(args[12])V1_hBN$(args[13])ψ_hBN$(args[14])V2_scalar$(args[15])ϕ$(args[16])pincof$(args[17])dedis$(args[18])depos" *
           "$(args[22])refreshevery" *
           "$(args[24])dt$(args[25])Emag$(args[26])Eag$(args[27])gamma$(args[28])temp$(args[29])workcutoff" *
           "$(stepnum)stepnum.jld2"
end

function hf_output_dir_from_args(args)
    scratch_dir = ENV["SCRATCH"]
    filepos = Int(args[19])

    return joinpath(
        scratch_dir,
        "triangle_R5G_contact_interpolation_withhBN_skyrmionexcitation_v6_time_evolution",
        "data_output$(filepos)"
    )
end

function seed_file_path_from_args(args)
    hf_output_dir = hf_output_dir_from_args(args)

    return joinpath(
        hf_output_dir,
        "seed",
        seed_filename_from_args(args)
    )
end

function tdhf_file_path_for_step(args, stepnum::Int)
    hf_output_dir = hf_output_dir_from_args(args)

    return joinpath(
        hf_output_dir,
        tdhf_filename_for_step(args, stepnum)
    )
end


function electric_field_from_magnitude_angle(
    electric_field_magnitude::Float64,
    electric_field_angle_degree::Float64,
    a1m::Vector{Float64}
)
    aag = electric_field_angle_degree / 180.0 * π

    xhat = a1m / norm(a1m)

    # +90 degree rotation of xhat in the current Cartesian coordinate system
    yhat = [-xhat[2], xhat[1]]

    return electric_field_magnitude * (cos(aag) .* xhat .+ sin(aag) .* yhat)
end





function load_tdhf_starting_point(args)
    # -------------------------
    # Fixed model parameters
    # -------------------------
    NL = Int(args[1])
    constq = args[3]
    ϵr = args[4]
    uD = args[5]
    filling = Int(args[6])
    λ = args[8]
    enlarge_factor = Int(args[10])

    V0_hBN = args[11]
    V1_hBN = args[12]
    ψ_hBN = args[13]
    V2_scalar = args[14]
    ϕ = args[15] / 180 * π

    pin_coeff = args[16]
    dedis = args[17]
    defec_pos = Int(args[18])

    # -------------------------
    # TDHF parameters
    # -------------------------
    total_steps = Int(args[20])
    save_every = Int(args[21])
    refresh_every = Int(args[22])
    plot_every = Int(args[23])

    dt = args[24]

    electric_field_magnitude = args[25]
    electric_field_angle_degree = args[26]

    gamma = args[27]
    temp = args[28]
    gcutoff_work = args[29]

    start_step = Int(args[30])


    # -------------------------
    # Load from seed or checkpoint
    # -------------------------
    if start_step == 0
        seed_file_path = seed_file_path_from_args(args)

        if !isfile(seed_file_path)
            error("Seed file does not exist: $seed_file_path")
        end

        println("Loading seed: ", seed_file_path)
        seed = load(seed_file_path)

        Ashift = [0.0, 0.0]
        time_now = 0.0
        diagnostics_record = NamedTuple[]

        P_seed = seed["densitymatrix"]
        wave_seed = seed["wave"]

        T1 = seed["T1"]
        T2 = seed["T2"]
        b1 = seed["b1"]
        b2 = seed["b2"]

        # b1, b2 are physical reciprocal vectors.
        # b1T, b2T are their integer coordinates in the T1,T2 basis.
        b1T = Int64.(round.(inv([T1 T2]) * b1))
        b2T = Int64.(round.(inv([T1 T2]) * b2))

        wave_work = make_wave_A(
            T1,
            T2,
            b1T,
            b2T,
            gcutoff_work,
            Ashift
        )

        old_index = Dict{Vector{Int64}, Int64}()
        for wave_index in eachindex(wave_seed)
            old_index[wave_seed[wave_index]] = wave_index
        end

        Prj = zeros(ComplexF64, length(wave_work), length(wave_work))

        for new_row in eachindex(wave_work), new_col in eachindex(wave_work)
            old_row = get(old_index, wave_work[new_row], 0)
            old_col = get(old_index, wave_work[new_col], 0)

            if old_row != 0 && old_col != 0
                Prj[new_row, new_col] = P_seed[old_row, old_col]
            end
        end

        Prj = (Prj + Prj') / 2

        

    else
        checkpoint_path = tdhf_file_path_for_step(args, start_step)

        if !isfile(checkpoint_path)
            error("TDHF checkpoint does not exist: $checkpoint_path")
        end

        println("Loading TDHF checkpoint: ", checkpoint_path)
        checkpoint = load(checkpoint_path)

        Prj = checkpoint["Prj"]
        Ashift = checkpoint["Ashift"]
        wave_work = checkpoint["wave_work"]
        time_now = checkpoint["time_now"]
        diagnostics_record = checkpoint["diagnostics_record"]

        seed_file_path = checkpoint["seed_file_path"]

        if !isfile(seed_file_path)
            error("Original seed file does not exist: $seed_file_path")
        end

        seed = load(seed_file_path)

        
    end

    # -------------------------
    # Reconstruct geometry and params
    # -------------------------
    T1 = seed["T1"]
    T2 = seed["T2"]
    a1m = seed["a1m"]
    a2m = seed["a2m"]
    b1 = seed["b1"]
    b2 = seed["b2"]

    b1T = Int64.(round.(inv([T1 T2]) * b1))
    b2T = Int64.(round.(inv([T1 T2]) * b2))

    Area = abs(a1m[1] * a2m[2] - a1m[2] * a2m[1])


            Efield = electric_field_from_magnitude_angle(
            electric_field_magnitude,
            electric_field_angle_degree,
            a1m
        )

        deltaA_step = -Efield .* dt
        dAshift_dt = deltaA_step ./ dt

        tdhf_args = Float64.([
            total_steps,
            save_every,
            refresh_every,
            plot_every,
            dt,
            electric_field_magnitude,
            electric_field_angle_degree,
            Efield[1],
            Efield[2],
            gamma,
            temp,
            gcutoff_work,
            start_step
        ])






    params = TDHFParams(
        Int(NL),
        uD,
        λ,
        enlarge_factor,

        V0_hBN,
        V1_hBN,
        ψ_hBN,
        V2_scalar,
        ϕ,

        pin_coeff,
        ϵr,
        dedis,
        defec_pos,
        constq,

        filling,
        gamma,
        temp,
        dt,
        gcutoff_work,
        dAshift_dt,

        T1,
        T2,
        a1m,
        a2m,
        b1T,
        b2T,
        Area
    )

    book = build_basis_bookkeeping(
        wave_work,
        params.T1,
        params.T2,
        params.ϵr,
        params.constq
    )

    proj = refresh_projected_data(
        wave_work,
        book,
        Ashift,
        params
    )

    work = TDHFWork(length(wave_work), book.csr.ix.NSHIFT)

    println("TDHF starting point:")
    println("  start_step = ", start_step)
    println("  total_steps = ", total_steps)
    println("  save_every = ", save_every)
    println("  refresh_every = ", refresh_every)
    println("  plot_every = ", plot_every)
    println("  time_now = ", time_now)
    println("  Ashift = ", Ashift)
    println("  electric_field_magnitude = ", electric_field_magnitude)
    println("  electric_field_angle_degree = ", electric_field_angle_degree)
    println("  Efield Cartesian = ", Efield)
    println("  deltaA_step = ", deltaA_step)
    println("  dAshift_dt = ", dAshift_dt)
    println("  gamma = ", gamma)
    println("  temp = ", temp)
    println("  gcutoff_work = ", gcutoff_work)
    println("  dimension = ", length(wave_work))
    println("  Tr(P) = ", real(tr(Prj)))

    return (
        Prj = Prj,
        Ashift = Ashift,
        wave_work = wave_work,
        book = book,
        proj = proj,
        work = work,
        params = params,
        total_steps = total_steps,
        save_every = save_every,
        refresh_every = refresh_every,
        plot_every = plot_every,
        time_now = time_now,

        diagnostics_record = diagnostics_record,
        seed_file_path = seed_file_path,

        args = args,
        tdhf_args = tdhf_args,
        deltaA_step = deltaA_step,
        dAshift_dt = dAshift_dt,
        Efield = Efield
    )
end


function save_tdhf_file!(
    args,
    Prj,
    Ashift,
    wave_work,
    step_index,
    time_now,
    diagnostics_record,
    tdhf_args,
    seed_file_path,
    deltaA_step,
    dAshift_dt,
    Efield,
    work,
    proj
)
    save_file_path = tdhf_file_path_for_step(args, step_index)

    mkpath(dirname(save_file_path))
    
    JLD2.jldsave(
        save_file_path;
        Prj = Prj,
        Ashift = Ashift,
        wave_work = wave_work,
        step_index = step_index,
        time_now = time_now,
        diagnostics_record = diagnostics_record,
        args = args,
        tdhf_args = tdhf_args,
        seed_file_path = seed_file_path,
        deltaA_step = deltaA_step,
        dAshift_dt = dAshift_dt,
        Efield = Efield,

        H_HF = copy(work.H_phys),
        HartreeMatrix = copy(work.HartreeMatrix),
        FockMatrix = copy(work.FockMatrix),
        spinor_set = proj.spinor_set
    )

    println("Saved TDHF file: ", save_file_path)

    return save_file_path
end






function run_tdhf_from_args!(args)
    state = load_tdhf_starting_point(args)

    Prj = state.Prj
    Ashift = state.Ashift
    wave_work = state.wave_work
    book = state.book
    proj = state.proj
    work = state.work
    params = state.params

    time_now = state.time_now
    diagnostics_record = state.diagnostics_record

    current_step = Int(args[30])

    while current_step <= state.total_steps

        # ------------------------------------------------------------
        # Current convention:
        # current_step = n means:
        #   Prj      = P_n
        #   Ashift   = A_n
        #   time_now = n * dt
        #
        # Build H_HF[P_n, A_n] before saving/evolving.
        # ------------------------------------------------------------
        Build_HHF!(
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
            proj.Coulomb_matrix
        )

        # ------------------------------------------------------------
        # Save checkpoint at the current step, before evolution.
        # This includes step 0 naturally.
        # ------------------------------------------------------------
        if current_step == 0 ||
           current_step % state.save_every == 0 ||
           current_step == state.total_steps

            save_tdhf_file!(
                state.args,
                Prj,
                Ashift,
                wave_work,
                current_step,
                time_now,
                diagnostics_record,
                state.tdhf_args,
                state.seed_file_path,
                state.deltaA_step,
                state.dAshift_dt,
                state.Efield,
                work,
                proj
            )
        end

        # ------------------------------------------------------------
        # Save charge density at the current step, before evolution.
        # This also includes step 0 naturally.
        # ------------------------------------------------------------
        if current_step == 0 ||
           current_step % state.plot_every == 0 ||
           current_step == state.total_steps

            save_charge_density_all_layers!(
                state.args,
                Prj,
                wave_work,
                proj,
                params,
                current_step
            )
        end

        # If this is the final requested state, stop after saving it.
        if current_step == state.total_steps
            break
        end

        # ------------------------------------------------------------
        # Evolve current_step -> current_step + 1
        # tdhf_one_step! should use the already-built work.H_phys.
        # ------------------------------------------------------------
        trace_before_step = real(tr(Prj))

        Prj, proj, Ashift, eps, bath_eigenvectors,
        chemical_potential, particle_number_check, step_transport_diagnostics =
            tdhf_one_step!(
                work,
                Prj,
                wave_work,
                book,
                proj,
                Ashift,
                params
            )

        time_now += params.dt
        next_step = current_step + 1

        trace_after_step = real(tr(Prj))
        trace_change_step = trace_after_step - trace_before_step

        # ------------------------------------------------------------
        # Refresh cutoff basis after arriving at next_step.
        # Then diagnostics for next_step refer to the refreshed basis.
        # ------------------------------------------------------------
        trace_change_refresh = 0.0

        if next_step % state.refresh_every == 0
            trace_before_refresh = real(tr(Prj))

            Prj, wave_work, book, proj, work =
                refresh_cutoff_basis(
                    Prj,
                    wave_work,
                    book,
                    proj,
                    Ashift,
                    params
                )

            trace_after_refresh = real(tr(Prj))
            trace_change_refresh = trace_after_refresh - trace_before_refresh
        end

        eigs_Prj = eigvals(Hermitian(Prj))

        diagnostics = (
            step_index = next_step,
            time_now = time_now,
            trace = real(tr(Prj)),
            trace_error = real(tr(Prj)) - params.filling,
            trace_change_step = trace_change_step,
            trace_change_refresh = trace_change_refresh,

            trace_change_old_basis = step_transport_diagnostics.trace_change_old_basis,
            transport_trace_change = step_transport_diagnostics.transport_trace_change,

            min_eigenvalue = minimum(eigs_Prj),
            max_eigenvalue = maximum(eigs_Prj),
            hermiticity_error = norm(Prj - Prj'),
            idempotency_error = norm(Prj * Prj - Prj),
            frobenius_norm_squared = sum(abs2, Prj),
            Ashift_x = Ashift[1],
            Ashift_y = Ashift[2],
            dimension = length(wave_work),
            chemical_potential = chemical_potential,
            particle_number_check = particle_number_check
        )

        push!(diagnostics_record, diagnostics)

        println(
            "step = ", next_step,
            " time = ", time_now,
            " Tr(P)-N = ", diagnostics.trace_error,
            " min/max eig = ", diagnostics.min_eigenvalue, " / ", diagnostics.max_eigenvalue,
            " ||P-P†|| = ", diagnostics.hermiticity_error,
            " ||P²-P|| = ", diagnostics.idempotency_error,
            " Ashift = ", Ashift,
            " dim = ", length(wave_work)
        )

        current_step = next_step
    end

    return nothing
end


function Densitymap_customize(
    xrange::Vector{Float64},
    yrange::Vector{Float64},
    DM::Matrix{ComplexF64},
    wave::Vector{Vector{Int64}},
    T1::Vector{Float64},
    T2::Vector{Float64}
)
    waveS = SVector{2, Int64}.(wave)

    ρ = DM
    N = length(waveS)

    coeff = Dict{SVector{2, Int64}, ComplexF64}()

    @inbounds for j in 1:N, i in 1:N
        d = waveS[i] - waveS[j]
        coeff[d] = get(coeff, d, 0.0 + 0.0im) + ρ[i, j]
    end

    deltas = collect(keys(coeff))
    c = ComplexF64[coeff[d] for d in deltas]

    M = SMatrix{2, 2, Float64}(hcat(T1, T2))

    kx = Vector{Float64}(undef, length(deltas))
    ky = Vector{Float64}(undef, length(deltas))

    @inbounds for n in eachindex(deltas)
        k = M * SVector{2, Float64}(deltas[n])
        kx[n] = k[1]
        ky[n] = k[2]
    end

    zvec = zeros(Float64, length(xrange), length(yrange))

    Threads.@threads for jb in eachindex(yrange)
        y = yrange[jb]

        @inbounds for ja in eachindex(xrange)
            x = xrange[ja]
            s = 0.0 + 0.0im

            @simd for n in eachindex(c)
                phase = kx[n] * x + ky[n] * y
                s += c[n] * cis(phase)
            end

            zvec[ja, jb] = real(s)
        end
    end

    return zvec
end
















function charge_density_file_path_for_step(args, step_index::Int, layer_index::Int)
    checkpoint_path = tdhf_file_path_for_step(args, step_index)

    charge_density_dir = joinpath(dirname(checkpoint_path), "CD_data")
    mkpath(charge_density_dir)

    return joinpath(
        charge_density_dir,
        "CD_layer$(layer_index)" * basename(checkpoint_path)
    )
end


function save_charge_density_one_layer!(
    args,
    Prj::Matrix{ComplexF64},
    wave_work::Vector{Vector{Int64}},
    proj,
    params::TDHFParams,
    step_index::Int,
    layer_index::Int,
    xrange::Vector{Float64},
    yrange::Vector{Float64}
)
    N = length(wave_work)

    densitymatrix_projected = zeros(ComplexF64, N, N)

    orbital_start = 2 * layer_index - 1
    orbital_stop = 2 * layer_index

    Threads.@threads for g1 in 1:N
        for li in orbital_start:orbital_stop
            s1a = proj.spinor_set[g1][li]

            @inbounds for g2 in 1:N
                ov = s1a * conj(proj.spinor_set[g2][li])
                densitymatrix_projected[g1, g2] += Prj[g1, g2] * ov
            end
        end
    end

    zvec = Densitymap_customize(
        xrange,
        yrange,
        densitymatrix_projected,
        wave_work,
        params.T1,
        params.T2
    )

    save_file_path = charge_density_file_path_for_step(args, step_index, layer_index)

    JLD2.jldsave(
        save_file_path;
        zvec = zvec,
        xrange = xrange,
        yrange = yrange,
        layer_index = layer_index,
        step_index = step_index,
        args = args
    )

    println("Saved charge density layer $(layer_index): ", save_file_path)

    return nothing
end


function save_charge_density_all_layers!(
    args,
    Prj::Matrix{ComplexF64},
    wave_work::Vector{Vector{Int64}},
    proj,
    params::TDHFParams,
    step_index::Int
)
    xrange = collect(range(0.0, stop = 1.1*norm(params.a1m), length = 160))
    yrange = collect(range(0.0, stop = 1.5 * norm(params.a1m), length = 160))

    for layer_index in params.NL-1:params.NL
        save_charge_density_one_layer!(
            args,
            Prj,
            wave_work,
            proj,
            params,
            step_index,
            layer_index,
            xrange,
            yrange
        )
    end

    return nothing
end