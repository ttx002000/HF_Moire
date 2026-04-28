
using LinearAlgebra
#To do
#(1) Check the units of m_{TMD} in constructing the single particle Hamiltonian

struct BasisSelection
    # RMG blocks are explicit (spin, valley) pairs.
    # Empty means no RMG.
    rmg_blocks::Vector{Tuple{Int,Int}}

    # TMD spins to include.
    # Empty means no TMD.
    tmd_spins::Vector{Int}
end

function BasisSelection(;
    rmg_blocks::Vector{Tuple{Int,Int}} = [(spin, valley) for valley in 1:2 for spin in 1:2],
    tmd_spins::Vector{Int} = [1, 2],
)
    isempty(rmg_blocks) && isempty(tmd_spins) &&
        error("BasisSelection cannot have both no RMG and no TMD.")
    return BasisSelection(rmg_blocks, tmd_spins)
end




struct BasisState
    material::Symbol
    spin::Int
    valley::Int
    layer::Int
    sublat::Int
    z::Float64
    channel::Int
end



function build_basis_RMG_TMD(
    NL::Int,
    z0_RMG::Float64,
    z_TMD::Float64;
    selection::BasisSelection = BasisSelection(),
)
    basis = BasisState[]

    # RMG channels: 1, ..., 2NL
    for (spin, valley) in selection.rmg_blocks
        for layer in 1:NL
            for sublat in 1:2
                z = z0_RMG + 0.335 * (layer - 1)

                # This matches the order of your RMG Hamiltonian
                channel = 2 * (layer - 1) + sublat

                push!(basis, BasisState(
                    :RMG,
                    spin,
                    valley,
                    layer,
                    sublat,
                    z,
                    channel,
                ))
            end
        end
    end

    # TMD channel: 2NL + 1
    tmd_channel = 2 * NL + 1

    for spin in selection.tmd_spins
        push!(basis, BasisState(
            :TMD,
            spin,
            0,
            1,
            1,
            z_TMD,
            tmd_channel,
        ))
    end

    isempty(basis) && error("The basis is empty. Include at least one RMG or TMD state.")

    return basis
end



function rmg_block_indices(basis::Vector{BasisState}, spin::Int, valley::Int, NL::Int)
    inds = Int[]

    for layer in 1:NL
        for sublat in 1:2 # This must be in this order
            idx = findfirst(b -> b.material == :RMG &&
                                b.spin == spin &&
                                b.valley == valley &&
                                b.layer == layer &&
                                b.sublat == sublat,
                            basis)

            idx === nothing && error(
                "Missing RMG basis state: spin=$spin valley=$valley layer=$layer sublat=$sublat"
            )

            push!(inds, idx)
        end
    end

    return inds
end


function tmd_index(basis::Vector{BasisState}, spin::Int)
    idx = findfirst(b -> b.material == :TMD && b.spin == spin, basis)
    idx === nothing && error("Missing TMD basis state: spin=$spin")
    return idx
end




function rmg_blocks_present(basis::Vector{BasisState})
    blocks = Tuple{Int,Int}[]

    for b in basis
        if b.material == :RMG
            block = (b.spin, b.valley)
            if !(block in blocks)
                push!(blocks, block)
            end
        end
    end

    return blocks
end


function tmd_spins_present(basis::Vector{BasisState})
    spins = Int[]

    for b in basis
        if b.material == :TMD
            if !(b.spin in spins)
                push!(spins, b.spin)
            end
        end
    end

    return spins
end


function has_RMG(basis::Vector{BasisState})
    return any(b -> b.material == :RMG, basis)
end


function has_TMD(basis::Vector{BasisState})
    return any(b -> b.material == :TMD, basis)
end







mutable struct HFWorkRMG
    Fock_matrix::Array{ComplexF64,3}          # (dim, dim, Nk)
    Hartree_matrix::Matrix{ComplexF64}        # (dim, dim)
    density_matrix_new::Array{ComplexF64,3}   # (dim, dim, Nk)
    output_density_matrix::Array{ComplexF64,3}
    DeltaMatrix::Array{ComplexF64,3}

    HF_eigenvalues::Matrix{Float64}           # (dim, Nk)
    HF_eigenvectors::Array{ComplexF64,3}      # (dim, dim, Nk)

    # workspace for finite-T occupations / temporary vectors
    fermifactor::Vector{Float64}
    evals_flat::Vector{Float64}
    occ_matrix::Matrix{Float64}               # (dim, Nk)

    # thread-local scratch for diagonalization / DM reconstruction
    H_scratch::Vector{Matrix{ComplexF64}}     # one dense (dim,dim) matrix per Julia thread
    Vocc_scratch::Vector{Matrix{ComplexF64}}  # one dense (dim,dim) matrix per Julia thread

    # DIIS storage
    DIIS_input_density_matrix::Vector{Array{ComplexF64,3}}
    DIIS_input_DeltaMatrix::Vector{Array{ComplexF64,3}}
    diis_head::Int
    diis_len::Int

    # optional scratch for Hartree diagonal accumulation
    density_sum_k::Matrix{ComplexF64}         # (dim, dim)
end


function Coulomb(dis::Float64,kvec::Vector{Float64})::Float64
  
  if norm(kvec)==0.0
      return 9047.5636*(-dis)
    # return 9047.5636*0.0
    else
     return 9047.5636/norm(kvec)*exp(-norm(kvec)*dis)
     #return 9047.5636/norm(kvec)
  end
  



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
    target_density::Float64,
    val_s::Float64,
    val_e::Float64,
    temp::Float64,
    Area::Float64,
    bg_particle_density::Float64;
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

        fl = sum(fermifactor) / Area - bg_particle_density

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

function get_single_particle(radius::Float64,num_points::Int,
    uD::Float64,deltaE::Float64,NL::Int,m_TMD::Float64,z_TMD::Float64;
    selection::BasisSelection = BasisSelection(),)

  vset=[1,-1]
  stacking=1
  if m_TMD>0
       CNP_TMD=-abs((NL-1)*uD/2)+deltaE # Check this
  else
       CNP_TMD=abs((NL-1)*uD/2)-deltaE # Check this
  end
  kx_grid=collect(range(-radius/2, stop=+radius/2, length=num_points))
  ky_grid=collect(range(-radius/2, stop=+radius/2, length=num_points))
  Area=4*π^2/((kx_grid[2]-kx_grid[1])*(ky_grid[2]-ky_grid[1]))

  
  k_set=Vector{Float64}[]
  k_index=Vector{Int}[]
   for ix in eachindex(kx_grid), iy in eachindex(ky_grid)
        push!(k_set, [kx_grid[ix], ky_grid[iy]])
        push!(k_index, [ix, iy])
   end

    basis = build_basis_RMG_TMD(
        NL,
        0.0,
        z_TMD;
        selection = selection,
    )

    dimension = length(basis)
    Nk = length(k_set)

    single_matrix = zeros(ComplexF64, dimension, dimension, Nk)
    
   rmg_blocks = rmg_blocks_present(basis)
    tmd_spins = tmd_spins_present(basis)

    rmg_inds = Dict{Tuple{Int,Int}, Vector{Int}}()
    for block in rmg_blocks
        spin, valley = block
        rmg_inds[block] = rmg_block_indices(basis, spin, valley, NL)
    end

    tmd_inds = Dict{Int,Int}()
    for spin in tmd_spins
        tmd_inds[spin] = tmd_index(basis, spin)
    end

  



  for ik in eachindex(k_set)
        k = k_set[ik]

        # RMG blocks
        for block in rmg_blocks
            spin, valley = block
            inds = rmg_inds[block]

            H_RMG = Hamiltonian(k, uD, vset[valley], stacking, NL)

            @views single_matrix[inds, inds, ik] .= H_RMG
        end

        # TMD blocks
        for spin in tmd_spins
            i = tmd_inds[spin]

            ε_TMD = dot(k, k) / (2 * m_TMD) * 76.1996 + CNP_TMD

            single_matrix[i, i, ik] = ε_TMD
        end
    end




    # ------------------------------------------------------------
    # Save noninteracting single-particle eigenvalues/eigenvectors
    # blockwise
    # ------------------------------------------------------------

     rmg_eig_set = Dict{Tuple{Int,Int}, Matrix{Float64}}()
    rmg_eig_vec_set = Dict{Tuple{Int,Int}, Array{ComplexF64,3}}()

    for block in rmg_blocks
        inds = rmg_inds[block]
        dim_rmg = length(inds)

        vals = zeros(Float64, dim_rmg, Nk)
        vecs = zeros(ComplexF64, dim_rmg, dim_rmg, Nk)

        for ik in 1:Nk
            H_RMG = single_matrix[inds, inds, ik]
            F = eigen(Hermitian(H_RMG))

            vals[:, ik] .= real(F.values)
            vecs[:, :, ik] .= F.vectors
        end

        rmg_eig_set[block] = vals
        rmg_eig_vec_set[block] = vecs
    end

    tmd_eig_set = Dict{Int, Vector{Float64}}()
    tmd_eig_vec_set = Dict{Int, Vector{ComplexF64}}()

    for spin in tmd_spins
        i = tmd_inds[spin]

        vals = zeros(Float64, Nk)
        for ik in 1:Nk
            vals[ik] = real(single_matrix[i, i, ik])
        end

        tmd_eig_set[spin] = vals
        tmd_eig_vec_set[spin] = ComplexF64[1.0 + 0.0im]
    end

    single_particle_wf = (
        rmg_eig_set = rmg_eig_set,
        rmg_eig_vec_set = rmg_eig_vec_set,
        tmd_eig_set = tmd_eig_set,
        tmd_eig_vec_set = tmd_eig_vec_set,
    )









   return basis, k_set, k_index, Area, single_matrix, single_particle_wf

end






function Construct_projector!(
    work::HFWorkRMG,
    k_set::Vector{Vector{Float64}},
    ϵr::Float64,
    density_matrix::Array{ComplexF64,3},
    single_matrix::Array{ComplexF64,3},
    energy_input::Float64,
    BG_density_matrix::Array{ComplexF64,3},
    hfcmatrix::Matrix{Float64},
    fcmatrix::Matrix{Matrix{Float64}},
    Area::Float64,
    dimension::Int,
    target_density::Float64,
    bg_particle_density::Float64,
    temp::Float64,
    mixing::Float64,
)
  


  Nk = length(k_set)

    fill!(work.Fock_matrix, 0)
    fill!(work.Hartree_matrix, 0)

tic=time()

 Threads.@threads :greedy for ja in 1:dimension
  for jb in 1:(ja - 1)
      mul!(
          @view(work.Fock_matrix[ja, jb, :]),
          @view(fcmatrix[ja, jb][:, :]),
          @view(density_matrix[ja, jb, :]),
          1,
          1,
      )
  end
 end





 Threads.@threads for ja in 1:dimension
  mul!(
      @view(work.Fock_matrix[ja, ja, :]),
      @view(fcmatrix[ja, ja][:, :]),
      @view(density_matrix[ja, ja, :]),
      1,
      1,
  )
 end
  symmetrize_from_lower_3d!(work.Fock_matrix)
    work.Fock_matrix .*= 1 / (ϵr * Area)


 toc=time()
 println("Focktime",toc-tic)
 

 tic=time()
   fill!(work.density_sum_k, 0)
   work.density_sum_k .= dropdims(sum(density_matrix, dims=3), dims=3)

   
    diag_density = diag(work.density_sum_k)
    hartree_diag = hfcmatrix * diag_density
    @inbounds for ja in 1:dimension
        work.Hartree_matrix[ja, ja] = hartree_diag[ja] / (ϵr * Area)
    end

toc=time()
 println("Hartreetime",toc-tic)



            tic = time()

    # keep Julia threading for the k-loop, but avoid BLAS oversubscription
    old_blas_threads = BLAS.get_num_threads()
    BLAS.set_num_threads(1)

    Threads.@threads for ja in eachindex(k_set)
        tid = Threads.threadid()
        Htmp = work.H_scratch[tid]

        copy!(Htmp, work.Hartree_matrix)
        @views Htmp .-= work.Fock_matrix[:, :, ja]
        @views Htmp .+= single_matrix[:, :, ja]

        FFF = eigen!(Hermitian(Htmp))
        @views copy!(work.HF_eigenvectors[:, :, ja], FFF.vectors)
        @views copy!(work.HF_eigenvalues[:, ja], real(FFF.values))
    end

    BLAS.set_num_threads(old_blas_threads)

    work.evals_flat .= vec(work.HF_eigenvalues)
    val_s = minimum(work.evals_flat)
    val_e = maximum(work.evals_flat)

    fermi_level, renormalized_density = find_FL_iterative!(
        work.fermifactor,
        work.evals_flat,
        target_density,
        val_s,
        val_e,
        temp,
        Area,
        bg_particle_density,
    )

    @inbounds for ja in 1:Nk, jb in 1:dimension
        work.occ_matrix[jb, ja] = fermi_occ((work.HF_eigenvalues[jb, ja] - fermi_level) / temp)
    end

    Threads.@threads for ja in 1:Nk
        tid = Threads.threadid()
        Vocc = work.Vocc_scratch[tid]

        @views copy!(Vocc, work.HF_eigenvectors[:, :, ja])

        @inbounds for band in 1:dimension
            occ = work.occ_matrix[band, ja]
            @views Vocc[:, band] .*= occ
        end

        mul!(
            @view(work.density_matrix_new[:, :, ja]),
            Vocc,
            adjoint(@view(work.HF_eigenvectors[:, :, ja])),
            1.0,
            0.0,
        )

        @views work.density_matrix_new[:, :, ja] .-= BG_density_matrix[:, :, ja]
    end

    toc = time()
    println("constructing DMtime", toc - tic)

   tic=time()
    @. work.output_density_matrix = mixing * density_matrix + (1 - mixing) * work.density_matrix_new
    @. work.DeltaMatrix = work.density_matrix_new - density_matrix



    eout = real(sum(abs2, work.DeltaMatrix) / Nk)
    energy = 0.0
    @inbounds for ja in eachindex(k_set)
        energy += real(tr(
            density_matrix[:,:,ja] *
            (work.Hartree_matrix/2 - work.Fock_matrix[:,:,ja]/2 + single_matrix[:,:,ja])
        )) / Nk
    end

    energy_change = energy - energy_input

 energy_change=real(energy-energy_input)

 toc=time()
  println("othertime",toc-tic)

    return eout, energy_change, real(energy), fermi_level, renormalized_density

end



@inline function symmetrize_from_lower_3d!(A::Array{ComplexF64,3})
    n1, n2, nk = size(A)
    @assert n1 == n2
    @inbounds for k in 1:nk
        for i in 1:n1
            A[i,i,k] = complex(real(A[i,i,k]), 0.0)
            for j in (i+1):n1
                A[i,j,k] = conj(A[j,i,k])
            end
        end
    end
    return A
end




@inline function diis_push!(
    DIIS_input_density_matrix::Vector{Array{ComplexF64,3}},
    DIIS_input_DeltaMatrix::Vector{Array{ComplexF64,3}},
    input_density_matrix::Array{ComplexF64,3},
    DeltaMatrix::Array{ComplexF64,3},
    diis_head::Int,
    diis_len::Int
)
    copy!(DIIS_input_density_matrix[diis_head], input_density_matrix)
    copy!(DIIS_input_DeltaMatrix[diis_head],  DeltaMatrix)

    diis_head = (diis_head == length(DIIS_input_density_matrix)) ? 1 : (diis_head + 1)
    diis_len  = min(diis_len + 1, length(DIIS_input_density_matrix))
    return diis_head, diis_len
end




function get_initial_proj(
    k_set::Vector{Vector{Float64}},
    basis::Vector{BasisState},
    single_matrix::Array{ComplexF64,3},
    NL::Int,
    m_TMD::Float64,
)
    dimension = length(basis)
    Nk = length(k_set)

    BG_density_matrix = zeros(ComplexF64, dimension, dimension, Nk)

    rmg_blocks = rmg_blocks_present(basis)
    tmd_spins = tmd_spins_present(basis)

    rmg_inds = Dict{Tuple{Int,Int}, Vector{Int}}()
    for block in rmg_blocks
        spin, valley = block
        rmg_inds[block] = rmg_block_indices(basis, spin, valley, NL)
    end

    tmd_inds = Dict{Int,Int}()
    for spin in tmd_spins
        tmd_inds[spin] = tmd_index(basis, spin)
    end

 # RMG background:
    # For every included RMG spin/valley block, fill the NL valence bands.
    for ik in 1:Nk
        for block in rmg_blocks
            inds = rmg_inds[block]

            H_RMG = single_matrix[inds, inds, ik]
            F = eigen(Hermitian(H_RMG))

            for band in 1:NL
                v = F.vectors[:, band]
                BG_density_matrix[inds, inds, ik] .+= v * v'
            end
        end
    end

    if m_TMD < 0
        for ik in 1:Nk
            for spin in tmd_spins
                i = tmd_inds[spin]
                BG_density_matrix[i, i, ik] = 1.0
            end
        end
    end


   

    # Random initial fluctuation relative to background
    initial_density_matrix = zeros(ComplexF64, dimension, dimension, Nk)

    for ik in 1:Nk
        A = randn(dimension, dimension) .+ im .* randn(dimension, dimension)
        initial_density_matrix[:, :, ik] .=1.0 .* (A + A')
    end

    return initial_density_matrix, BG_density_matrix
end



function background_density(BG_density_matrix::Array{ComplexF64,3}, Area::Float64)
    Nk = size(BG_density_matrix, 3)
    total = 0.0

    for ik in 1:Nk
        total += real(tr(BG_density_matrix[:, :, ik]))
    end

    return total / Area
end


function HFWorkRMG(dimension::Int, Nk::Int, DIIS_size::Int)
    Z3 = zeros(ComplexF64, dimension, dimension, Nk)
    Z2 = zeros(ComplexF64, dimension, dimension)

    nscratch = Threads.maxthreadid()

    return HFWorkRMG(
        copy(Z3),  # Fock_matrix
        copy(Z2),  # Hartree_matrix
        copy(Z3),  # density_matrix_new
        copy(Z3),  # output_density_matrix
        copy(Z3),  # DeltaMatrix

        zeros(Float64, dimension, Nk),  # HF_eigenvalues
        copy(Z3),                       # HF_eigenvectors

        zeros(Float64, dimension * Nk), # fermifactor
        zeros(Float64, dimension * Nk), # evals_flat
        zeros(Float64, dimension, Nk),  # occ_matrix

        [zeros(ComplexF64, dimension, dimension) for _ in 1:nscratch], # H_scratch
        [zeros(ComplexF64, dimension, dimension) for _ in 1:nscratch], # Vocc_scratch

        [zeros(ComplexF64, dimension, dimension, Nk) for _ in 1:DIIS_size], # DIIS_input_density_matrix
        [zeros(ComplexF64, dimension, dimension, Nk) for _ in 1:DIIS_size], # DIIS_input_DeltaMatrix
        1,  # diis_head
        0,  # diis_len

        copy(Z2),  # density_sum_k
    )
end










function build_smaller_fc_RMG_TMD(
    k_set::Vector{Vector{Float64}},
    NL::Int,
    z0_RMG::Float64,
    z_TMD::Float64,
)
    Nk = length(k_set)

    nchannels = 2 * NL + 1

    channel_z = zeros(Float64, nchannels)

    # RMG channels
    for layer in 1:NL, sublat in 1:2
        channel = 2 * (layer - 1) + sublat
        channel_z[channel] = z0_RMG + 0.335 * (layer - 1)
    end

    # TMD channel
    channel_z[2 * NL + 1] = z_TMD

    smaller_fc = Matrix{Matrix{Float64}}(undef, nchannels, nchannels)

    Threads.@threads :greedy for a in 1:nchannels
        for b in 1:a
            fab = zeros(Float64, Nk, Nk)

            dz = abs(channel_z[a] - channel_z[b])

            for k1 in 1:Nk, k2 in 1:k1
                fab[k1, k2] = Coulomb(dz, k_set[k1] - k_set[k2])
            end

            for k1 in 1:Nk, k2 in k1+1:Nk
                fab[k1, k2] = fab[k2, k1]
            end

            smaller_fc[a, b] = fab
        end
    end

    return smaller_fc, channel_z
end







function build_fcmatrix_from_smaller(
    basis::Vector{BasisState},
    smaller_fc::Matrix{Matrix{Float64}},
)
    dimension = length(basis)

    fcmatrix = Matrix{Matrix{Float64}}(undef, dimension, dimension)

    Threads.@threads :greedy for i in 1:dimension
        for j in 1:i
            ci = basis[i].channel
            cj = basis[j].channel

            cmax = max(ci, cj)
            cmin = min(ci, cj)

            # This is a reference, not a copy
            fcmatrix[i, j] = smaller_fc[cmax, cmin]
        end
    end

    return fcmatrix
end


function build_hfcmatrix_from_fcmatrix(
    fcmatrix::Matrix{Matrix{Float64}},
)
    dimension = size(fcmatrix, 1)

    hfcmatrix = zeros(Float64, dimension, dimension)

    for i in 1:dimension
        for j in 1:i
            hfcmatrix[i, j] = fcmatrix[i, j][1, 1]
        end

        for j in i+1:dimension
            hfcmatrix[i, j] = fcmatrix[j, i][1, 1]
        end
    end

    return hfcmatrix
end




function iteration(
            initial_density_matrix::Array{ComplexF64,3},
            BG_density_matrix::Array{ComplexF64,3},
            ϵr::Float64,
            k_set::Vector{Vector{Float64}},
            single_matrix::Array{ComplexF64,3},
            basis::Vector{BasisState},
            Area::Float64,
            NL::Int,
            z_TMD::Float64,
            target_density::Float64,
            temp::Float64
)
    dimension = length(basis)
    Nk = length(k_set)

    bg_particle_density=background_density(BG_density_matrix, Area)


    @assert size(single_matrix, 1) == dimension
    @assert size(single_matrix, 2) == dimension
    @assert size(single_matrix, 3) == Nk

    @assert size(initial_density_matrix) == (dimension, dimension, Nk)
    @assert size(BG_density_matrix) == (dimension, dimension, Nk)

    DIIS_size = 5
    mixing = 0.5

    eout = 1.0
    itcount = 0
    bad_count = 0
    energy = 0.0
    energy_change = 0.0
    fermi_level = 0.0
    renormalized_density = 0.0

    input_density_matrix = copy(initial_density_matrix)

    work = HFWorkRMG(dimension, Nk, DIIS_size)

    eout_hist = Float64[]
    PLATEAU_N = 10
    PLATEAU_FRAC = 0.10
    E_EPS = 1e-30

    diis_fire_once = false
    diis_cooldown = 0
    DIIS_COOLDOWN = 10
  



   
  tic = time()

    smaller_fc, channel_z = build_smaller_fc_RMG_TMD(
        k_set,
        NL,
        0.0,
        z_TMD,
    )

    fcmatrix = build_fcmatrix_from_smaller(
        basis,
        smaller_fc,
    )
    smaller_fc=nothing
    hfcmatrix = build_hfcmatrix_from_fcmatrix(fcmatrix)

    toc = time()
    println("formfactorstime", toc - tic)






while (eout > 1e-14) || (bad_count < DIIS_size + 2) || (abs(energy_change) > 1e-8)
   
 
        if eout < 1e-14
            bad_count += 1
        else
            bad_count = 0
        end

        dmk_used = input_density_matrix
        tic = time()

        use_diis = ((itcount > 60 && abs(eout) > 1e-2) ||
                    (itcount > 50 && abs(eout) < 1e-6) ||
                    diis_fire_once) && (work.diis_len >= 4)

        if use_diis
            dmk = implement_DIIS(
                work.DIIS_input_density_matrix,
                work.DIIS_input_DeltaMatrix,
                k_set,
                work.diis_len
            )

            if dmk === nothing
                A = randn(dimension, dimension, Nk) .+ im * randn(dimension, dimension, Nk)
                dmk = similar(input_density_matrix)
                @inbounds for k in 1:Nk
                    dmk[:,:,k] .= 0.01 .* (A[:,:,k] + A[:,:,k]')
                end
                itcount = 0
                work.diis_head = 1
                work.diis_len = 0
                empty!(eout_hist)
                diis_fire_once = false
                diis_cooldown = 0
                println("DIIS failed, random restart")
            end

            dmk_used = dmk
            eout, energy_change, energy, fermi_level, renormalized_density = Construct_projector!(
                work,
                k_set, ϵr,
                dmk,
                single_matrix,
                energy,
                BG_density_matrix,
                hfcmatrix,
                fcmatrix,
                Area,
                dimension,
                target_density,
                bg_particle_density,
                temp,
                mixing,
            )
            println("using DIIS")

            if diis_fire_once
                diis_fire_once = false
                empty!(eout_hist)
                diis_cooldown = DIIS_COOLDOWN
            end
        else
            eout, energy_change, energy, fermi_level, renormalized_density = Construct_projector!(
                work,
                k_set, ϵr,
                input_density_matrix,
                single_matrix,
                energy,
                BG_density_matrix,
                hfcmatrix,
                fcmatrix,
                Area,
                dimension,
                target_density,
                bg_particle_density,
                temp,
                mixing,
            )
        end

        work.diis_head, work.diis_len = diis_push!(
            work.DIIS_input_density_matrix,
            work.DIIS_input_DeltaMatrix,
            dmk_used,
            work.DeltaMatrix,
            work.diis_head,
            work.diis_len
        )

        copy!(input_density_matrix, work.output_density_matrix)

        itcount += 1
        println(time() - tic, "eout=$eout", "energy_change=$energy_change", "itcount=$itcount")
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
           work.DIIS_input_density_matrix,
           fermi_level,
           work.Hartree_matrix,
           work.Fock_matrix,
           eout,
           renormalized_density
end




function implement_DIIS(DIIS_input_projector::Vector{Array{ComplexF64,3}},DIIS_input_DeltaMatrix::Vector{Array{ComplexF64,3}},k_set::Vector{Vector{Float64}},DIIS_size::Int)



      Bmatrix=zeros(ComplexF64,DIIS_size+1,DIIS_size+1)
      for ja in 1:DIIS_size
       Bmatrix[ja,DIIS_size+1]=1
       Bmatrix[DIIS_size+1,ja]=1
      end
  
      for ja in 1:DIIS_size,jb in 1:DIIS_size
          for jc in eachindex(k_set)
             Bmatrix[ja,jb]+=real(tr((DIIS_input_DeltaMatrix[ja][:,:,jc])'*(DIIS_input_DeltaMatrix[jb][:,:,jc])))
          end
      end

      inB=safe_inverse(Bmatrix)
      if inB!==nothing
         onh=zeros(Float64,DIIS_size+1)
         onh[end]=1
         coeff=inB* onh

        dmk = zero(DIIS_input_projector[1])           # alloc once
        @inbounds for ja in 1:DIIS_size
            BLAS.axpy!(coeff[ja], DIIS_input_projector[ja], dmk)  # dmk += coeff[ja] * projector[ja]
       
            BLAS.axpy!(coeff[ja], DIIS_input_DeltaMatrix[ja], dmk)  # dmk += coeff[ja] * projector[ja]
        end
         #dmk=coeff[1]*(DIIS_input_projector[1]+DIIS_input_DeltaMatrix[1])+coeff[2]*(DIIS_input_projector[2]+DIIS_input_DeltaMatrix[2])+coeff[3]*(DIIS_input_projector[3]+DIIS_input_DeltaMatrix[3])
         return dmk
      else
        return nothing
      end
end



function safe_inverse(A)
  try
      return inv(A)  # Attempt to compute inverse
  catch e
      if isa(e, SingularException)
          println("Matrix is singular, doing pseudoinverse.")
          return pinv(A, 10^(-8))  # Use pseudoinverse as an alternative
      else
          return nothing  # If another error occurs, propagate it
      end
  end
end

function get_selection(sel::Int)
    if sel==1
     return BasisSelection(
        rmg_blocks = [(1, 1)],
        tmd_spins = [1],
        )
    elseif sel==2

        return BasisSelection(
            rmg_blocks = [(1, 1),(1, 2),(2, 1),(2, 2)],
            tmd_spins = [1],
        )
    elseif sel==3

        return BasisSelection(
            rmg_blocks = [(1, 1),(1, 2),(2, 1),(2, 2)],
            tmd_spins = [1,2],
        )   
    end

end