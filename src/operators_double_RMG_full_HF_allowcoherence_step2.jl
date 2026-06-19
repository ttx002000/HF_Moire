
using LinearAlgebra

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
#=
function find_FL(quasi_particle_energy::Vector{Float64},target_density::Float64,val_s::Float64,val_e::Float64,temp::Float64,Area::Float64,bg_particle_density::Float64)
  try_FL=(val_s+val_e)/2

  #stan=10^(-5)*target_density
  if target_density==0.0
    stan=10^(-9)
  else
    stan=abs(10^(-8)*target_density)
  end
  fermifactor=[1/(exp((quasi_particle_energy[ja]-try_FL)/temp)+1) for ja in eachindex(quasi_particle_energy)]
 
  fl=sum(fermifactor)/Area-bg_particle_density



 if abs(fl-target_density)<stan
  
    return try_FL,fl
  elseif fl-target_density>=stan
 
    return find_FL(quasi_particle_energy,target_density,val_s, try_FL,temp,Area,bg_particle_density)
  elseif fl-target_density<=-stan
 
    return find_FL(quasi_particle_energy,target_density,try_FL,val_e,temp,Area,bg_particle_density)
   end
 

end
=#



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


function Coulomb_matrix(z_pos::Vector{Vector{Float64}},kvec::Vector{Float64},NL::Int)::Matrix{Float64}
 
    spin_num=2
    layer_num=2
    valley_num=2
    sublattice_num=2*NL
    dimension=spin_num*layer_num*valley_num*sublattice_num

  s1=[Coulomb(abs(z_pos[1][ja]-z_pos[1][jb]),kvec) for ja in 1:sublattice_num, jb in 1:sublattice_num]
  s2=[Coulomb(abs(z_pos[1][ja]-z_pos[2][jb]),kvec) for ja in 1:sublattice_num, jb in 1:sublattice_num]

  cmatrix=zeros(Float64,valley_num,spin_num,layer_num,sublattice_num,valley_num,spin_num,layer_num,sublattice_num)
  for sindex1 in 1:spin_num, sindex2 in 1:spin_num,vindex1 in 1:valley_num, vindex2 in 1:valley_num
    cmatrix[vindex1,sindex1,1,:,vindex2,sindex2,1,:]+=s1
    cmatrix[vindex1,sindex1,1,:,vindex2,sindex2,2,:]+=s2
    cmatrix[vindex1,sindex1,2,:,vindex2,sindex2,1,:]+=s2'
    cmatrix[vindex1,sindex1,2,:,vindex2,sindex2,2,:]+=s1
  end
  
  return reshape(cmatrix,(dimension,dimension))

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

function get_single_particle(radius::Float64,num_points::Int,uD::Float64,CNP::Float64,NL::Int)

  vset=[1,-1]
  offset=[0.0,-(NL-1)*uD+CNP]

  kx_grid=collect(range(-radius/2, stop=+radius/2, length=num_points))
  ky_grid=collect(range(-radius/2, stop=+radius/2, length=num_points))
  Area=4*π^2/((kx_grid[2]-kx_grid[1])*(ky_grid[2]-ky_grid[1]))

  valley_num=2
  spin_num=2
  layer_num=2
  sublattice_num=2*NL
  dimension=valley_num*spin_num*layer_num*sublattice_num
  
  eig_set=[Vector{Float64}[] for _ in 1:spin_num, _ in 1:layer_num,_ in 1:valley_num]
  k_set=Vector{Float64}[]
  k_index=Vector{Int}[]
  eig_vec_set=[Matrix{ComplexF64}[] for _ in 1:spin_num, _ in 1:layer_num, _ in 1:valley_num]
  Ham_set=[Matrix{ComplexF64}[] for _ in 1:spin_num, _ in 1:layer_num, _ in 1:valley_num]
      
  for ja in eachindex(kx_grid),jb in eachindex(ky_grid)
    push!(k_set,[kx_grid[ja],ky_grid[jb]])
    push!(k_index,[ja,jb])

    for vi in 1:valley_num, si in 1:spin_num, li in 1:layer_num
      Ham=Hamiltonian([kx_grid[ja],ky_grid[jb]],uD,vset[vi],1,NL)
      FFF=eigen(Ham)
      push!(eig_set[si,li,vi],real(FFF.values).+offset[li])
      push!(Ham_set[si,li,vi],Ham+offset[li]*Matrix{Float64}(I,sublattice_num,sublattice_num))
      push!(eig_vec_set[si,li,vi],FFF.vectors)

    end
  end
  
   single_matrix_complex=zeros(ComplexF64,valley_num,spin_num,layer_num,sublattice_num,valley_num,spin_num,layer_num,sublattice_num,length(k_set))
  

  for ja in eachindex(k_set), vi in 1:valley_num, si in 1:spin_num, li in 1:layer_num
    single_matrix_complex[vi,si,li,:,vi,si,li,:,ja]+=Ham_set[si,li,vi][ja]
  end

  single_matrix=reshape(single_matrix_complex,(dimension,dimension,length(k_set)))

  return eig_set,k_set,k_index,eig_vec_set,Area,single_matrix

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
    temp::Float64,
    mixing::Float64,
    NL::Int
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

        FFF = eigen(Hermitian(Htmp))
        @views copy!(work.HF_eigenvectors[:, :, ja], FFF.vectors)
        @views copy!(work.HF_eigenvalues[:, ja], real(FFF.values))
    end

    BLAS.set_num_threads(old_blas_threads)

    work.evals_flat .= vec(work.HF_eigenvalues)
    val_s = minimum(work.evals_flat)
    val_e = maximum(work.evals_flat)
    bg_particle_density = dimension / (2 * Area) * Nk

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
    kill_intersheet_coherence!(work.output_density_matrix, NL)
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




function kill_intersheet_coherence!(
    A::Array{ComplexF64,3},
    NL::Int;
    valley_num::Int = 2,
    spin_num::Int = 2,
    sheet_num::Int = 2,
)
    sublattice_num = 2 * NL
    Nk = size(A, 3)

    A8 = reshape(
        A,
        valley_num, spin_num, sheet_num, sublattice_num,
        valley_num, spin_num, sheet_num, sublattice_num,
        Nk,
    )

    # sheet 1 <-> sheet 2 coherence
    @views A8[:, :, 1, :, :, :, 2, :, :] .= 0.0 + 0.0im
    @views A8[:, :, 2, :, :, :, 1, :, :] .= 0.0 + 0.0im

    return A
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

@inline function symmetrize_from_lower_2d!(A::Matrix{ComplexF64})
    n = size(A,1)
    @assert n == size(A,2)
    @inbounds for i in 1:n
        A[i,i] = complex(real(A[i,i]), 0.0)
        for j in (i+1):n
            A[i,j] = conj(A[j,i])
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

function pick_random_jld2_path(folder::AbstractString)
    isdir(folder) || return nothing

    jld2_files = filter(readdir(folder; join=true)) do p
        isfile(p) && endswith(lowercase(p), ".jld2")
    end

    isempty(jld2_files) && return nothing
    return rand(jld2_files)
end



function get_initial_proj(k_set::Vector{Vector{Float64}},eig_vec_set::Array{Vector{Matrix{ComplexF64}}},NL::Int)
 
 valley_num=2
 spin_num=2
 layer_num=2
 sublattice_num=2*NL
 dimension=valley_num*spin_num*layer_num*sublattice_num

   BG_density_matrix_complex=zeros(ComplexF64,valley_num,spin_num,layer_num,sublattice_num,valley_num,spin_num,layer_num,sublattice_num,length(k_set))

   for ja in eachindex(k_set), bandindex in 1:NL, vi in 1:valley_num, si in 1:spin_num, li in 1:layer_num
    BG_density_matrix_complex[vi,si,li,:,vi,si,li,:,ja]+=eig_vec_set[si,li,vi][ja][:,bandindex]*(eig_vec_set[si,li,vi][ja][:,bandindex])'
   end

   BG_density_matrix=reshape(BG_density_matrix_complex,(dimension,dimension,length(k_set)))

 initial_density_matrix=zeros(ComplexF64,dimension,dimension,length(k_set))

 for ja in eachindex(k_set)
   A=randn(dimension,dimension)+im*randn(dimension,dimension)
   #A=0.5*ones(dimension,dimension)
   initial_density_matrix[:,:,ja]+=(A+A')*0.1
 end



 return initial_density_matrix, BG_density_matrix
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




function iteration(initial_density_matrix::Array{ComplexF64},BG_density_matrix::Array{ComplexF64},
                           ϵr::Float64,k_set::Vector{Vector{Float64}},single_matrix::Array{ComplexF64},ildis::Float64,Area::Float64,NL::Int,target_density::Float64,temp::Float64)
   valley_num = 2
    spin_num = 2
    layer_num = 2
    sublattice_num = 2 * NL
    dimension = valley_num * spin_num * layer_num * sublattice_num
    Nk = length(k_set)

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
  



   
  z_pos=zeros(Float64,layer_num,sublattice_num)
  z_pos[1,:]=0.335*[i for i in 0:NL-1 for _ in 1:2]
  z_pos[2,:]=0.335*[i for i in -NL+1:0 for _ in 1:2].-ildis
  z_pos=reshape(z_pos,layer_num*sublattice_num)
  sublayer_i=zeros(Int,valley_num,spin_num,sublattice_num*layer_num)
  for ja in 1:valley_num,jb in 1:spin_num, jc in 1:sublattice_num*layer_num
    sublayer_i[ja,jb,jc]=jc
  end
  sublayer_i=reshape(sublayer_i,dimension)

  fcmatrix=Matrix{Matrix{Float64}}(undef,dimension,dimension)
  smaller_fc=Matrix{Matrix{Float64}}(undef,layer_num*sublattice_num,layer_num*sublattice_num)

  tic=time()

  Threads.@threads :greedy for ja in 1:sublattice_num*layer_num
   for jb in 1:ja
      smaller_fc[ja,jb]=zeros(Float64,length(k_set),length(k_set))
        for jc in 1:length(k_set), jd in 1:jc
         smaller_fc[ja,jb][jc,jd]+=Coulomb(abs(z_pos[ja]-z_pos[jb]),k_set[jc]-k_set[jd])
        end

        for jc in 1:length(k_set), jd in jc+1:length(k_set)
         smaller_fc[ja,jb][jc,jd]+=smaller_fc[ja,jb][jd,jc]
        end
   end
 end

 Threads.@threads :greedy for ja in 1:dimension
   for jb in 1:ja

      sub_one=max(sublayer_i[ja],sublayer_i[jb])
      sub_two=min(sublayer_i[ja],sublayer_i[jb])
     fcmatrix[ja,jb]=smaller_fc[sub_one,sub_two]
       
   end
 end

smaller_fc=nothing


      hfcmatrix=zeros(Float64,dimension,dimension)
  for ja in 1:dimension
    for jb in 1:ja
        hfcmatrix[ja,jb]+=fcmatrix[ja,jb][1,1]
    end
    for jb in ja+1:dimension
      hfcmatrix[ja,jb]+=fcmatrix[jb,ja][1,1]
    end
  end


 toc=time()
 println("formfactorstime",toc-tic)






while (eout > 1e-14) || (bad_count < DIIS_size + 2) || (abs(energy_change) > 1e-8)
      
 
        if eout < 1e-14
            bad_count += 1
        else
            bad_count = 0
        end

        dmk_used = input_density_matrix
        tic = time()

        use_diis = ((itcount > 60 && abs(eout) > 1e-2) ||
                    (itcount > 50 && abs(eout) < 1e-7) ||
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
                temp,
                mixing,
                NL
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
                temp,
                mixing,
                NL
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
      if inB≠0
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
        return 0
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
          rethrow(e)  # If another error occurs, propagate it
      end
  end
end