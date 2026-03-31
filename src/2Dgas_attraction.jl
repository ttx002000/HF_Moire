using LinearAlgebra
using Arpack
using Combinatorics
using Random
BLAS.set_num_threads(1)





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







function triangle_initial_Densitymatrix(scale::Float64,gcutoff::Float64)
  
  
    b1=scale*[1,0]
    b2=scale*[-1/2,√3/2]
    T1=b1
    T2=b2

    a1m=inv([b1';b2'])*[2π,0]
    a2m=inv([b1';b2'])*[0,2π]
    am=norm(a1m);
    
    Area=am^2
    b1T=Int.(round.(inv([T1 T2])*b1))
    b2T=Int.(round.(inv([T1 T2])*b2))




    
    
    
    wave=Vector{Int64}[]
    cutoffstandard=gcutoff*norm(b1)
    cutoff=Int(round(gcutoff))*5
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        gtest=ja*b1+jb*b2;
        if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
            push!(wave,ja*b1T+jb*b2T)
        end
    end   
    
    dimension=length(wave)

    single_Ham=zeros(ComplexF64,dimension,dimension)

    
      for jb in eachindex(wave)     
          single_Ham[jb,jb]=norm([T1 T2]*wave[jb])^4
      end


     

      A=randn(ComplexF64,dimension,dimension)
      input_DensityMatrix=(A+A')*1.0
   

    

    
    return wave, input_DensityMatrix, single_Ham, T1, T2, a1m, a2m, b1,b2,Area
      

       
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


@inline function diis_push!(
    DIIS_input_DensityMatrix::Vector{Matrix{ComplexF64}},
    DIIS_input_DeltaMatrix::Vector{Matrix{ComplexF64}},
    input_DensityMatrix::Matrix{ComplexF64},
    DeltaMatrix::Matrix{ComplexF64},
    diis_head::Int,
    diis_len::Int
)
    copy!(DIIS_input_DensityMatrix[diis_head], input_DensityMatrix)
    copy!(DIIS_input_DeltaMatrix[diis_head],  DeltaMatrix)

    diis_head = (diis_head == length(DIIS_input_DensityMatrix)) ? 1 : (diis_head + 1)
    diis_len  = min(diis_len + 1, length(DIIS_input_DensityMatrix))
    return diis_head, diis_len
end



mutable struct HFWork
    HartreeMatrix::Matrix{ComplexF64}
    FockMatrix::Matrix{ComplexF64}
    H_phys::Matrix{ComplexF64}
    NewDensityMatrix::Matrix{ComplexF64}
    DeltaMatrix::Matrix{ComplexF64}
    output_DensityMatrix::Matrix{ComplexF64}

    HF_eigenvalue::Vector{Float64}
    HF_eigenvector::Matrix{ComplexF64}


    DIIS_input_DensityMatrix::Vector{Matrix{ComplexF64}}
    DIIS_input_DeltaMatrix::Vector{Matrix{ComplexF64}}
    diis_head::Int
    diis_len::Int
    HartreeAccShift::Vector{ComplexF64}
end

function HFWork(dimension::Int, DIIS_size::Int, NSHIFT::Int)
    Z = zeros(ComplexF64, dimension, dimension)
    HFWork(
        copy(Z), copy(Z), copy(Z), copy(Z), copy(Z), copy(Z),
        zeros(Float64, dimension),
        copy(Z),
        [zeros(ComplexF64, dimension, dimension) for _ in 1:DIIS_size],
        [zeros(ComplexF64, dimension, dimension) for _ in 1:DIIS_size],
        1, 0,
         zeros(ComplexF64, NSHIFT)
    )
end


function Coulomb(k::Vector{Int},T1::Vector{Float64},T2::Vector{Float64},gatedis::Float64)::Float64
  
   return k==[0,0] ? gatedis : tanh(norm([T1 T2]*k*gatedis))/norm(k[1]*T1+k[2]*T2)
end

function Construct_DensityMatrix(work::HFWork,csr::ShiftCSR,wave::Vector{Vector{Int64}},wave_n1::Vector{Int64},wave_n2::Vector{Int64},
                               input_DensityMatrix::Matrix{ComplexF64},single_Ham::Matrix{ComplexF64},
                               energy_input::Float64,density::Float64,temp::Float64,Area::Float64,Coulomb_matrix::Matrix{Float64})
  
 
   dimension=length(wave)

 

    fill!(work.HartreeMatrix, 0)
    fill!(work.FockMatrix, 0)
   fill!(work.NewDensityMatrix, 0)


  


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
              
              
                  tmp=Coulomb_matrix[g1,g3]
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

    @inbounds for sid in 1:ix.NSHIFT
        acc = 0.0 + 0.0im
        lo = Int(offsets[sid])
        hi = Int(offsets[sid+1]) - 1
        for k in lo:hi
            g2 = Int(g2_list[k])
            g4 = Int(g3_list[k])
            acc += input_DensityMatrix[g4, g2]
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

    copy!(work.H_phys, single_Ham)
    work.H_phys .+= work.HartreeMatrix
    work.H_phys .-= work.FockMatrix

 
   

     
   
      old = BLAS.get_num_threads()
    
       BLAS.set_num_threads(min(8,Threads.nthreads()))   # or some smaller number like 4/8
       FFF = eigen(Hermitian(work.H_phys))        # or eigen(Hermitian(...), 1:filling)
    
      BLAS.set_num_threads(old)
   

    copy!(work.HF_eigenvalue, real(FFF.values))
    copy!(work.HF_eigenvector, FFF.vectors)

    fermi_level,_=find_FL(work.HF_eigenvalue,density,work.HF_eigenvalue[1]-10*temp,work.HF_eigenvalue[end]+10*temp,temp,Area,0.0)
    fermifactor=fermi_function.((work.HF_eigenvalue.-fermi_level)/temp)
    entropy=sum(s_function.((work.HF_eigenvalue.-fermi_level)/temp))

 
   
    mul!(work.NewDensityMatrix,
     work.HF_eigenvector * Diagonal(fermifactor),
     adjoint(work.HF_eigenvector))
    
    symmetrize_from_lower!(work.NewDensityMatrix)

  @. work.DeltaMatrix = work.NewDensityMatrix - input_DensityMatrix
       mixing=0.5
   @. work.output_DensityMatrix = mixing*input_DensityMatrix + (1-mixing)*work.NewDensityMatrix
    
  
  
   

  




  eout = sum(abs2, work.DeltaMatrix)
  
  

 

   energy =  real(dot(input_DensityMatrix, single_Ham)) + 0.5*real(dot(input_DensityMatrix, work.HartreeMatrix)) -0.5*real(dot(input_DensityMatrix, work.FockMatrix))
         
  
   freeenergy= energy-temp*entropy

   energy_change=real(energy-energy_input)



 return eout, energy_change, energy,freeenergy,fermi_level,fermifactor
end


function s_function(x::Float64)
    absx=abs(x)
   return log(1+exp(-absx))+absx/(1+exp(absx))

end



function fermi_function(x::Float64)
  if x>0
    return exp(-x)/(1+exp(-x))
  else
    return 1/(1+exp(x))

  end
end


function find_FL(quasi_particle_energy::Vector{Float64},target_density::Float64,val_s::Float64,val_e::Float64,temp::Float64,Area::Float64,bg_particle_density::Float64)
  try_FL=(val_s+val_e)/2

  
  if target_density==0.0
    stan=10^(-9)
  else
    stan=abs(10^(-8)*target_density)
  end
  #fermifactor=[1/(exp((quasi_particle_energy[ja]-try_FL)/temp)+1) for ja in eachindex(quasi_particle_energy)]
  fermifactor=fermi_function.((quasi_particle_energy.-try_FL)/temp)
 
  fl=sum(fermifactor)/Area-bg_particle_density



 if abs(fl-target_density)<stan
  
    return try_FL,fl
  elseif fl-target_density>=stan
 
    return find_FL(quasi_particle_energy,target_density,val_s, try_FL,temp,Area,bg_particle_density)
  elseif fl-target_density<=-stan
 
    return find_FL(quasi_particle_energy,target_density,try_FL,val_e,temp,Area,bg_particle_density)
   end
 

end



function attraction_po(k::Vector{Int},T1::Vector{Float64},T2::Vector{Float64},lpo::Float64)::Float64
   knorm=norm(k[1]*T1+k[2]*T2)
   return -exp(-knorm^2*lpo^2)
end


function iteration_loop(initial_DensityMatrix::Matrix{ComplexF64},
                       T1::Vector{Float64},T2::Vector{Float64},
                      wave::Vector{Vector{Int64}},single_Ham::Matrix{ComplexF64},
                     constq::Float64,gatedis::Float64,density::Float64,temp::Float64,Area::Float64,lpo::Float64,attstr::Float64)
    eout=1.0
    itcount=0
    bad_count=0
    energy=0.0
    freeenergy=0.0
    energy_change=0.0
    DIIS_size=5
    input_DensityMatrix=copy(initial_DensityMatrix)
  
    fermi_level=0.0
    fermifactor=zeros(Float64, length(wave))


    dimension=length(wave)
    wave_n1 = Vector{Int64}(undef, length(wave))
    wave_n2 = Vector{Int64}(undef, length(wave))
     for g in 1:length(wave)
        wave_n1[g] = Int(wave[g][1])
        wave_n2[g] = Int(wave[g][2])
    end
    Coulomb_matrix=[(Coulomb(wave[g1]-wave[g2],T1,T2,gatedis)*constq+attstr*attraction_po(wave[g1]-wave[g2],T1,T2,lpo)) for g1 in eachindex(wave), g2 in eachindex(wave)]
    
    wl = build_wave_lookup(wave)
    csr = build_shiftcsr(wl, wave_n1, wave_n2)
    
 
    work = HFWork(dimension, DIIS_size, csr.ix.NSHIFT)

  

  


     eout_hist = Float64[]
     PLATEAU_N = 10
     PLATEAU_FRAC = 0.10
     E_EPS = 1e-30

    diis_fire_once = false
    diis_cooldown = 0              # prevent immediate re-trigger after DIIS
     DIIS_COOLDOWN = 10 
   
  
   println(Threads.nthreads())
   while eout>10^(-16) || abs(energy_change)>10^(-9) || bad_count< DIIS_size+4
      if eout<1*10^(-16)
       bad_count+=1
      else
        bad_count=0
      end
     dmk_used = input_DensityMatrix

    if (itcount > 1000 && abs(eout) > 10)
                itcount = 0

                # forget DIIS/plateau state completely
                bad_count = 0
                energy = 0.0
                energy_change = 0.0
                empty!(eout_hist)
                diis_fire_once = false
                diis_cooldown = 0
                work.diis_head = 1
                work.diis_len  = 0

                # random restart DM (your style)
                A = randn(dimension, dimension) + im*randn(dimension, dimension)
                dmk = (A + A') * 0.01
    end
      
      tic=time()

      if (itcount>150 && abs(eout)>10) || (itcount>30 && abs(eout)<10^(-7)) || diis_fire_once

      
        dmk=implement_DIIS(work.DIIS_input_DensityMatrix,work.DIIS_input_DeltaMatrix,DIIS_size)
    
        dmk_used=dmk
      
      

        eout, energy_change, energy,freeenergy,fermi_level,fermifactor=Construct_DensityMatrix(work,csr,wave,wave_n1,wave_n2,
                                                                                dmk,single_Ham,                                                                           
                                                                               energy,density,temp,Area,Coulomb_matrix)
       
     
        println("using DIIS")
        

           if diis_fire_once
            diis_fire_once = false
            empty!(eout_hist)          # <-- yes: clear history after firing
            diis_cooldown = DIIS_COOLDOWN
           end
      else
  
      

       eout, energy_change, energy,freeenergy,fermi_level,fermifactor=Construct_DensityMatrix(work,csr,wave,wave_n1,wave_n2,
                                                                                input_DensityMatrix,single_Ham,                                                       
                                                                               energy,density,temp,Area,Coulomb_matrix)
   
        
       
         
     

      end
         work.diis_head, work.diis_len = diis_push!(
        work.DIIS_input_DensityMatrix,
        work.DIIS_input_DeltaMatrix,
        dmk_used,              # <--- this is the fix: store the DM that was used
        work.DeltaMatrix,
        work.diis_head,
        work.diis_len
          )

    
        copy!(input_DensityMatrix, work.output_DensityMatrix)


      itcount+=1
  
      toc=time()
      println(toc-tic,"eout=$eout","energy_change=$energy_change","itcount=$itcount")
      flush(stdout)


          # update cooldown
        if diis_cooldown > 0
            diis_cooldown -= 1
        end

        # update history
        push!(eout_hist, eout)  # keep sign; we'll use abs where needed
        if length(eout_hist) > PLATEAU_N
            popfirst!(eout_hist)
        end

        # plateau detection (only if not cooling down and not already scheduled)
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





    return work.DIIS_input_DensityMatrix,
       work.DIIS_input_DeltaMatrix,
       work.HF_eigenvalue,
       work.HF_eigenvector,
       energy,freeenergy, eout,
       work.HartreeMatrix,
       work.FockMatrix,fermi_level,fermifactor
  
end

















function implement_DIIS(DIIS_input_projector::Vector{Matrix{ComplexF64}},DIIS_input_DeltaMatrix::Vector{Matrix{ComplexF64}},DIIS_size::Int)



      Bmatrix=zeros(Float64,DIIS_size+1,DIIS_size+1)
      for ja in 1:DIIS_size
       Bmatrix[ja,DIIS_size+1]=1
       Bmatrix[DIIS_size+1,ja]=1
      end
  
      for ja in 1:DIIS_size,jb in 1:DIIS_size
   
             Bmatrix[ja,jb]+=real(dot(DIIS_input_DeltaMatrix[ja],DIIS_input_DeltaMatrix[jb]))
          
      end

      inB=safe_inverse(Bmatrix)
      if inB≠0
         onh=zeros(Float64,DIIS_size+1)
         onh[end]=1
         coeff=inB* onh
         #dmk=coeff[1]*(DIIS_input_projector[1])
         #for ja in 2:DIIS_size
          # dmk+=coeff[ja]*(DIIS_input_projector[ja])
         #end
        dmk = zero(DIIS_input_projector[1])           # alloc once
        @inbounds for ja in 1:DIIS_size
            BLAS.axpy!(coeff[ja], DIIS_input_projector[ja], dmk)  # dmk += coeff[ja] * projector[ja]
       
            BLAS.axpy!(coeff[ja], DIIS_input_DeltaMatrix[ja], dmk)  # dmk += coeff[ja] * projector[ja]
        end
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
          println("Matrix is singular, doing randomstart again.")
          return pinv(A, 10^(-8))  # Use pseudoinverse as an alternative
      else
          rethrow(e)  # If another error occurs, propagate it
      end
  end
end