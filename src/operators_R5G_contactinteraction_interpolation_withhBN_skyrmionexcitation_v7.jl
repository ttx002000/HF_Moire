using LinearAlgebra
using Arpack
using Combinatorics
using Random
BLAS.set_num_threads(1)

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







function triangle_initial_Densitymatrix(NL::Int,θ::Float64,gcutoff::Float64,
              uD::Float64,λ::Float64,enlarge_factor::Int,V0_hBN::Float64,V1_hBN::Float64,
              ψ_hBN::Float64,V2_scalar::Float64,ϕ::Float64,pin_coeff::Float64,ϵr::Float64,dedis::Float64,defec_pos::Int)
   
    aGr=0.246
    ϵ=0.2504/aGr-1 #This is the normal one
  
    G1=2π/aGr*[1,-1/√3]
    G2=2π/aGr*[0,2/√3]
    
    Rθ=[cos(θ) -sin(θ);sin(θ) cos(θ)]
    b1=(G1-(1+ϵ)^(-1)*Rθ*G1)/enlarge_factor
    b2=(G2-(1+ϵ)^(-1)*Rθ*G2)/enlarge_factor
    T1=b1
    T2=b2

    a1m=inv([b1';b2'])*[2π,0]
    a2m=inv([b1';b2'])*[0,2π]
    am=norm(a1m);
    
    Area=√3/2*am^2
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
    
    wave_dict=Dict{Vector{Int64},Int64}()
    for ja in eachindex(wave)
     wave_dict[wave[ja]]=ja
    end
 
    wave_diff=Vector{Int64}[]
    for ja in eachindex(wave), jb in eachindex(wave)
        push!(wave_diff,wave[ja]-wave[jb])
       
    end
    wave_diff=unique(wave_diff)

    dimension=length(wave)
   
    
    
    spinor_set=Vector{Vector{ComplexF64}}(undef,length(wave))
    for jb in eachindex(wave)
      hh=(1-λ)*get_Ham_Holomorphic([T1 T2]*(wave[jb]),uD,1,1,NL)+(λ)*get_Ham([T1 T2]*(wave[jb]),uD,1,1,NL)
      v1=eigen(hh).vectors[:,NL+1]
      if uD>0.0
        v1angle=angle(v1[2*NL])
      else
        throw("there is an error")
      end
      spinor_set[jb]=v1*exp(-im*v1angle)/norm(v1)
    end



     overlapmatrix=zeros(ComplexF64,length(wave),length(wave))
   
      for jb in eachindex(wave), jd in eachindex(wave)
       overlapmatrix[jb,jd]=spinor_set[jb]'*spinor_set[jd]
      end
 



    single_Ham=zeros(ComplexF64,dimension,dimension)
    single_MoirePo=zeros(ComplexF64,dimension,dimension)
     

    single_eigenvalue=zeros(Float64,length(wave))
    single_eigenvector=zeros(ComplexF64,length(wave),length(wave))

         ω=exp(im*2π/3)
       op_1=zeros(ComplexF64,2*NL,2*NL)
       op_1[1:2,1:2]=[1 1;ω ω]

       op_2=zeros(ComplexF64,2*NL,2*NL)
       op_2[1:2,1:2]=[1 ω^2;ω^2 ω]

       op_3=zeros(ComplexF64,2*NL,2*NL)
       op_3[1:2,1:2]=[1 ω;1 ω]

       
       op_4=zeros(ComplexF64,2*NL,2*NL)
       op_4[1:2,1:2]=[1 0; 0 1]

       op_5=Matrix{Float64}(I,2*NL,2*NL)
    
   
      
      
     
    
      for jb in eachindex(wave)
         
         
          single_Ham[jb,jb]=norm([T1 T2]*(wave[jb]))^2*200

         
      end


     
     for jc in eachindex(wave)
    
        pos=findfirst(item->item==wave[jc]-enlarge_factor*b1T,wave)
        if pos≠nothing
       
          single_MoirePo[jc,pos]+=V1_hBN*exp(-im*ψ_hBN)*(spinor_set[jc]'*op_1*spinor_set[pos])
          single_MoirePo[jc,pos]+=V2_scalar*exp(-im*ϕ)*(spinor_set[jc]'*op_5*spinor_set[pos])
          
        end
        
        pos=findfirst(item->item==wave[jc]-enlarge_factor*b2T,wave)
        if pos≠nothing

            single_MoirePo[jc,pos]+=V1_hBN*exp(-im*ψ_hBN)*(spinor_set[jc]'*op_2*spinor_set[pos])
            single_MoirePo[jc,pos]+=V2_scalar*exp(-im*ϕ)*(spinor_set[jc]'*op_5*spinor_set[pos])
        end


        pos=findfirst(item->item==wave[jc]+(b2T+b1T)*enlarge_factor,wave)
        if pos≠nothing
   
            single_MoirePo[jc,pos]+=V1_hBN*exp(-im*ψ_hBN)*(spinor_set[jc]'*op_3*spinor_set[pos])
            single_MoirePo[jc,pos]+=V2_scalar*exp(-im*ϕ)*(spinor_set[jc]'*op_5*spinor_set[pos])
        end
    
        single_MoirePo[jc,jc]+=V0_hBN/2*(spinor_set[jc]'*op_4*spinor_set[jc])
     end


    
     single_MoirePo=single_MoirePo+single_MoirePo'
     FFF=eigen(single_MoirePo+single_Ham)
    
      single_eigenvalue=real(FFF.values)
      single_eigenvector=FFF.vectors
    
     # more localized potential for small d
     
      pinning_po=zeros(ComplexF64,dimension,dimension)

      if defec_pos==1
         defec_pos_vec=(2*a2m-a1m)/enlarge_factor*2/3
      elseif defec_pos==2
         defec_pos_vec=(2*a2m-a1m)/enlarge_factor*1/3
      elseif defec_pos==3
        defec_pos_vec=(2*a2m-a1m)/enlarge_factor*0
      elseif defec_pos==4
        defec_pos_vec=(2*a2m-a1m)/enlarge_factor*(-2/3)
      end


      for jc in eachindex(wave)
        for jd in eachindex(wave_diff)
          if haskey(wave_dict,wave[jc]+wave_diff[jd])
              pos=wave_dict[wave[jc]+wave_diff[jd]]
              pinning_po[pos,jc]+=pin_coeff/(Area*ϵr)*get_Fourier_potential(wave_diff[jd],T1,T2,dedis)*(spinor_set[pos]'*op_5*spinor_set[jc])*exp(im*dot([T1 T2]*wave_diff[jd],defec_pos_vec))
          end
        end
      
      end

      pinning_po=(pinning_po+pinning_po')/2
     


    
   
 
    
     
      A=randn(ComplexF64,dimension,dimension)
      input_DensityMatrix=(A+A')*1.0
   

    

    
    return overlapmatrix, wave, input_DensityMatrix, single_MoirePo, pinning_po, single_Ham, single_eigenvalue,single_eigenvector, T1, T2, a1m, a2m, b1,b2,spinor_set,Area
      

       
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

    HF_eigenvalue_occ::Vector{Float64}        # length = filling
    HF_eigenvector_occ::Matrix{ComplexF64}

    DIIS_input_DensityMatrix::Vector{Matrix{ComplexF64}}
    DIIS_input_DeltaMatrix::Vector{Matrix{ComplexF64}}
    diis_head::Int
    diis_len::Int
    HartreeAccShift::Vector{ComplexF64}
end

function HFWork(dimension::Int, DIIS_size::Int, NSHIFT::Int,filling::Int)
    Z = zeros(ComplexF64, dimension, dimension)
    HFWork(
        copy(Z), copy(Z), copy(Z), copy(Z), copy(Z), copy(Z),
        zeros(Float64, dimension),
        copy(Z),
        zeros(Float64,filling+5),
        zeros(ComplexF64,dimension,filling+5),
        [zeros(ComplexF64, dimension, dimension) for _ in 1:DIIS_size],
        [zeros(ComplexF64, dimension, dimension) for _ in 1:DIIS_size],
        1, 0,
         zeros(ComplexF64, NSHIFT)
    )
end


function Coulomb(k::Vector{Int},T1::Vector{Float64},T2::Vector{Float64})::Float64
   D=25
   return k==[0,0] ? D*9047.5636 : tanh(norm([T1 T2]*k*D))/norm(k[1]*T1+k[2]*T2)*9047.5636
end

function Construct_DensityMatrix(work::HFWork,csr::ShiftCSR,wave::Vector{Vector{Int64}},wave_n1::Vector{Int64},wave_n2::Vector{Int64},
                               input_DensityMatrix::Matrix{ComplexF64},single_Ham::Matrix{ComplexF64},
                               single_MoirePo::Matrix{ComplexF64},pinning_po::Matrix{ComplexF64},overlapmatrix::Matrix{ComplexF64},
                               energy_input::Float64,filling::Int,Area::Float64,Coulomb_matrix::Matrix{ComplexF64})
  
 
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


    #=
    tic=time()
    
    Threads.@threads :greedy  for g1 in eachindex(wave)
        @inbounds begin
            w1n1=wave_n1[g1]; 
            w1n2=wave_n2[g1]
        for g3 in 1:g1
        
            w3n1 = wave_n1[g3];
            w3n2 = wave_n2[g3]
            
            dg_n1=w1n1- w3n1
            dg_n2=w1n2- w3n2
            lo, hi = csr_range(csr,  dg_n1,  dg_n2)
        
            acc = 0.0 + 0.0im
            @inbounds for k in lo:hi
        
            g2 = Int(csr.g2_list[k])
            g4 = Int(csr.g3_list[k])
            
                acc +=input_DensityMatrix[g4, g2] * overlapmatrix[g2, g4]
            
            end
            val = acc * Coulomb_matrix[g1,g3]
        
            work.HartreeMatrix[g1,g3] = val
        
        end
        end
    end
    toc=time()
    println(toc-tic,"Hartree")
    =#
    
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

 
   

     
   
      old = BLAS.get_num_threads()
    
       BLAS.set_num_threads(min(8,Threads.nthreads()))   # or some smaller number like 4/8
       FFF = eigen(Hermitian(work.H_phys),1:filling+5)        # or eigen(Hermitian(...), 1:filling)
    
      BLAS.set_num_threads(old)
   

    copy!(work.HF_eigenvalue_occ, real(FFF.values))
    copy!(work.HF_eigenvector_occ, FFF.vectors)

 
     @views Vocc = work.HF_eigenvector_occ[:, 1:filling]
    mul!(work.NewDensityMatrix, Vocc, adjoint(Vocc))
    symmetrize_from_lower!(work.NewDensityMatrix)

  @. work.DeltaMatrix = work.NewDensityMatrix - input_DensityMatrix
       mixing=0.5
   @. work.output_DensityMatrix = mixing*input_DensityMatrix + (1-mixing)*work.NewDensityMatrix
    
  
  
   

  




  eout = sum(abs2, work.DeltaMatrix)
  
  

 

   energy = real(dot(input_DensityMatrix, single_MoirePo)) +real(dot(input_DensityMatrix, single_Ham)) +  real(dot(input_DensityMatrix, pinning_po))+ 0.5*real(dot(input_DensityMatrix, work.HartreeMatrix)) -0.5*real(dot(input_DensityMatrix, work.FockMatrix))
         

   energy_change=real(energy-energy_input)



 return eout, energy_change, energy
end








function pick_random_jld2_path(folder::AbstractString)
    isdir(folder) || return nothing

    jld2_files = filter(readdir(folder; join=true)) do p
        isfile(p) && endswith(lowercase(p), ".jld2")
    end

    isempty(jld2_files) && return nothing
    return rand(jld2_files)
end



function transform_dm(dm_old::Matrix{ComplexF64},wave::Vector{Vector{Int}},seed_spinor_set::Vector{Vector{ComplexF64}},spinor_set::Vector{Vector{ComplexF64}})
    dm_new=zeros(ComplexF64,length(wave),length(wave))
    for ja in eachindex(wave), jb in eachindex(wave)
       dm_new[ja,jb]=(spinor_set[ja]'*seed_spinor_set[ja])*dm_old[ja,jb]*(seed_spinor_set[jb]'*spinor_set[jb])
   end

   return dm_new
end






function iteration_loop(initial_DensityMatrix::Matrix{ComplexF64},
                       T1::Vector{Float64},T2::Vector{Float64},
                      wave::Vector{Vector{Int64}},single_Ham::Matrix{ComplexF64},
                       single_MoirePo::Matrix{ComplexF64},pinning_po::Matrix{ComplexF64},constq::Float64,ϵr::Float64,overlapmatrix::Matrix{ComplexF64},filling::Int,Area::Float64)
    eout=1.0
    itcount=0
    bad_count=0
    energy=0.0
    energy_change=0.0
    DIIS_size=5
    input_DensityMatrix=copy(initial_DensityMatrix)
  
    dimension=length(wave)
    wave_n1 = Vector{Int64}(undef, length(wave))
    wave_n2 = Vector{Int64}(undef, length(wave))
     for g in 1:length(wave)
        wave_n1[g] = Int(wave[g][1])
        wave_n2[g] = Int(wave[g][2])
    end
    Coulomb_matrix=[(Coulomb(wave[g1]-wave[g2],T1,T2)/ϵr+constq)*overlapmatrix[g1,g2] for g1 in eachindex(wave), g2 in eachindex(wave)]
    wl = build_wave_lookup(wave)
    csr = build_shiftcsr(wl, wave_n1, wave_n2)
    
 
    work = HFWork(dimension, DIIS_size, csr.ix.NSHIFT,filling)

  

  


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
                input_DensityMatrix = (A + A') * 0.01
                dmk_used=input_DensityMatrix
                println("random start again")
    end
      
      tic=time()

      if (itcount>150 && abs(eout)>10) || (itcount>30 && abs(eout)<10^(-7)) || diis_fire_once

      
        dmk=implement_DIIS(work.DIIS_input_DensityMatrix,work.DIIS_input_DeltaMatrix,DIIS_size)
    
        dmk_used=dmk
      
      

        eout, energy_change, energy=Construct_DensityMatrix(work,csr,wave,wave_n1,wave_n2,
                                                                                dmk,single_Ham,
                                                                                single_MoirePo, pinning_po,overlapmatrix,
                                                                               energy,filling,Area,Coulomb_matrix)
       
     
        println("using DIIS")
        

           if diis_fire_once
            diis_fire_once = false
            empty!(eout_hist)          # <-- yes: clear history after firing
            diis_cooldown = DIIS_COOLDOWN
           end
      else
  
      

       eout, energy_change, energy=Construct_DensityMatrix(work,csr,wave,wave_n1,wave_n2,
                                                                                input_DensityMatrix,single_Ham,
                                                                                single_MoirePo, pinning_po,overlapmatrix,
                                                                               energy,filling,Area,Coulomb_matrix)
   
        
       
         
     

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

   Ffull = eigen(Hermitian(work.H_phys))   # or eigen(Hermitian(work.H_phys)) if you prefer non-mutating
   copy!(work.HF_eigenvalue, real(Ffull.values)) 
   copy!(work.HF_eigenvector, Ffull.vectors)




    return work.DIIS_input_DensityMatrix,
       work.DIIS_input_DeltaMatrix,
       work.HF_eigenvalue,
       work.HF_eigenvector,
       energy, eout,
       work.HartreeMatrix,
       work.FockMatrix
  
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