using LinearAlgebra
using Arpack
using Combinatorics
using Random


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
    offsets = Vector{Int64}(undef, NSHIFT + 1)
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
                    g2_list[p] = Int32(g2)
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
              ψ_hBN::Float64,V2_scalar::Float64,ϕ::Float64)
   
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
          hh=(1-λ)*get_Ham_Holomorphic([T1 T2]*(wave[jb]),uD,1,1,NL)+(λ)*get_Ham([T1 T2]*(wave[jb]),uD,1,1,NL)
         
         
          single_Ham[jb,jb]=real(eigen(hh).values[NL+1])

         
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
    
    


    
   
 
    
     
      A=randn(ComplexF64,dimension,dimension)
      input_DensityMatrix=(A+A')*1.0
    
      if filepos==6 && rand()>0.5
        
      scratch_dir = ENV["SCRATCH"]
        
        seed_path=joinpath(scratch_dir,"triangle_R5G_contact_interpolation_withhBN_skyrmionexcitation_v4/data_output$(Int(args[16]))/seed/$(args[2])theta0.0constq$(args[10])enlarge$(args[7])cutoff$(args[6])filling_hBNMoire.jld2")
          st=load(seed_path)
          
          input_DensityMatrix=st["densitymatrix"]
          println("teaking seed")
      end

    

    
    return overlapmatrix, wave, input_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,single_eigenvector, T1, T2, a1m, a2m, b1,b2,spinor_set
      

       
end










function Coulomb(k::Vector{Int},T1::Vector{Float64},T2::Vector{Float64})::Float64
   D=25
   return k==[0,0] ? D*9047.5636 : tanh(norm([T1 T2]*k*D))/norm(k[1]*T1+k[2]*T2)*9047.5636
end

function Construct_DensityMatrix(csr::ShiftCSR,wave::Vector{Vector{Int64}},wave_n1::Vector{Int64},wave_n2::Vector{Int64},
                               input_DensityMatrix::Matrix{ComplexF64},single_Ham::Matrix{ComplexF64},
                               single_MoirePo::Matrix{ComplexF64},overlapmatrix::Matrix{ComplexF64},
                               energy_input::Float64,filling::Int,Area::Float64,Coulomb_matrix::Matrix{ComplexF64},
                               whether_DIIS::Int64,DIIS_H::Matrix{ComplexF64})
  
 
   dimension=length(wave)
  HartreeMatrix=zeros(ComplexF64,length(wave),length(wave))
  FockMatrix=zeros(ComplexF64,length(wave),length(wave))
  output_DensityMatrix=zeros(ComplexF64,dimension,dimension)
  DeltaMatrix=zeros(ComplexF64,dimension,dimension)
  NewDensityMatrix=zeros(ComplexF64,dimension,dimension)
  HF_eigenvalue=zeros(Float64,dimension)
  HF_eigenvector=zeros(ComplexF64,dimension,dimension)

  input_DM=input_DensityMatrix
 

  Threads.@threads :greedy for g1 in eachindex(wave)
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
                  acc += input_DM[g3,g2] * tmp
                 
              end
                
          

              FockMatrix[g1,g4] = acc
        
              
          
  

      end
    end
  end


  Threads.@threads :greedy for g1 in eachindex(wave)
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
         
            acc += input_DM[g4, g2] * overlapmatrix[g2, g4]
        
        end
        val = acc * Coulomb_matrix[g1,g3]
     
        HartreeMatrix[g1,g3] = val
      
     end
    end
  end


  FockMatrix=reshape(FockMatrix,dimension,dimension)
  HartreeMatrix=reshape(HartreeMatrix,dimension,dimension)


    HartreeMatrix=(HartreeMatrix+HartreeMatrix'-real(Diagonal(HartreeMatrix)))/Area
    FockMatrix=(FockMatrix+FockMatrix'-real(Diagonal(FockMatrix)))/Area
 
    H_phys = single_MoirePo + single_Ham + HartreeMatrix - FockMatrix

    #H_diag = (whether_DIIS==1 ? DIIS_H : H_phys)


   #H = single_MoirePo + single_Ham + HartreeMatrix - FockMatrix
   FFF = eigen!(H_phys)
   #FFF = eigen!(Hermitian(H_diag))
   HF_eigenvalue=real(FFF.values)
   HF_eigenvector=FFF.vectors

   #DeltaMatrix=H_phys*input_DensityMatrix-input_DensityMatrix*H_phys
 
      NewDensityMatrix=HF_eigenvector[:,1:filling]*HF_eigenvector[:,1:filling]'
      NewDensityMatrix=0.5*(NewDensityMatrix'+ NewDensityMatrix)
       DeltaMatrix=NewDensityMatrix-input_DensityMatrix
       mix_ratio=0.5
       output_DensityMatrix=mix_ratio*input_DensityMatrix+(1-mix_ratio)*NewDensityMatrix

  
  
   

  




  eout=sum(abs2,DeltaMatrix)
  
  

 
  ss=single_MoirePo+single_Ham+0.5*HartreeMatrix-0.5*FockMatrix
  energy = real(sum(ss .* transpose(input_DensityMatrix)))


   energy_change=real(energy-energy_input)



 return  eout,energy_change,output_DensityMatrix,DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,HartreeMatrix,FockMatrix,H_phys
end



function get_manybodyoverlap(args)

    NL=args[1]
    θ=args[2]/180*pi;
    constq=args[3]
    ϵr=args[4]
    uD=args[5]
    filling=Int(args[6])
    gcutoff=args[7]
    λ=args[8]
    trytimes=Int(args[9])
    enlarge_factor=Int(args[10])
    V0_hBN=args[11]
    V1_hBN=args[12]
    ψ_hBN=args[13]
    V2_scalar=args[14]
    ϕ=args[15]/180*π
    filepos=Int(args[16])


     scratch_dir = ENV["SCRATCH"]
     seed_path=joinpath(scratch_dir, "triangle_R5G_contact_interpolation_withhBN_skyrmionexcitation_v4/data_output$(Int(args[16]))/dispersion/seed/$(args[1])NL$(args[2])theta$(args[3])constq$(args[4])ϵr$(args[5])uD$(args[6])filling$(args[7])cutoff$(args[8])lambda$(args[9])trytime$(args[10])enlarge$(args[11])V0_hBN$(args[12])V1_hBN$(args[13])ψ_hBN$(args[14])V2_scalar$(args[15])ϕ.jld2")




    st=load(seed_path)
    densitymatrix=st["densitymatrix"]
    a1m=st["a1m"]
    a2m=st["a2m"]
    wave=st["wave"]
    T1=st["T1"]
    T2=st["T2"]
    single_Ham=st["single_Ham"]
    single_MoirePo=st["single_MoirePo"]
    HF_eigenvector=st["HF_eigenvector"]
    spinor_set=st["spinor_set"]
    b1T=[1,0]
    b2T=[0,1]
    b1=[T1 T2]*b1T
    b2=[T1 T2]*b2T
    Area=√3/2*norm(a1m)^2

    overlapmatrix=zeros(ComplexF64,length(wave),length(wave))
    
    for jb in eachindex(wave), jd in eachindex(wave)
    overlapmatrix[jb,jd]=spinor_set[jb]'*spinor_set[jd]
    end

    wave_n1 = Vector{Int64}(undef, length(wave))
    wave_n2 = Vector{Int64}(undef, length(wave))
    for g in eachindex(wave)
        wave_n1[g] = Int32(wave[g][1])
        wave_n2[g] = Int32(wave[g][2])
    end
    Coulomb_matrix=[(Coulomb(wave[g1]-wave[g2],T1,T2)/ϵr+constq)*overlapmatrix[g1,g2] for g1 in eachindex(wave), g2 in eachindex(wave)]
    wl = build_wave_lookup(wave)
    csr = build_shiftcsr(wl, wave_n1, wave_n2)


    shift_set=vec([a1m/enlarge_factor*ja+a2m/enlarge_factor*jb for ja in 0:Int(enlarge_factor)-1, jb in 0:Int(enlarge_factor)-1])
    shift_Hartree_Fockvector=Vector{Matrix{ComplexF64}}(undef,length(shift_set))
    for ja in eachindex(shift_set)
        hh=zeros(ComplexF64,length(wave),filling)
        for jb in eachindex(wave)
        hh[jb,:]=HF_eigenvector[jb,1:filling]*exp(-im*dot([T1 T2]*wave[jb],shift_set[ja])) #I think should be a minus sign, but double check
        end

        shift_Hartree_Fockvector[ja]=copy(hh)
    end



    H_matrixelement=zeros(ComplexF64,length(shift_set),length(shift_set))
    manybody_overlap=zeros(ComplexF64,length(shift_set),length(shift_set))


    for ja in eachindex(shift_set),jb in eachindex(shift_set)
        Smatrix=shift_Hartree_Fockvector[ja]'*shift_Hartree_Fockvector[jb]
        manybody_overlap[ja,jb]= det(Smatrix)
    end


        for ja in eachindex(shift_set),jb in 1:ja
            println(ja,jb)
            flush(stdout)
                Smatrix=shift_Hartree_Fockvector[ja]'*shift_Hartree_Fockvector[jb]
                X = Smatrix \ (shift_Hartree_Fockvector[ja]')     # (filling × Ng)
                Projector = shift_Hartree_Fockvector[jb] * X 
                #Projector=shift_Hartree_Fockvector[jb]*inv(Smatrix)*shift_Hartree_Fockvector[ja]'

                H_matrixelement[ja,jb]=manybody_overlap[ja,jb]* Construct_DispersionMatrixelement(csr,wave,wave_n1,wave_n2,
                                            Projector,single_Ham,
                                        single_MoirePo,overlapmatrix,
                                        Area,Coulomb_matrix)
                H_matrixelement[jb,ja]= conj(H_matrixelement[ja,jb])
       end


    return manybody_overlap,H_matrixelement,shift_set

end



function Construct_DispersionMatrixelement(csr::ShiftCSR,wave::Vector{Vector{Int64}},wave_n1::Vector{Int64},wave_n2::Vector{Int64},
                               input_DensityMatrix::Matrix{ComplexF64},single_Ham::Matrix{ComplexF64},
                               single_MoirePo::Matrix{ComplexF64},overlapmatrix::Matrix{ComplexF64},
                               Area::Float64,Coulomb_matrix::Matrix{ComplexF64})
  
 
   dimension=length(wave)
  HartreeMatrix=zeros(ComplexF64,length(wave),length(wave))
  FockMatrix=zeros(ComplexF64,length(wave),length(wave))
  input_DM=input_DensityMatrix
 

  Threads.@threads for g1 in eachindex(wave)
    @inbounds begin
    
        w1n1=wave_n1[g1]; 
        w1n2=wave_n2[g1]
      for g4 in eachindex(wave)
          
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
                  acc += input_DM[g3,g2] * tmp
                 
              end
                
          

              FockMatrix[g1,g4] = acc
        
              
          
  

      end
    end
  end


  Threads.@threads for g1 in eachindex(wave)
    @inbounds begin
         w1n1=wave_n1[g1]; 
         w1n2=wave_n2[g1]
     for g3 in eachindex(wave)
     
        w3n1 = wave_n1[g3];
        w3n2 = wave_n2[g3]
         
        dg_n1=w1n1- w3n1
        dg_n2=w1n2- w3n2
        lo, hi = csr_range(csr,  dg_n1,  dg_n2)
     
        acc = 0.0 + 0.0im
        @inbounds for k in lo:hi
       
        g2 = Int(csr.g2_list[k])
        g4 = Int(csr.g3_list[k])
         
            acc += input_DM[g4, g2] * overlapmatrix[g2, g4]
        
        end
        val = acc * Coulomb_matrix[g1,g3]
     
        HartreeMatrix[g1,g3] = val
      
     end
    end
  end





   
    H_phys = single_MoirePo + single_Ham + 0.5*(HartreeMatrix/Area - FockMatrix/Area)


  Hmatrix_element =(sum(H_phys .* transpose(input_DensityMatrix)))


   



 return  Hmatrix_element 
end








function iteration_loop(initial_DensityMatrix::Matrix{ComplexF64},
                       T1::Vector{Float64},T2::Vector{Float64},
                      wave::Vector{Vector{Int64}},single_Ham::Matrix{ComplexF64},
                       single_MoirePo::Matrix{ComplexF64},constq::Float64,ϵr::Float64,overlapmatrix::Matrix{ComplexF64},filling::Int,Area::Float64)
    eout=1.0
    itcount=0
  
    dimension=length(wave)
    wave_n1 = Vector{Int64}(undef, length(wave))
    wave_n2 = Vector{Int64}(undef, length(wave))
     for g in 1:length(wave)
        wave_n1[g] = Int32(wave[g][1])
        wave_n2[g] = Int32(wave[g][2])
    end
    Coulomb_matrix=[(Coulomb(wave[g1]-wave[g2],T1,T2)/ϵr+constq)*overlapmatrix[g1,g2] for g1 in eachindex(wave), g2 in eachindex(wave)]
    wl = build_wave_lookup(wave)
    csr = build_shiftcsr(wl, wave_n1, wave_n2)
    
    HF_eigenvalue=zeros(Float64,dimension)
    HF_eigenvector=zeros(ComplexF64,dimension,dimension)
    HartreeMatrix=zeros(ComplexF64,dimension,dimension)
    FockMatrix=zeros(ComplexF64,dimension,dimension)
  
    DIIS_size=5
    DIIS_input_DensityMatrix=Vector{Matrix{ComplexF64}}(undef,DIIS_size)
    DIIS_input_DeltaMatrix=Vector{Matrix{ComplexF64}}(undef,DIIS_size)
    DIIS_output_HFHam=Vector{Matrix{ComplexF64}}(undef,DIIS_size)
    input_DensityMatrix=initial_DensityMatrix
    bad_count=0
    energy=0.0
    energy_change=0.0

     eout_hist = Float64[]
     PLATEAU_N = 10
     PLATEAU_FRAC = 0.10
     E_EPS = 1e-30

    diis_fire_once = false
    diis_cooldown = 0              # prevent immediate re-trigger after DIIS
     DIIS_COOLDOWN = 10 
   
  
   println(Threads.nthreads())
   while (eout>1*10^(-16)) || (bad_count<DIIS_size) || (abs(energy_change)>1*10^(-9))
      if eout<1*10^(-16)
       bad_count+=1
      else
        bad_count=0
      end
      
      tic=time()

      if (itcount>150 && abs(eout)>10) || (itcount>30 && abs(eout)<10^(-6)) || diis_fire_once

      
        dmk=implement_DIIS(DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,DIIS_size)
        if (itcount>1000 && abs(eout)>10)
           itcount=0
           
              A=randn(dimension,dimension)+im*randn(dimension,dimension)
              dmk=(A+A')*0.01
         
        end

          #fk= implement_DIIS(DIIS_output_HFHam,DIIS_input_DeltaMatrix,DIIS_size)
      

          eout,energy_change,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,DIIS_size)+1],HF_eigenvalue,HF_eigenvector,energy,HartreeMatrix,FockMatrix,DIIS_output_HFHam[mod(itcount,DIIS_size)+1]=Construct_DensityMatrix(csr,wave,wave_n1,wave_n2,
                                                                                                                                                                             dmk,single_Ham,
                                                                                                                                                                            single_MoirePo,overlapmatrix,
                                                                                                                                                                            energy,filling,Area,Coulomb_matrix,0,zeros(ComplexF64,dimension,dimension))
       
        DIIS_input_DensityMatrix[mod(itcount,DIIS_size)+1]=dmk
        #DIIS_input_DensityMatrix[mod(itcount,DIIS_size)+1]=input_DensityMatrix
        input_DensityMatrix=output_DensityMatrix
        println("using DIIS")
        

           if diis_fire_once
            diis_fire_once = false
            empty!(eout_hist)          # <-- yes: clear history after firing
            diis_cooldown = DIIS_COOLDOWN
           end
      else
  
      

        eout,energy_change,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,DIIS_size)+1],HF_eigenvalue,HF_eigenvector,energy,HartreeMatrix,FockMatrix,DIIS_output_HFHam[mod(itcount,DIIS_size)+1]=Construct_DensityMatrix(csr,wave,wave_n1,wave_n2,
                                                                                                                                                                            input_DensityMatrix,single_Ham,
                                                                                                                                                                            single_MoirePo,overlapmatrix,
                                                                                                                                                                            energy,filling,Area,Coulomb_matrix,0,zeros(ComplexF64,dimension,dimension))
        DIIS_input_DensityMatrix[mod(itcount,DIIS_size)+1]=input_DensityMatrix
        input_DensityMatrix=output_DensityMatrix
        
       
         
     

      end

    



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







  return DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,eout,HartreeMatrix,FockMatrix

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