using LinearAlgebra
using Arpack
using Combinatorics
using Random


# The goal of this piece of code is to do the calculation for specifically Nq=1 case.
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
    pos::Matrix{Int32}   # 0 means “missing”
    n1min::Int
    n2min::Int
end

@inline function lookup(wl::WaveLookup, n1::Int, n2::Int)::Int64
    i1 = n1 - wl.n1min + 1
    i2 = n2 - wl.n2min + 1
    # bounds check + read
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






function triangle_initial_Densitymatrix(NL::Int,θ::Float64,gcutoff::Float64,
              uD::Float64,λ::Float64,enlarge_factor::Int,V0_hBN::Float64,V1_hBN::Float64,ψ_hBN::Float64,V2_scalar::Float64,ϕ::Float64)
   
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
    cutoff=18
    cutoffstandard=gcutoff*norm(b1)
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        gtest=ja*b1+jb*b2;
        if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
            push!(wave,ja*b1T+jb*b2T)
        end
    end
    num_spin=2
    dimension=num_spin*length(wave)
   
    
    
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
     
    single_Ham_spinless=zeros(ComplexF64,length(wave),length(wave))
    single_MoirePo_spinless=zeros(ComplexF64,length(wave),length(wave))
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
         
         
          single_Ham_spinless[jb,jb]=real(eigen(hh).values[NL+1])

         
      end


     
     for jc in eachindex(wave)
    
        pos=findfirst(item->item==wave[jc]-enlarge_factor*b1T,wave)
        if pos≠nothing
       
          single_MoirePo_spinless[jc,pos]+=V1_hBN*exp(-im*ψ_hBN)*(spinor_set[jc]'*op_1*spinor_set[pos])
          single_MoirePo_spinless[jc,pos]+=V2_scalar*exp(-im*ϕ)*(spinor_set[jc]'*op_5*spinor_set[pos])
          
        end
        
        pos=findfirst(item->item==wave[jc]-enlarge_factor*b2T,wave)
        if pos≠nothing

            single_MoirePo_spinless[jc,pos]+=V1_hBN*exp(-im*ψ_hBN)*(spinor_set[jc]'*op_2*spinor_set[pos])
            single_MoirePo_spinless[jc,pos]+=V2_scalar*exp(-im*ϕ)*(spinor_set[jc]'*op_5*spinor_set[pos])
        end


        pos=findfirst(item->item==wave[jc]+(b2T+b1T)*enlarge_factor,wave)
        if pos≠nothing
   
            single_MoirePo_spinless[jc,pos]+=V1_hBN*exp(-im*ψ_hBN)*(spinor_set[jc]'*op_3*spinor_set[pos])
            single_MoirePo_spinless[jc,pos]+=V2_scalar*exp(-im*ϕ)*(spinor_set[jc]'*op_5*spinor_set[pos])
        end
    
        single_MoirePo_spinless[jc,jc]+=V0_hBN/2*(spinor_set[jc]'*op_4*spinor_set[jc])
     end


    
     single_MoirePo_spinless=single_MoirePo_spinless+single_MoirePo_spinless'
     FFF=eigen(single_MoirePo_spinless+single_Ham_spinless)
    
      single_eigenvalue=real(FFF.values)
      single_eigenvector=FFF.vectors
    
    


    
      dd=zeros(ComplexF64,num_spin,length(wave),num_spin,length(wave))
      dd[1,:,1,:]=single_Ham_spinless
      dd[2,:,2,:]=single_Ham_spinless
      single_Ham=reshape(dd,dimension,dimension)

      ee=zeros(ComplexF64,num_spin,length(wave),num_spin,length(wave))
      ee[1,:,1,:]=single_MoirePo_spinless
      ee[2,:,2,:]=single_MoirePo_spinless
      single_MoirePo=reshape(ee,dimension,dimension)
    


   
 
    
     
      A=randn(ComplexF64,dimension,dimension)
      input_DensityMatrix=(A+A')*1.0
   

    

    
    return overlapmatrix, wave, input_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,single_eigenvector, T1, T2, a1m, a2m, b1,b2,spinor_set
      

       
end










function Coulomb(k::Vector{Int},T1::Vector{Float64},T2::Vector{Float64})::Float64
   D=25
   return k==[0,0] ? D*9047.5636 : tanh(norm([T1 T2]*k*D))/norm(k[1]*T1+k[2]*T2)*9047.5636
end

function Construct_DensityMatrix(wl::WaveLookup,wave::Vector{Vector{Int64}},wave_n1::Vector{Int64},wave_n2::Vector{Int64},
                               input_DensityMatrix::Matrix{ComplexF64},single_Ham::Matrix{ComplexF64},
                               single_MoirePo::Matrix{ComplexF64},overlapmatrix::Matrix{ComplexF64},
                               energy_input::Float64,filling::Int,Area::Float64,Coulomb_matrix::Matrix{ComplexF64})
  
    num_spin=2
   dimension=num_spin*length(wave)
  HartreeMatrix=zeros(ComplexF64,num_spin,length(wave),num_spin,length(wave))
  FockMatrix=zeros(ComplexF64,num_spin,length(wave),num_spin,length(wave))
  output_DensityMatrix=zeros(ComplexF64,dimension,dimension)
  DeltaMatrix=zeros(ComplexF64,dimension,dimension)
  NewDensityMatrix=zeros(ComplexF64,dimension,dimension)
  HF_eigenvalue=zeros(Float64,dimension)
  HF_eigenvector=zeros(ComplexF64,dimension,dimension)

  input_DM_reshaped=reshape(input_DensityMatrix,num_spin,length(wave),num_spin,length(wave))
  input_DM_traced=input_DM_reshaped[1,:,1,:]+input_DM_reshaped[2,:,2,:]

  tic=time()
  Threads.@threads for g1 in eachindex(wave)
    @inbounds begin
      #w1 = wave[g1]; 
      w1n1=wave_n1[g1]; 
      w1n2=wave_n2[g1]
    for g4 in 1:g1
          #w4 = wave[g4]; 
          w4n1=wave_n1[g4]; 
          w4n2=wave_n2[g4]
          acc11 = 0.0 + 0.0im
          acc12 = 0.0 + 0.0im
          acc21 = 0.0 + 0.0im
          acc22 = 0.0 + 0.0im

        if g1≠g4
            for g2 in eachindex(wave)
             # w2 = wave[g2]
             
              g3=lookup(wl,w1n1 + wave_n1[g2] - w4n1, w1n2 + wave_n2[g2] - w4n2)
              if g3≠0
                #tmp=Coulomb_matrix[g1,g3]*overlapmatrix[g1,g3]*overlapmatrix[g2,g4]
                tmp=Coulomb_matrix[g1,g3]*overlapmatrix[g2,g4]
                acc11 += input_DM_reshaped[1,g3,1,g2] * tmp
                acc12 += input_DM_reshaped[1,g3,2,g2] * tmp
                acc21 += input_DM_reshaped[2,g3,1,g2] * tmp
                acc22 += input_DM_reshaped[2,g3,2,g2] * tmp
            
              end
            end

            FockMatrix[1,g1,1,g4] = acc11
            FockMatrix[1,g1,2,g4] = acc12
            FockMatrix[2,g1,1,g4] = acc21
            FockMatrix[2,g1,2,g4] = acc22
            
        else
          for g2 in eachindex(wave)
              #w2 = wave[g2]
             
              g3=lookup(wl,w1n1 +wave_n1[g2]  - w4n1, w1n2 + wave_n2[g2] - w4n2)
              if g3≠0
                #tmp=Coulomb_matrix[g1,g3]*overlapmatrix[g1,g3]*overlapmatrix[g2,g4]
                tmp=Coulomb_matrix[g1,g3]*overlapmatrix[g2,g4]
                acc11 += input_DM_reshaped[1,g3,1,g2] * tmp
                acc21 += input_DM_reshaped[2,g3,1,g2] * tmp
                acc22 += input_DM_reshaped[2,g3,2,g2] * tmp
            
              end
            end

            FockMatrix[1,g1,1,g4] = acc11
            FockMatrix[2,g1,1,g4] = acc21
            FockMatrix[2,g1,2,g4] = acc22
        end
 

    end
   end
  end
  toc=time()
  println(toc-tic,"focktime")
  tic=time()

  Threads.@threads for g1 in eachindex(wave)
    @inbounds begin
        
        w1n1 = wave_n1[g1]; 
        w1n2 = wave_n2[g1]
     for g3 in 1:g1
        
        w3n1 = wave_n1[g3]; 
        w3n2 = wave_n2[g3]
        #tmp=Coulomb_matrix[g1,g3]*overlapmatrix[g1,g3]
        acc = 0.0 + 0.0im
        for g2 in eachindex(wave)
       
          g4 = lookup(wl,
                        w1n1 + wave_n1[g2] - w3n1,
                        w1n2 + wave_n2[g2] - w3n2)
          if g4≠0
            acc += input_DM_traced[g4, g2] * overlapmatrix[g2, g4]
          end
        end
        val = acc * Coulomb_matrix[g1,g3]
        HartreeMatrix[1, g1, 1, g3] = val
        HartreeMatrix[2, g1, 2, g3] = val
     end
    end
  end

 toc=time()
  println(toc-tic,"Hartreetime")

  tic=time()
  FockMatrix=reshape(FockMatrix,dimension,dimension)
  HartreeMatrix=reshape(HartreeMatrix,dimension,dimension)


    HartreeMatrix=(HartreeMatrix+HartreeMatrix'-real(Diagonal(HartreeMatrix)))/Area
    FockMatrix=(FockMatrix+FockMatrix'-real(Diagonal(FockMatrix)))/Area

   toc=time() 
     println(toc-tic,"the rest 1")

     tic=time()

   H = single_MoirePo + single_Ham + HartreeMatrix - FockMatrix
   FFF = eigen!(Hermitian(H))
   HF_eigenvalue=real(FFF.values)
   HF_eigenvector=FFF.vectors
   toc=time()
    println(toc-tic,"the rest 2")

       tic=time()
       NewDensityMatrix=HF_eigenvector[:,1:filling]*HF_eigenvector[:,1:filling]'
        NewDensityMatrix=0.5*(NewDensityMatrix'+ NewDensityMatrix)
       DeltaMatrix=NewDensityMatrix-input_DensityMatrix
       mix_ratio=0.5
       output_DensityMatrix=mix_ratio*input_DensityMatrix+(1-mix_ratio)*NewDensityMatrix


   

  




  eout=sum(abs2,DeltaMatrix)
  
  

 
    ss=single_MoirePo+single_Ham+0.5*HartreeMatrix-0.5*FockMatrix
    energy = real(sum(ss .* transpose(input_DensityMatrix)))


   energy_change=real(energy-energy_input)
   toc=time()
  println(toc-tic,"the rest 3")

 return  eout,energy_change,output_DensityMatrix,DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,HartreeMatrix,FockMatrix
end





function iteration_loop(initial_DensityMatrix::Matrix{ComplexF64},
                       T1::Vector{Float64},T2::Vector{Float64},
                      wave::Vector{Vector{Int64}},single_Ham::Matrix{ComplexF64},
                       single_MoirePo::Matrix{ComplexF64},constq::Float64,ϵr::Float64,overlapmatrix::Matrix{ComplexF64},filling::Int,Area::Float64)
    eout=1.0
    itcount=0
    num_spin=2
    dimension=num_spin*length(wave)
    wave_n1 = Vector{Int64}(undef, length(wave))
    wave_n2 = Vector{Int64}(undef, length(wave))
     for g in 1:length(wave)
        wave_n1[g] = Int32(wave[g][1])
        wave_n2[g] = Int32(wave[g][2])
    end
    Coulomb_matrix=[(Coulomb(wave[g1]-wave[g2],T1,T2)/ϵr+constq)*overlapmatrix[g1,g2] for g1 in eachindex(wave), g2 in eachindex(wave)]
    wl = build_wave_lookup(wave)
    
    HF_eigenvalue=zeros(Float64,dimension)
    HF_eigenvector=zeros(ComplexF64,dimension,dimension)
    HartreeMatrix=zeros(ComplexF64,dimension,dimension)
    FockMatrix=zeros(ComplexF64,dimension,dimension)
  

    DIIS_input_DensityMatrix=Vector{Matrix{ComplexF64}}(undef,6)
    DIIS_input_DeltaMatrix=Vector{Matrix{ComplexF64}}(undef,6)
    input_DensityMatrix=initial_DensityMatrix
    bad_count=0
    energy=0.0
    energy_change=0.0
  
   println(Threads.nthreads())
   while (eout>1*10^(-22)) || (bad_count<4) || (abs(energy_change)>1*10^(-10))
      if eout<1*10^(-22)
       bad_count+=1
      end
      
      tic=time()

      if (itcount>100 && abs(eout)>10^(-2)) || (itcount>30 && abs(eout)<10^(-6))
      
        dmk=implement_DIIS(DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix)
        if dmk==0
            itcount=0
           
              A=randn(dimension,dimension)+im*randn(dimension,dimension)
              dmk+=(A+A')*0.01
         
        end
      

       eout,energy_change,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,6)+1],HF_eigenvalue,HF_eigenvector,energy,HartreeMatrix,FockMatrix=Construct_DensityMatrix(wl,wave,wave_n1,wave_n2,
                                                                                                                                                                             dmk,single_Ham,
                                                                                                                                                                            single_MoirePo,overlapmatrix,
                                                                                                                                                                            energy,filling,Area,Coulomb_matrix)
       
        DIIS_input_DensityMatrix[mod(itcount,6)+1]=dmk
        input_DensityMatrix=output_DensityMatrix
        println("using DIIS")
       
      else
  
      

        eout,energy_change,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,6)+1],HF_eigenvalue,HF_eigenvector,energy,HartreeMatrix,FockMatrix=Construct_DensityMatrix(wl,wave,wave_n1,wave_n2,
                                                                                                                                                                            input_DensityMatrix,single_Ham,
                                                                                                                                                                            single_MoirePo,overlapmatrix,
                                                                                                                                                                            energy,filling,Area,Coulomb_matrix)
        DIIS_input_DensityMatrix[mod(itcount,6)+1]=input_DensityMatrix
        input_DensityMatrix=output_DensityMatrix
        
       
         
     

      end

    



      itcount+=1
     
      toc=time()
      println(toc-tic,"eout=$eout","energy_change=$energy_change","itcount=$itcount")
      flush(stdout)
     
    
  end







  return DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,eout,HartreeMatrix,FockMatrix

end


















function implement_DIIS(DIIS_input_projector::Vector{Matrix{ComplexF64}},DIIS_input_DeltaMatrix::Vector{Matrix{ComplexF64}})



      Bmatrix=zeros(ComplexF64,7,7)
      for ja in 1:6
       Bmatrix[ja,7]=1
       Bmatrix[7,ja]=1
      end
  
      for ja in 1:6,jb in 1:6
   
             Bmatrix[ja,jb]+=real(tr((DIIS_input_DeltaMatrix[ja])'*(DIIS_input_DeltaMatrix[jb])))
          
      end

      inB=safe_inverse(Bmatrix)
      if inB≠0
         coeff=inB*[0;0;0;0;0;0;1]
         dmk=coeff[1]*(DIIS_input_projector[1]+DIIS_input_DeltaMatrix[1])
         for ja in 2:6
           dmk+=coeff[ja]*(DIIS_input_projector[ja]+DIIS_input_DeltaMatrix[ja])
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
          return pinv(A, 0.1)  # Use pseudoinverse as an alternative
      else
          rethrow(e)  # If another error occurs, propagate it
      end
  end
end