using LinearAlgebra
using Arpack
using Combinatorics
using Random



 

function overlap(k::Vector{Float64},q::Vector{Float64},spin::Float64,M::Float64)::ComplexF64
    v=(M^2+norm(k)^2+k[1]*q[1]+k[2]*q[2]-im*(k[1]*q[2]-k[2]*q[1]))^(Int(2*spin))/((M^2+norm(k)^2)^spin*(M^2+norm(k+q)^2)^spin)
  
    return v
end
 











function triangle_initial_Densitymatrix(spin::Float64,vf::Float64,V0::Float64,scale::Float64,Nq::Int64,filling::Int,gcutoff::Float64)
    am=2*π/(scale);
    
    mass=0.5;
    
    b1=scale*[0,1]
    b2=scale*[1,0]


    a1m=am*[0,1]
    a2m=am*[1,0]
    
    
    T1=b1/(Nq)
    T2=b2/(Nq)
    
    
    
    b1T=Int.(round.(inv([T1 T2])*b1))
    b2T=Int.(round.(inv([T1 T2])*b2))
    
    
    
    allowedq=Vector{Int64}[]
    for ja in 0:Nq-1,jb in 0:Nq-1
        push!(allowedq,[ja,jb])
    end

    
    
    
    wave=Vector{Int64}[]
    cutoff=18
    cutoffstandard=gcutoff*scale
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        gtest=ja*b1+jb*b2;
        if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
            push!(wave,ja*b1T+jb*b2T)
        end
    end
    dimension=length(wave)
   
    
    
    overlapmatrix=zeros(ComplexF64,Nq^2,length(wave),Nq^2,length(wave))
    for ja in 1:Nq^2, jb in eachindex(wave), jc in 1:Nq^2, jd in eachindex(wave)
       overlapmatrix[ja,jb,jc,jd]=overlap([T1 T2]*(allowedq[ja]+wave[jb]),[T1 T2]*(allowedq[jc]+wave[jd]-wave[jb]-allowedq[ja]),spin,mass*vf/2)
    end


    
     
    single_Ham=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
    single_MoirePo=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
    single_eigenvalue=[zeros(Float64,dimension) for _ in 1:Nq^2]
    single_eigenvector=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
    
    
    Threads.@threads for ja in 1:Nq^2
      
      
      k=allowedq[ja][1]*T1+allowedq[ja][2]*T2
    
      for jb in eachindex(wave)
       single_Ham[ja][jb,jb]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2)^2/(2*mass)
      end



     for jc in eachindex(wave)
        k1=k+wave[jc][1]*T1+wave[jc][2]*T2
        pos=findfirst(item->item==wave[jc]-b1T,wave)
        if pos≠nothing
          single_MoirePo[ja][jc,pos]=V0*overlap(k1,-b1,spin,mass*vf/2)
        end
        
      
        pos=findfirst(item->item==wave[jc]-b2T,wave)
        if pos≠nothing
            single_MoirePo[ja][jc,pos]=V0*overlap(k1,-b2,spin,mass*vf/2)
        end
     end
    
     single_MoirePo[ja]=single_MoirePo[ja]+single_MoirePo[ja]'
     FFF=eigen(single_MoirePo[ja]+single_Ham[ja])
    
      single_eigenvalue[ja]=real(FFF.values)
      single_eigenvector[ja]=FFF.vectors
    
    end


    input_DensityMatrix=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
    for ja in 1:Nq^2,jb in 1:filling
     
       input_DensityMatrix[ja]+=single_eigenvector[ja][:,jb]*(single_eigenvector[ja][:,jb])'
    end    
    
     

    

    
    return overlapmatrix, wave, input_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,single_eigenvector,allowedq, T1, T2, a1m, a2m
      

       
end















function construct_loop_dic(wave::Vector{Vector{Int}})::Dict{Vector{Int},Any}
    g_dic=Dict{Vector{Int},Int}()
    for ja in eachindex(wave)
      g_dic[wave[ja]]=ja
    end

    loop_dic=Dict{Vector{Int},Any}()

    for gindex in eachindex(wave), g2index in eachindex(wave)
        
        if !haskey(loop_dic,wave[g2index]-wave[gindex])
           loop_dic[wave[g2index]-wave[gindex]]=Dict{Vector{Int},Vector{Vector{Int}}}()
        end
        
        loop_dic[wave[g2index]-wave[gindex]][[gindex,g2index]]=Vector{Int64}[]
        
         
    end
    
    for gindex in eachindex(wave), g1index in eachindex(wave),g2index in eachindex(wave)
       
        if haskey(g_dic,wave[gindex]+wave[g1index]-wave[g2index])
            push!(loop_dic[wave[g2index]-wave[gindex]][[gindex,g2index]],[g1index,g_dic[wave[gindex]+wave[g1index]-wave[g2index]]])
        end
        
    end
    
   return loop_dic

end






function Construct_DensityMatrix(loop_dic::Dict{Vector{Int},Any},allowedq::Vector{Vector{Int}},
                               T1::Vector{Float64},T2::Vector{Float64},Nq::Int64,wave::Vector{Vector{Int64}},
                               input_DensityMatrix::Vector{Matrix{ComplexF64}},single_Ham::Vector{Matrix{ComplexF64}},
                               single_MoirePo::Vector{Matrix{ComplexF64}},constq::Float64,overlapmatrix::Array{ComplexF64,4},
                               energy_input::Float64,filling::Int)
  
 
   dimension=length(wave)
  HartreeMatrix=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
  FockMatrix=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
  output_DensityMatrix=Vector{Matrix{ComplexF64}}(undef,Nq^2)
  DeltaMatrix=Vector{Matrix{ComplexF64}}(undef,Nq^2)
  NewDensityMatrix=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
  HF_eigenvalue=Vector{Vector{Float64}}(undef,Nq^2)
  HF_eigenvector=Vector{Matrix{ComplexF64}}(undef,Nq^2)

 
  Threads.@threads for jk in 1:Nq^2
    Fk = FockMatrix[jk]
    for jk1 in 1:Nq^2
        dmk = input_DensityMatrix[jk1]
      
       for (dg,loop_dic_dg) in loop_dic
          
           for (gg2,loop_dic_dg_gg2) in loop_dic_dg          
               CoulF=overlapmatrix[jk1,gg2[2],jk,gg2[1]]
           for g1g3 in loop_dic_dg_gg2
               Fk[g1g3[2],gg2[1]]+=dmk[g1g3[1],gg2[2]]*CoulF*overlapmatrix[jk,g1g3[2],jk1,g1g3[1]]
           end 
           end
   
       end    
    end
  end
 
  Hartree_Density=zeros(ComplexF64,dimension,dimension)
  for ja in 1:Nq^2
   Hartree_Density+=input_DensityMatrix[ja] .* transpose(overlapmatrix[ja,:,ja,:])
  end

  


  Threads.@threads for jk in 1:Nq^2
      for dg in keys(loop_dic)
         
          for gg2 in keys(loop_dic[dg])          
              CoulH=overlapmatrix[jk,gg2[2],jk,gg2[1]]            
          for g1g3 in loop_dic[dg][gg2]           
              HartreeMatrix[jk][gg2[2],gg2[1]]+=Hartree_Density[g1g3[1],g1g3[2]]*CoulH               
          end 
          end
  
      end    
 end


 

 for ja in 1:Nq^2
   FFF=eigen(single_MoirePo[ja]+single_Ham[ja]+constq*HartreeMatrix[ja]-constq*FockMatrix[ja])
   HF_eigenvalue[ja]=real(FFF.values)
   HF_eigenvector[ja]=FFF.vectors
   
 end

 bound=(sort(reduce(vcat,HF_eigenvalue))[filling*Nq^2+1]+sort(reduce(vcat,HF_eigenvalue))[filling*Nq^2])/2

  for ja in 1:Nq^2
    
       for jd in eachindex(HF_eigenvalue[ja])
          if HF_eigenvalue[ja][jd]<bound
             NewDensityMatrix[ja]+=HF_eigenvector[ja][:,jd]*(HF_eigenvector[ja][:,jd])'
          end
       end
       DeltaMatrix[ja]=NewDensityMatrix[ja]-input_DensityMatrix[ja]
       output_DensityMatrix[ja]=0.0*input_DensityMatrix[ja]+1.0*NewDensityMatrix[ja]
  end


  
  e1=0.0
  for ja in 1:Nq^2
    e1+=tr(DeltaMatrix[ja]'*DeltaMatrix[ja])
  end
  eout=real(e1)/Nq^2
  
  
  energy=0
   for ja in 1:Nq^2
       ss=single_MoirePo[ja]+single_Ham[ja]+0.5*constq*HartreeMatrix[ja]-0.5*constq*FockMatrix[ja]
       energy+=real(tr(ss*output_DensityMatrix[ja]))
   end

   energy_change=real(energy-energy_input)



 return  eout,energy_change,output_DensityMatrix,DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,HartreeMatrix,FockMatrix
end




function iteration_loop(initial_DensityMatrix::Vector{Matrix{ComplexF64}},
                       allowedq::Vector{Vector{Int64}},T1::Vector{Float64},T2::Vector{Float64},
                       Nq::Int,wave::Vector{Vector{Int64}},single_Ham::Vector{Matrix{ComplexF64}},
                       single_MoirePo::Vector{Matrix{ComplexF64}},constq::Float64,overlapmatrix::Array{ComplexF64,4},filling::Int)
    eout=1.0
    itcount=0
    dimension=length(wave)
    loop_dic=construct_loop_dic(wave)
    HF_eigenvalue=Vector{Vector{Float64}}(undef,Nq^2)
    HF_eigenvector=Vector{Matrix{ComplexF64}}(undef,Nq^2)
    HartreeMatrix=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
    FockMatrix=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
  

    DIIS_input_DensityMatrix=Vector{Vector{Matrix{ComplexF64}}}(undef,3)
    DIIS_input_DeltaMatrix=Vector{Vector{Matrix{ComplexF64}}}(undef,3)
    input_DensityMatrix=initial_DensityMatrix
    bad_count=0
    energy=0.0
    energy_change=0.0
  

    while (eout>1*10^(-16)) || (bad_count<4) || (energy_change>1*10^(-8))
      if eout<1*10^(-16)
       bad_count+=1
      end
      
      tic=time()

      if (itcount>200 && abs(eout)>10^(-2)) || (itcount>30 && abs(eout)<10^(-8))
      
        dmk=implement_DIIS(DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,Nq)
        if dmk==0
            itcount=0
            dmk=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
            for ja in 1:Nq^2
              A=randn(dimension,dimension)+im*randn(dimension,dimension)
              dmk[ja]+=(A+A')*0.01
            end
        end
      

       eout,energy_change,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue,HF_eigenvector,energy,HartreeMatrix,FockMatrix=Construct_DensityMatrix(loop_dic,allowedq,T1,T2,Nq,wave,dmk,single_Ham,single_MoirePo,constq,overlapmatrix,energy,filling)
       
        DIIS_input_DensityMatrix[mod(itcount,3)+1]=dmk
        input_DensityMatrix=output_DensityMatrix
        println("using DIIS")
       
      else
  
      

        eout,energy_change,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue,HF_eigenvector,energy,HartreeMatrix,FockMatrix=Construct_DensityMatrix(loop_dic,allowedq,T1,T2,Nq,wave,input_DensityMatrix,single_Ham,single_MoirePo,constq,overlapmatrix,energy,filling)
        DIIS_input_DensityMatrix[mod(itcount,3)+1]=input_DensityMatrix
        input_DensityMatrix=output_DensityMatrix
        
       
         
     

      end

    



      itcount+=1
     
      toc=time()
      println(toc-tic,"eout=$eout","energy_change=$energy_change","itcount=$itcount")
      flush(stdout)
     
    
  end







  return DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,eout,HartreeMatrix,FockMatrix

end








function metric(wavelist::Vector{Vector{Int64}},spin::Float64,vf::Float64,k::Vector{Float64},q::Vector{Float64},T1::Vector{Float64},T2::Vector{Float64})::Matrix{ComplexF64}
    Amatrix=zeros(ComplexF64,length(wavelist),length(wavelist))
    mass=0.5
    for ja in 1:length(wavelist)
    Amatrix[ja,ja]=overlap(k+wavelist[ja][1]*T1+wavelist[ja][2]*T2,q,spin,mass*vf/2)
    end
    return Amatrix
end






function shift_vector(vector_toshift::Vector{ComplexF64},shiftamount::Vector{Int},wave::Vector{Vector{Int}})
    shifted_vector=zeros(ComplexF64,length(wave))
    for ja in eachindex(wave)
        pos=findfirst(item->item==wave[ja]+shiftamount,wave)
        if pos≠nothing
        shifted_vector[ja]=vector_toshift[pos]
        end
    
    end
    return shifted_vector/norm(shifted_vector)

end



function triangle_chern(Nq::Int,wave::Vector{Vector{Int}},scale::Float64,spin::Float64,vf::Float64,allowedq::Vector{Vector{Int}},HF_eigenvector::Vector{Matrix{ComplexF64}},single_eigenvector::Vector{Matrix{ComplexF64}})

  
    
    mass=0.5;
    dimension=length(wave)
    
    b1=scale*[0,1]
    b2=scale*[1,0]


   
    
    
    T1=b1/(Nq)
    T2=b2/(Nq)
    b1T=Int.(round.(inv([T1 T2])*b1))
    b2T=Int.(round.(inv([T1 T2])*b2))
    
    
    chern_allowedq=Vector{Int64}[]
    for ja in 0:Nq,jb in 0:Nq
        push!(chern_allowedq,[ja,jb])
    end

  

    eigenvector_intermediate_bc=Vector{Vector{ComplexF64}}(undef,(Nq+1)^2)
    eigenvector_intermediate_single=Vector{Vector{ComplexF64}}(undef,(Nq+1)^2)

    Threads.@threads for ja in eachindex(chern_allowedq)
      
      kbraket=mod.(chern_allowedq[ja],Nq)
      kbraket_pos=findfirst(item->item==kbraket,allowedq)
      if kbraket==chern_allowedq[ja]
        eigenvector_intermediate_bc[ja]=HF_eigenvector[kbraket_pos][:,1] 
        eigenvector_intermediate_single[ja]=single_eigenvector[kbraket_pos][:,1] 
      else
        eigenvector_intermediate_bc[ja]=shift_vector(HF_eigenvector[kbraket_pos][:,1],chern_allowedq[ja]-kbraket,wave)
        eigenvector_intermediate_single[ja]=shift_vector(single_eigenvector[kbraket_pos][:,1],chern_allowedq[ja]-kbraket,wave)

      end  
    end

    eigenvector_bc=zeros(ComplexF64,dimension,Nq+1,Nq+1)
    eigenvector_bc_single=zeros(ComplexF64,dimension,Nq+1,Nq+1)
    for ja in eachindex(chern_allowedq)
       eigenvector_bc[:,chern_allowedq[ja][1]+1,chern_allowedq[ja][2]+1]=eigenvector_intermediate_bc[ja]
       eigenvector_bc_single[:,chern_allowedq[ja][1]+1,chern_allowedq[ja][2]+1]=eigenvector_intermediate_single[ja]
    end
    
    
    Uonelink=zeros(ComplexF64,Nq,Nq+1)
    Utwolink=zeros(ComplexF64,Nq+1,Nq)
    tra=zeros(ComplexF64,Nq,Nq)
    tra_single=zeros(ComplexF64,Nq,Nq)
    
    for ja in 1:Nq, jb in 1:Nq+1
        Amatrix=metric(wave,spin,vf,(ja-1)*T1+(jb-1)*T2,T1,T1,T2)
       Uonelink[ja,jb]=dot(eigenvector_bc[:,ja,jb],Amatrix*eigenvector_bc[:,ja+1,jb])/abs(dot(eigenvector_bc[:,ja,jb],Amatrix*eigenvector_bc[:,ja+1,jb]))
    end

    
    
    for ja in 1:Nq+1, jb in 1:Nq
        Amatrix=metric(wave,spin,vf,(ja-1)*T1+(jb-1)*T2,T2,T1,T2)
     Utwolink[ja,jb]=dot(eigenvector_bc[:,ja,jb],Amatrix*eigenvector_bc[:,ja,jb+1])/abs(dot(eigenvector_bc[:,ja,jb],Amatrix*eigenvector_bc[:,ja,jb+1]))
    end
    
    dG=norm(T2)
    for ja in 1:Nq, jb in 1:Nq
        Amatrix=metric(wave,spin,vf,(ja-1)*T1+(jb-1)*T2,T1,T1,T2)
        Bmatrix=metric(wave,spin,vf,(ja-1)*T1+(jb-1)*T2,T2,T1,T2)
        Cmatrix=metric(wave,spin,vf,(ja-1)*T1+(jb-1)*T2,T2+T1,T1,T2)
        A1=dot(eigenvector_bc[:,ja,jb],Amatrix*eigenvector_bc[:,ja+1,jb])
        B1=dot(eigenvector_bc[:,ja,jb],Bmatrix*eigenvector_bc[:,ja,jb+1])
    
       gyy=(1-abs(A1)^2)/dG^2
       gxx=(1-abs(B1)^2)/dG^2
       tra[ja,jb]+=(gxx+gyy)*dG^2

        A1=dot(eigenvector_bc_single[:,ja,jb],Amatrix*eigenvector_bc_single[:,ja+1,jb])
        B1=dot(eigenvector_bc_single[:,ja,jb],Bmatrix*eigenvector_bc_single[:,ja,jb+1])
       
       gyy=(1-abs(A1)^2)/dG^2
       gxx=(1-abs(B1)^2)/dG^2
       tra_single[ja,jb]+=(gxx+gyy)*dG^2
      
      
    end

   


    Flink=zeros(ComplexF64,Nq,Nq)
    for ja in 1:Nq, jb in 1:Nq
     Flink[ja,jb]=log(Uonelink[ja,jb]*Utwolink[ja+1,jb]/(Uonelink[ja,jb+1]*Utwolink[ja,jb]))
    end
    chern=sum(Flink)/(2*π*im)
    aveF=sum(Flink)/Nq^2
    uniform=0
    for ja in 1:Nq, jb in 1:Nq
        uniform+=(imag(Flink[ja,jb])-imag(aveF))^2*Nq^2/(2π)^2
    end


    Uonelink=zeros(ComplexF64,Nq,Nq+1)
    Utwolink=zeros(ComplexF64,Nq+1,Nq)
  
    
    for ja in 1:Nq, jb in 1:Nq+1
        Amatrix=metric(wave,spin,vf,(ja-1)*T1+(jb-1)*T2,T1,T1,T2)
       Uonelink[ja,jb]=dot(eigenvector_bc_single[:,ja,jb],Amatrix*eigenvector_bc_single[:,ja+1,jb])/abs(dot(eigenvector_bc_single[:,ja,jb],Amatrix*eigenvector_bc_single[:,ja+1,jb]))
    end
    
    for ja in 1:Nq+1, jb in 1:Nq
        Amatrix=metric(wave,spin,vf,(ja-1)*T1+(jb-1)*T2,T2,T1,T2)
     Utwolink[ja,jb]=dot(eigenvector_bc_single[:,ja,jb],Amatrix*eigenvector_bc_single[:,ja,jb+1])/abs(dot(eigenvector_bc_single[:,ja,jb],Amatrix*eigenvector_bc_single[:,ja,jb+1]))
    end
    

   

    Flink_single=zeros(ComplexF64,Nq,Nq)
    for ja in 1:Nq, jb in 1:Nq
     Flink_single[ja,jb]=log(Uonelink[ja,jb]*Utwolink[ja+1,jb]/(Uonelink[ja,jb+1]*Utwolink[ja,jb]))
    end
    chern_single=sum(Flink_single)/(2*π*im)
    aveF=sum(Flink_single)/Nq^2
    uniform_single=0
    for ja in 1:Nq, jb in 1:Nq
        uniform_single+=(imag(Flink_single[ja,jb])-imag(aveF))^2*Nq^2/(2π)^2
    end
    
    

    
  trace_condition=sum(tra)-sum(abs.(Flink))
  trace_condition_single=sum(tra_single)-sum(abs.(Flink_single))
      


    
      
   return chern,Flink,chern_single,Flink_single,trace_condition,trace_condition_single,uniform,uniform_single

end



function implement_DIIS(DIIS_input_projector::Vector{Vector{Matrix{ComplexF64}}},DIIS_input_DeltaMatrix::Vector{Vector{Matrix{ComplexF64}}},Nq::Int)



      Bmatrix=zeros(ComplexF64,4,4)
      for ja in 1:3
       Bmatrix[ja,4]=1
       Bmatrix[4,ja]=1
      end
  
      for ja in 1:3,jb in 1:3
          for jc in Nq^2
             Bmatrix[ja,jb]+=real(tr((DIIS_input_DeltaMatrix[ja][jc])'*(DIIS_input_DeltaMatrix[jb][jc])))
          end
      end

      inB=safe_inverse(Bmatrix)
      if inB≠0
         coeff=inB*[0;0;0;1]
         dmk=coeff[1]*(DIIS_input_projector[1]+DIIS_input_DeltaMatrix[1])+coeff[2]*(DIIS_input_projector[2]+DIIS_input_DeltaMatrix[2])+coeff[3]*(DIIS_input_projector[3]+DIIS_input_DeltaMatrix[3])
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
          return 0  # Use pseudoinverse as an alternative
      else
          rethrow(e)  # If another error occurs, propagate it
      end
  end
end