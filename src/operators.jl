using LinearAlgebra
using Arpack
using Combinatorics
using LinearAlgebra
using Random



function overlap(k::Vector{Float64},q::Vector{Float64},β::Float64)::ComplexF64
    v=q[1]^2+q[2]^2+2*im*(k[1]*q[2]-k[2]*q[1])
    #v=2*im*(k[1]*q[2]-k[2]*q[1])
    return exp(-β/4*v)
end
 
 
 
 
function Coulomb(k::Vector{Int64},T1::Vector{Float64},T2::Vector{Float64})::Float64
  
 
   return k==[0,0] ? 0.0 : 1/norm(k[1]*T1+k[2]*T2)
   
end


function square_initial_Densitymatrix(flux::Float64,V0::Float64,ϕ::Float64,scale::Float64,Nq::Int64)
    am=2*π/(scale);
    β=flux/(scale^2)
    
    mass=0.5;
    
    b1=scale*[1,0]
    b2=scale*[0,1]
    
    a1m=am*[1,0]
    a2m=am*[0,1]
    
    
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
    cutoffstandard=4.01*scale
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        gtest=ja*b1+jb*b2;
        if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
            push!(wave,ja*b1T+jb*b2T)
        end
    end
    dimension=length(wave)
   
    
    
    overlapmatrix=zeros(ComplexF64,Nq^2,length(wave),Nq^2,length(wave))
    for ja in 1:Nq^2, jb in eachindex(wave), jc in 1:Nq^2, jd in eachindex(wave)
       overlapmatrix[ja,jb,jc,jd]=overlap([T1 T2]*(allowedq[ja]+wave[jb]),[T1 T2]*(allowedq[jc]+wave[jd]-wave[jb]-allowedq[ja]),β)
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
          single_MoirePo[ja][jc,pos]=V0*exp(im*ϕ)*overlap(k1,-b1,β)
        end
        
    
        pos=findfirst(item->item==wave[jc]-b2T,wave)
        if pos≠nothing
            single_MoirePo[ja][jc,pos]=V0*exp(im*ϕ)*overlap(k1,-b2,β)
        end
     end
    
     single_MoirePo[ja]=single_MoirePo[ja]+single_MoirePo[ja]'
     FFF=eigen(single_MoirePo[ja]+single_Ham[ja])
    
      single_eigenvalue[ja]=real(FFF.values)
      single_eigenvector[ja]=FFF.vectors
    
    end


    input_DensityMatrix=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
    for ja in 1:Nq^2
     
       input_DensityMatrix[ja]+=(single_eigenvector[ja][:,1]*(single_eigenvector[ja][:,1])')
        
    end
     

    

    
    return overlapmatrix, wave, input_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,allowedq, T1, T2, a1m, a2m
      

       
end






function triangle_initial_Densitymatrix(flux::Float64,V0::Float64,ϕ::Float64,scale::Float64,Nq::Int64)
    am=4*π/(√3*scale);
    β=4*flux/(√3*scale^2)
    mass=0.5;
    
    b1=4*π/(√3*am)*[0,1]
    b2=4*π/(√3*am)*[√3/2,-1/2]


    a1m=am*[1/2,√3/2]
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
    cutoffstandard=4.01*scale
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        gtest=ja*b1+jb*b2;
        if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
            push!(wave,ja*b1T+jb*b2T)
        end
    end
    dimension=length(wave)
   
    
    
    overlapmatrix=zeros(ComplexF64,Nq^2,length(wave),Nq^2,length(wave))
    for ja in 1:Nq^2, jb in eachindex(wave), jc in 1:Nq^2, jd in eachindex(wave)
       overlapmatrix[ja,jb,jc,jd]=overlap([T1 T2]*(allowedq[ja]+wave[jb]),[T1 T2]*(allowedq[jc]+wave[jd]-wave[jb]-allowedq[ja]),β)
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
          single_MoirePo[ja][jc,pos]=V0*exp(im*ϕ)*overlap(k1,-b1,β)
        end
        
      
        pos=findfirst(item->item==wave[jc]+b2T+b1T,wave)
        if pos≠nothing
            single_MoirePo[ja][jc,pos]=V0*exp(im*ϕ)*overlap(k1,+b2+b1,β)
        end
    
        pos=findfirst(item->item==wave[jc]-b2T,wave)
        if pos≠nothing
            single_MoirePo[ja][jc,pos]=V0*exp(im*ϕ)*overlap(k1,-b2,β)
        end
     end
    
     single_MoirePo[ja]=single_MoirePo[ja]+single_MoirePo[ja]'
     FFF=eigen(single_MoirePo[ja]+single_Ham[ja])
    
      single_eigenvalue[ja]=real(FFF.values)
      single_eigenvector[ja]=FFF.vectors
    
    end


    input_DensityMatrix=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
    for ja in 1:Nq^2
     
       input_DensityMatrix[ja]+=(single_eigenvector[ja][:,1]*(single_eigenvector[ja][:,1])')
        
    end
     

    

    
    return overlapmatrix, wave, input_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,allowedq, T1, T2, a1m, a2m
      

       
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






function Construct_DensityMatrix(loop_dic::Dict{Vector{Int},Any},allowedq::Vector{Vector{Int}},T1::Vector{Float64},T2::Vector{Float64},Nq::Int64,wave::Vector{Vector{Int64}},input_DensityMatrix::Vector{Matrix{ComplexF64}},single_Ham::Vector{Matrix{ComplexF64}},single_MoirePo::Vector{Matrix{ComplexF64}},constq::Float64,overlapmatrix::Array{ComplexF64,4})::Tuple{Float64,Vector{Matrix{ComplexF64}},Vector{Matrix{ComplexF64}},Vector{Vector{Float64}}}
  
 
   dimension=length(wave)
  HartreeMatrix=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
  FockMatrix=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
  output_DensityMatrix=Vector{Matrix{ComplexF64}}(undef,Nq^2)
  DeltaMatrix=Vector{Matrix{ComplexF64}}(undef,Nq^2)
  NewDensityMatrix=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
  HF_eigenvalue=Vector{Vector{Float64}}(undef,Nq^2)
  HF_eigenvector=Vector{Any}(undef,Nq^2)

 
  Threads.@threads for jk in 1:Nq^2
    Fk = FockMatrix[jk]
    for jk1 in 1:Nq^2
        dmk = input_DensityMatrix[jk1]
        q=allowedq[jk1]-allowedq[jk]
       for (dg,loop_dic_dg) in loop_dic
           CoulF1=Coulomb(q+dg,T1,T2)     
           for (gg2,loop_dic_dg_gg2) in loop_dic_dg          
               CoulF=CoulF1*overlapmatrix[jk1,gg2[2],jk,gg2[1]]
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
          CoulH1=Coulomb(dg,T1,T2)
          for gg2 in keys(loop_dic[dg])          
              CoulH=CoulH1*overlapmatrix[jk,gg2[2],jk,gg2[1]]            
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

 bound=sort(reduce(vcat,HF_eigenvalue))[Nq^2+1]

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
  eout=real(e1)
  
  
 

 return  eout,output_DensityMatrix,DeltaMatrix,HF_eigenvalue
end




function iteration_loop(initial_DensityMatrix::Vector{Matrix{ComplexF64}},allowedq::Vector{Vector{Int64}},T1::Vector{Float64},T2::Vector{Float64},Nq::Int,wave::Vector{Vector{Int64}},single_Ham::Vector{Matrix{ComplexF64}},single_MoirePo::Vector{Matrix{ComplexF64}},constq::Float64,overlapmatrix::Array{ComplexF64,4})::Tuple{Vector{Vector{Matrix{ComplexF64}}},Vector{Vector{Matrix{ComplexF64}}},Vector{Vector{Float64}}}
    eout=1.0
    itcount=0
    bad_count=0
    loop_dic=construct_loop_dic(wave)
    HF_eigenvalue=Vector{Any}(undef,Nq^2)
    DIIS_input_DensityMatrix=Vector{Vector{Matrix{ComplexF64}}}(undef,3)
    DIIS_input_DeltaMatrix=Vector{Vector{Matrix{ComplexF64}}}(undef,3)
    input_DensityMatrix=initial_DensityMatrix
    while (eout>1*10^-13) || (bad_count<4)
      if eout<1*10^-13 
      bad_count+=1
      end
      tic=time()
      eout,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue=Construct_DensityMatrix(loop_dic,allowedq,T1,T2,Nq,wave,input_DensityMatrix,single_Ham,single_MoirePo,constq,overlapmatrix)
      DIIS_input_DensityMatrix[mod(itcount,3)+1]=input_DensityMatrix
      input_DensityMatrix=output_DensityMatrix
      itcount+=1
      toc=time()
      println(toc-tic,"eout=$eout")
      flush(stdout)
     
    
    end
    
    println("startDIIS",itcount)
    
    bad_count=0
    while (eout>10^-13) || (bad_count<4)
        if eout<1*10^-13 
            bad_count+=1
        end
        tic=time()
        Bmatrix=zeros(ComplexF64,4,4)
        for ja in 1:3
         Bmatrix[ja,4]=1
         Bmatrix[4,ja]=1
        end
    
        for ja in 1:3,jb in 1:3
            for jc in 1:Nq^2
               Bmatrix[ja,jb]+=tr((DIIS_input_DeltaMatrix[ja][jc])'*(DIIS_input_DeltaMatrix[jb][jc]))
            end
        end
        coeff=inv(Bmatrix)*[0;0;0;1]
       
        dmk=coeff[1]*(DIIS_input_DensityMatrix[1]+DIIS_input_DeltaMatrix[1])+coeff[2]*(DIIS_input_DensityMatrix[2]+DIIS_input_DeltaMatrix[2])+coeff[3]*(DIIS_input_DensityMatrix[3]+DIIS_input_DeltaMatrix[3])
         eout,_,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue=Construct_DensityMatrix(loop_dic,allowedq,T1,T2,Nq,wave,dmk,single_Ham,single_MoirePo,constq,overlapmatrix)
        DIIS_input_DensityMatrix[mod(itcount,3)+1]=dmk
        itcount+=1

        toc=time()
        println(toc-tic,"eout=$eout")
        flush(stdout)
    end

 


  return DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue

end


function Densitymap(a1m::Vector{Float64},a2m::Vector{Float64},overlapmatrix::Array{ComplexF64,4},wave::Vector{Vector{Int}},input_DensityMatrix::Vector{Matrix{ComplexF64}})
    N3=50
    dimension=length(wave)
    zgrid=zeros(Float64,N3,N3)
    Hartree_Density=zeros(ComplexF64,dimension,dimension)
    for jk1 in eachindex(input_DensityMatrix)
       Hartree_Density+=input_DensityMatrix[jk1] .* transpose(overlapmatrix[jk1,:,jk1,:])
    end
    
    for ja in 1:50, jb in 1:50
      rvec=ja/25*a1m+jb/25*a2m
      for jc in eachindex(wave), jd in eachindex(wave)
        gvec=[T1 T2]*(wave[jc]-wave[jd])
       zgrid[ja,jb]+=real(Hartree_Density[jc,jd]*exp(im*(gvec[1]*rvec[1]+gvec[2]*rvec[2])))
      end
    
    end
    
    return zgrid
end





function metric(wavelist::Vector{Vector{Int64}},β::Float64,k::Vector{Float64},q::Vector{Float64},T1::Vector{Float64},T2::Vector{Float64})::Matrix{ComplexF64}
    Amatrix=zeros(ComplexF64,length(wavelist),length(wavelist))
    for ja in 1:length(wavelist)
    Amatrix[ja,ja]=overlap(k+wavelist[ja][1]*T1+wavelist[ja][2]*T2,q,β)
    end
    return Amatrix
end



function Construct_HFmatrix(loop_dic::Dict{Vector{Int},Any},pathpointindex::Int64,pathpoint::Vector{Int64},allowedq::Vector{Vector{Int}},T1::Vector{Float64},T2::Vector{Float64},Nq::Int64,wave::Vector{Vector{Int64}},input_DensityMatrix::Vector{Matrix{ComplexF64}},constq::Float64,chern_overlapmatrix::Array{ComplexF64,4})::Matrix{ComplexF64}
  
  
  dimension=length(wave)
  HartreeMatrix=zeros(ComplexF64,dimension,dimension) 
  FockMatrix=zeros(ComplexF64,dimension,dimension) 
 
  


    
    for jk1 in 1:Nq^2
        dmk = input_DensityMatrix[jk1]
        q=allowedq[jk1]-pathpoint
       for (dg,loop_dic_dg) in loop_dic
           CoulF1=Coulomb(q+dg,T1,T2)     
           for (gg2,loop_dic_dg_gg2) in loop_dic_dg          
               CoulF=CoulF1*chern_overlapmatrix[jk1,gg2[2],pathpointindex,gg2[1]]
           for g1g3 in loop_dic_dg_gg2
               FockMatrix[g1g3[2],gg2[1]]+=dmk[g1g3[1],gg2[2]]*CoulF*chern_overlapmatrix[pathpointindex,g1g3[2],jk1,g1g3[1]]
           end 
           end
   
       end    
    end

  

  Hartree_Density=zeros(ComplexF64,dimension,dimension)
  for jk1 in 1:Nq^2
   Hartree_Density+=input_DensityMatrix[jk1] .* transpose(chern_overlapmatrix[jk1,:,jk1,:])
  end

 

    
 
    for dg in keys(loop_dic)
        CoulH1=Coulomb(dg,T1,T2)
        for gg2 in keys(loop_dic[dg])          
            CoulH=CoulH1*chern_overlapmatrix[pathpointindex,gg2[2],pathpointindex,gg2[1]]            
        for g1g3 in loop_dic[dg][gg2]           
            HartreeMatrix[gg2[2],gg2[1]]+=Hartree_Density[g1g3[1],g1g3[2]]*CoulH               
        end 
        end

    end    
 


 return  constq*(HartreeMatrix-FockMatrix)
end


function square_chern(Nq::Int,wave::Vector{Vector{Int}},scale::Float64,ϕ::Float64,flux::Float64,input_DensityMatrix::Vector{Matrix{ComplexF64}},constq::Float64)

  
    β=flux/(scale^2)
    mass=0.5;
    dimension=length(wave)
    
    b1=scale*[1,0]
    b2=scale*[0,1]
    
    
    T1=b1/(Nq)
    T2=b2/(Nq)
    b1T=Int.(round.(inv([T1 T2])*b1))
    b2T=Int.(round.(inv([T1 T2])*b2))
    
    
    chern_allowedq=Vector{Int64}[]
    for ja in 0:Nq,jb in 0:Nq
        push!(chern_allowedq,[ja,jb])
    end
    
    chern_overlapmatrix=zeros(ComplexF64,(Nq+1)^2,length(wave),(Nq+1)^2,length(wave))
    for ja in 1:(Nq+1)^2, jb in eachindex(wave), jc in 1:(Nq+1)^2, jd in eachindex(wave)
       chern_overlapmatrix[ja,jb,jc,jd]=overlap([T1 T2]*(chern_allowedq[ja]+wave[jb]),[T1 T2]*(chern_allowedq[jc]+wave[jd]-wave[jb]-chern_allowedq[ja]),β)
    end
    loop_dic=construct_loop_dic(wave)
    
    
    eigenvector_bc=zeros(ComplexF64,dimension,Nq+1,Nq+1)
    eigenvector_bc_single=zeros(ComplexF64,dimension,Nq+1,Nq+1)
    eigenvector_intermediate_bc=Vector{Vector{ComplexF64}}(undef,(Nq+1)^2)
    eigenvector_intermediate_single=Vector{Vector{ComplexF64}}(undef,(Nq+1)^2)


    Threads.@threads for ja in eachindex(chern_allowedq)
        k=[T1 T2]*chern_allowedq[ja]
        chern_Ham=zeros(ComplexF64,dimension,dimension)
        chern_MoirePo=zeros(ComplexF64,dimension,dimension)
        for jb in eachindex(wave)
            chern_Ham[jb,jb]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2)^2/(2*mass)
        end
    
        for jc in eachindex(wave)
            k1=k+wave[jc][1]*T1+wave[jc][2]*T2
            pos=findfirst(item->item==wave[jc]-b1T,wave)
            if pos≠nothing
           chern_MoirePo[jc,pos]=V0*exp(im*ϕ)*overlap(k1,-b1,β)
            end
            
        
            pos=findfirst(item->item==wave[jc]-b2T,wave)
            if pos≠nothing
                chern_MoirePo[jc,pos]=V0*exp(im*ϕ)*overlap(k1,-b2,β)
            end
         end
        
      chern_MoirePo=chern_MoirePo+chern_MoirePo'
      HFmatrix=Construct_HFmatrix(loop_dic,ja,chern_allowedq[ja],allowedq,T1,T2,Nq,wave,input_DensityMatrix,constq,chern_overlapmatrix)
      eigenvector_intermediate_bc[ja]=eigvecs(chern_Ham+chern_MoirePo+HFmatrix)[:,1] 
      eigenvector_intermediate_single[ja]=eigvecs(chern_Ham+chern_MoirePo)[:,1] 
    
    end

    for ja in eachindex(chern_allowedq)
        eigenvector_bc[:,chern_allowedq[ja][1]+1,chern_allowedq[ja][2]+1]=eigenvector_intermediate_bc[ja]
        eigenvector_bc_single[:,chern_allowedq[ja][1]+1,chern_allowedq[ja][2]+1]=eigenvector_intermediate_single[ja]
     end
    
    Uonelink=zeros(ComplexF64,Nq,Nq+1)
    Utwolink=zeros(ComplexF64,Nq+1,Nq)
    
    for ja in 1:Nq, jb in 1:Nq+1
        Amatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T1,T1,T2)
       Uonelink[ja,jb]=dot(eigenvector_bc[:,ja,jb],Amatrix*eigenvector_bc[:,ja+1,jb])/abs(dot(eigenvector_bc[:,ja,jb],Amatrix*eigenvector_bc[:,ja+1,jb]))
    end
    
    for ja in 1:Nq+1, jb in 1:Nq
        Amatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T2,T1,T2)
     Utwolink[ja,jb]=dot(eigenvector_bc[:,ja,jb],Amatrix*eigenvector_bc[:,ja,jb+1])/abs(dot(eigenvector_bc[:,ja,jb],Amatrix*eigenvector_bc[:,ja,jb+1]))
    end
    
    Flink=zeros(ComplexF64,Nq,Nq)
    for ja in 1:Nq, jb in 1:Nq
     Flink[ja,jb]=log(Uonelink[ja,jb]*Utwolink[ja+1,jb]/(Uonelink[ja,jb+1]*Utwolink[ja,jb]))
    end
    chern=sum(Flink)/(2*π*im)



    Uonelink=zeros(ComplexF64,Nq,Nq+1)
    Utwolink=zeros(ComplexF64,Nq+1,Nq)
    
    for ja in 1:Nq, jb in 1:Nq+1
        Amatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T1,T1,T2)
       Uonelink[ja,jb]=dot(eigenvector_bc_single[:,ja,jb],Amatrix*eigenvector_bc_single[:,ja+1,jb])/abs(dot(eigenvector_bc_single[:,ja,jb],Amatrix*eigenvector_bc_single[:,ja+1,jb]))
    end
    
    for ja in 1:Nq+1, jb in 1:Nq
        Amatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T2,T1,T2)
     Utwolink[ja,jb]=dot(eigenvector_bc_single[:,ja,jb],Amatrix*eigenvector_bc_single[:,ja,jb+1])/abs(dot(eigenvector_bc_single[:,ja,jb],Amatrix*eigenvector_bc_single[:,ja,jb+1]))
    end
    
    Flink_single=zeros(ComplexF64,Nq,Nq)
    for ja in 1:Nq, jb in 1:Nq
     Flink_single[ja,jb]=log(Uonelink[ja,jb]*Utwolink[ja+1,jb]/(Uonelink[ja,jb+1]*Utwolink[ja,jb]))
    end
    chern_single=sum(Flink_single)/(2*π*im)
    

    
    
      
   return chern,Flink,chern_single,Flink_single

end





function triangle_chern(Nq::Int,wave::Vector{Vector{Int}},scale::Float64,ϕ::Float64,flux::Float64,input_DensityMatrix::Vector{Matrix{ComplexF64}},constq::Float64)

  
    β=4*flux/(√3*scale^2)
    mass=0.5;
    dimension=length(wave)
    
    b1=scale*[0,1]
    b2=scale*[√3/2,-1/2]

    
    
    T1=b1/(Nq)
    T2=b2/(Nq)
    b1T=Int.(round.(inv([T1 T2])*b1))
    b2T=Int.(round.(inv([T1 T2])*b2))
    
    
    chern_allowedq=Vector{Int64}[]
    for ja in 0:Nq,jb in 0:Nq
        push!(chern_allowedq,[ja,jb])
    end
    
    chern_overlapmatrix=zeros(ComplexF64,(Nq+1)^2,length(wave),(Nq+1)^2,length(wave))
    for ja in 1:(Nq+1)^2, jb in eachindex(wave), jc in 1:(Nq+1)^2, jd in eachindex(wave)
       chern_overlapmatrix[ja,jb,jc,jd]=overlap([T1 T2]*(chern_allowedq[ja]+wave[jb]),[T1 T2]*(chern_allowedq[jc]+wave[jd]-wave[jb]-chern_allowedq[ja]),β)
    end
    loop_dic=construct_loop_dic(wave)
    
    
    eigenvector_intermediate_bc=Vector{Vector{ComplexF64}}(undef,(Nq+1)^2)
    eigenvector_intermediate_single=Vector{Vector{ComplexF64}}(undef,(Nq+1)^2)

    Threads.@threads for ja in eachindex(chern_allowedq)
        k=[T1 T2]*chern_allowedq[ja]
        chern_Ham=zeros(ComplexF64,dimension,dimension)
        chern_MoirePo=zeros(ComplexF64,dimension,dimension)
        for jb in eachindex(wave)
            chern_Ham[jb,jb]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2)^2/(2*mass)
        end
    
        for jc in eachindex(wave)
            k1=k+wave[jc][1]*T1+wave[jc][2]*T2
            pos=findfirst(item->item==wave[jc]-b1T,wave)
            if pos≠nothing
           chern_MoirePo[jc,pos]=V0*exp(im*ϕ)*overlap(k1,-b1,β)
            end
            
          
            pos=findfirst(item->item==wave[jc]+b2T+b1T,wave)
            if pos≠nothing
                chern_MoirePo[jc,pos]=V0*exp(im*ϕ)*overlap(k1,b2+b1,β)
            end
        
            pos=findfirst(item->item==wave[jc]-b2T,wave)
            if pos≠nothing
                chern_MoirePo[jc,pos]=V0*exp(im*ϕ)*overlap(k1,-b2,β)
            end
         end


      chern_MoirePo=chern_MoirePo+chern_MoirePo'
      HFmatrix=Construct_HFmatrix(loop_dic,ja,chern_allowedq[ja],allowedq,T1,T2,Nq,wave,input_DensityMatrix,constq,chern_overlapmatrix)
      eigenvector_intermediate_bc[ja]=eigvecs(chern_Ham+chern_MoirePo+HFmatrix)[:,1] 
      eigenvector_intermediate_single[ja]=eigvecs(chern_Ham+chern_MoirePo)[:,1] 
    
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
        Amatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T1,T1,T2)
       Uonelink[ja,jb]=dot(eigenvector_bc[:,ja,jb],Amatrix*eigenvector_bc[:,ja+1,jb])/abs(dot(eigenvector_bc[:,ja,jb],Amatrix*eigenvector_bc[:,ja+1,jb]))
    end

    
    
    for ja in 1:Nq+1, jb in 1:Nq
        Amatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T2,T1,T2)
     Utwolink[ja,jb]=dot(eigenvector_bc[:,ja,jb],Amatrix*eigenvector_bc[:,ja,jb+1])/abs(dot(eigenvector_bc[:,ja,jb],Amatrix*eigenvector_bc[:,ja,jb+1]))
    end
    
    dG=norm(T2)
    for ja in 1:Nq, jb in 1:Nq
        Amatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T1,T1,T2)
        Bmatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T2,T1,T2)
        Cmatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T2+T1,T1,T2)
        A1=dot(eigenvector_bc[:,ja,jb],Amatrix*eigenvector_bc[:,ja+1,jb])
        B1=dot(eigenvector_bc[:,ja,jb],Bmatrix*eigenvector_bc[:,ja,jb+1])
        C1=dot(eigenvector_bc[:,ja,jb],Cmatrix*eigenvector_bc[:,ja+1,jb+1])
       gyy=(1-abs(A1)^2)/dG^2
       gxx=(2-1/2*gyy*dG^2-abs(B1)^2-abs(C1)^2)/(1.5*dG^2)
       tra[ja,jb]+=(gxx+gyy)*3^(1/2)/2*dG^2

        A1=dot(eigenvector_bc_single[:,ja,jb],Amatrix*eigenvector_bc_single[:,ja+1,jb])
        B1=dot(eigenvector_bc_single[:,ja,jb],Bmatrix*eigenvector_bc_single[:,ja,jb+1])
        C1=dot(eigenvector_bc_single[:,ja,jb],Cmatrix*eigenvector_bc_single[:,ja+1,jb+1])
       gyy=(1-abs(A1)^2)/dG^2
       gxx=(2-1/2*gyy*dG^2-abs(B1)^2-abs(C1)^2)/(1.5*dG^2)
       tra_single[ja,jb]+=(gxx+gyy)*3^(1/2)/2*dG^2
      
      
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
        Amatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T1,T1,T2)
       Uonelink[ja,jb]=dot(eigenvector_bc_single[:,ja,jb],Amatrix*eigenvector_bc_single[:,ja+1,jb])/abs(dot(eigenvector_bc_single[:,ja,jb],Amatrix*eigenvector_bc_single[:,ja+1,jb]))
    end
    
    for ja in 1:Nq+1, jb in 1:Nq
        Amatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T2,T1,T2)
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







function calculate_energy(Nq::Int,wave::Vector{Vector{Int}},scale::Float64,ϕ::Float64,flux::Float64,input_DensityMatrix::Vector{Matrix{ComplexF64}},constq::Float64,overlapmatrix::Array{ComplexF64,4})

  
    β=4*flux/(√3*scale^2)
    mass=0.5;
    dimension=length(wave)
    
    b1=scale*[0,1]
    b2=scale*[√3/2,-1/2]

    
    
    T1=b1/(Nq)
    T2=b2/(Nq)
    b1T=Int.(round.(inv([T1 T2])*b1))
    b2T=Int.(round.(inv([T1 T2])*b2))
    
    
    chern_allowedq=Vector{Int64}[]
    for ja in 0:Nq-1,jb in 0:Nq-1
        push!(chern_allowedq,[ja,jb])
    end
    

    loop_dic=construct_loop_dic(wave)
    
    Energy_Matrix=[zeros(ComplexF64,length(wave),length(wave)) for _ in 1:Nq^2]
    
    Threads.@threads for ja in eachindex(chern_allowedq)
        k=[T1 T2]*chern_allowedq[ja]
        chern_Ham=zeros(ComplexF64,dimension,dimension)
        chern_MoirePo=zeros(ComplexF64,dimension,dimension)
        for jb in eachindex(wave)
            chern_Ham[jb,jb]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2)^2/(2*mass)
        end
    
        for jc in eachindex(wave)
            k1=k+wave[jc][1]*T1+wave[jc][2]*T2
            pos=findfirst(item->item==wave[jc]-b1T,wave)
            if pos≠nothing
           chern_MoirePo[jc,pos]=V0*exp(im*ϕ)*overlap(k1,-b1,β)
            end
            
          
            pos=findfirst(item->item==wave[jc]+b2T+b1T,wave)
            if pos≠nothing
                chern_MoirePo[jc,pos]=V0*exp(im*ϕ)*overlap(k1,b2+b1,β)
            end
        
            pos=findfirst(item->item==wave[jc]-b2T,wave)
            if pos≠nothing
                chern_MoirePo[jc,pos]=V0*exp(im*ϕ)*overlap(k1,-b2,β)
            end
         end


      chern_MoirePo=chern_MoirePo+chern_MoirePo'
      HFmatrix=Construct_HFmatrix(loop_dic,ja,chern_allowedq[ja],allowedq,T1,T2,Nq,wave,input_DensityMatrix,constq,overlapmatrix)
      Energy_Matrix[ja]=1/2*HFmatrix+chern_MoirePo+chern_Ham
    
    end
    
    
   energy=0
   for ja in 1:Nq^2
       energy+=tr(Energy_Matrix[ja]*input_DensityMatrix[ja])
   end


   return energy

end