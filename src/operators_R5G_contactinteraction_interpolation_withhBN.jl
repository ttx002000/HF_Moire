using LinearAlgebra
using Arpack
using Combinatorics
using Random



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








function triangle_initial_Densitymatrix(NL::Int,θ::Float64,Nq::Int64,gcutoff::Float64,uD::Float64,λ::Float64)
   
    aGr=0.246
    ϵ=0.2504/aGr-1 #This is the normal one
  
    G1=2π/aGr*[1,-1/√3]
    G2=2π/aGr*[0,2/√3]
    
    Rθ=[cos(θ) -sin(θ);sin(θ) cos(θ)]
    b1=(G1-(1+ϵ)^(-1)*Rθ*G1)
    b2=(G2-(1+ϵ)^(-1)*Rθ*G2)
    T1=b1/Nq
    T2=b2/Nq

    a1m=inv([b1';b2'])*[2π,0]
    a2m=inv([b1';b2'])*[0,2π]
    am=norm(a1m);
    

    b1T=Int.(round.(inv([T1 T2])*b1))
    b2T=Int.(round.(inv([T1 T2])*b2))



    V0=28.9
    V1=0.0
    ψ=-0.29

      
    
 
    
    

    
    

    
    
    
    allowedq=Vector{Int64}[]
    for ja in 0:Nq-1,jb in 0:Nq-1
        push!(allowedq,[ja,jb])
    end

    
    
    
    wave=Vector{Int64}[]
    cutoff=18
    cutoffstandard=gcutoff*norm(b1)
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        gtest=ja*b1+jb*b2;
        if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
            push!(wave,ja*b1T+jb*b2T)
        end
    end
    dimension=length(wave)
   
    
    
    spinor_set=Matrix{Vector{ComplexF64}}(undef,Nq^2,length(wave))
    for ja in 1:Nq^2, jb in eachindex(wave)
      hh=(1-λ)*get_Ham_Holomorphic([T1 T2]*(allowedq[ja]+wave[jb]),uD,1,1,NL)+(λ)*get_Ham([T1 T2]*(allowedq[ja]+wave[jb]),uD,1,1,NL)
      v1=eigen(hh).vectors[:,NL+1]
      if uD>0.0
        v1angle=angle(v1[2*NL])
      else
        throw("there is an error")
      end
      spinor_set[ja,jb]=v1*exp(-im*v1angle)/norm(v1)
    end



     overlapmatrix=zeros(ComplexF64,Nq^2,length(wave),Nq^2,length(wave))
    for ja in 1:Nq^2
      for jb in eachindex(wave), jc in 1:Nq^2, jd in eachindex(wave)
       overlapmatrix[ja,jb,jc,jd]=spinor_set[ja,jb]'*spinor_set[jc,jd]
      end
    end



    
     
    single_Ham=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
    single_MoirePo=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
    single_eigenvalue=[zeros(Float64,dimension) for _ in 1:Nq^2]
    single_eigenvector=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]

         ω=exp(im*2π/3)
       op_1=zeros(ComplexF64,2*NL,2*NL)
       op_1[1:2,1:2]=[1 1;ω ω]

       op_2=zeros(ComplexF64,2*NL,2*NL)
       op_2[1:2,1:2]=[1 ω^2;ω^2 ω]

       op_3=zeros(ComplexF64,2*NL,2*NL)
       op_3[1:2,1:2]=[1 ω;1 ω]

       
       op_4=zeros(ComplexF64,2*NL,2*NL)
       op_4[1:2,1:2]=[1 0; 0 1]
    
    Threads.@threads for ja in 1:Nq^2
      
      
      k=allowedq[ja][1]*T1+allowedq[ja][2]*T2
    
      for jb in eachindex(wave)
          hh=(1-λ)*get_Ham_Holomorphic([T1 T2]*(allowedq[ja]+wave[jb]),uD,1,1,NL)+(λ)*get_Ham([T1 T2]*(allowedq[ja]+wave[jb]),uD,1,1,NL)
          single_Ham[ja][jb,jb]=real(eigen(hh).values[NL+1])
      end


     
     for jc in eachindex(wave)
    
        pos=findfirst(item->item==wave[jc]-b1T,wave)
        if pos≠nothing
       
          single_MoirePo[ja][jc,pos]=V1*exp(-im*ψ)*(spinor_set[ja,jc]'*op_1*spinor_set[ja,pos])
        end
        
        pos=findfirst(item->item==wave[jc]-b2T,wave)
        if pos≠nothing
   
            single_MoirePo[ja][jc,pos]=V1*exp(-im*ψ)*(spinor_set[ja,jc]'*op_2*spinor_set[ja,pos])
        end


        pos=findfirst(item->item==wave[jc]+b2T+b1T,wave)
        if pos≠nothing
   
            single_MoirePo[ja][jc,pos]=V1*exp(-im*ψ)*(spinor_set[ja,jc]'*op_3*spinor_set[ja,pos])
        end
    
        single_MoirePo[ja][jc,jc]+=V0/2*(spinor_set[ja,jc]'*op_4*spinor_set[ja,jc])
     end


    
     single_MoirePo[ja]=single_MoirePo[ja]+single_MoirePo[ja]'
     FFF=eigen(single_MoirePo[ja]+single_Ham[ja])
    
      single_eigenvalue[ja]=real(FFF.values)
      single_eigenvector[ja]=FFF.vectors
    
    end


    input_DensityMatrix=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
 
    
     for ja in 1:Nq^2
      A=randn(ComplexF64,length(wave),length(wave))
      input_DensityMatrix[ja]=(A+A')*1.0
     end

    

    
    return overlapmatrix, wave, input_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,single_eigenvector,allowedq, T1, T2, a1m, a2m, b1,b2,spinor_set
      

       
end










function construct_loop_dic(wave::Vector{Vector{Int}},Nq,T1,T2,ϵr,constq)


    loop_dic_Fock=Vector{Int}[]

    for ja in eachindex(wave), jb in eachindex(wave), jc in eachindex(wave), jd in 1:ja
      if wave[ja]+wave[jb]==wave[jc]+wave[jd]
      
        push!(loop_dic_Fock,[ja,jb,jc,jd])
        
      end
    end

    loop_dic_Fock_val=zeros(Float64,Nq^2,Nq^2,length(loop_dic_Fock))

    Threads.@threads for k1 in 1:Nq^2
        for k2 in 1:Nq^2
      for waveset in eachindex(loop_dic_Fock)
      cc=Coulomb(allowedq[k2]+wave[loop_dic_Fock[waveset][3]]-allowedq[k1]-wave[loop_dic_Fock[waveset][1]],T1,T2)/ϵr+constq
      loop_dic_Fock_val[k1,k2,waveset]=cc
    end
    end
   end
   


  loop_dic_Hartree=Vector{Int}[]
    for ja in eachindex(wave), jb in eachindex(wave), jc in 1:ja, jd in eachindex(wave)
      if wave[ja]+wave[jb]==wave[jc]+wave[jd]
      
        push!(loop_dic_Hartree,[ja,jb,jc,jd])
      
      end
    end

      loop_dic_Hartree_val=zeros(Float64,length(loop_dic_Hartree))
      
      for waveset in eachindex(loop_dic_Hartree)
      cc=Coulomb(wave[loop_dic_Hartree[waveset][3]]-wave[loop_dic_Hartree[waveset][1]],T1,T2)/ϵr+constq
      loop_dic_Hartree_val[waveset]=cc
      end
 


    
   return loop_dic_Fock,loop_dic_Fock_val,loop_dic_Hartree,loop_dic_Hartree_val

end





function Coulomb(k::Vector{Int},T1::Vector{Float64},T2::Vector{Float64})::Float64
   D=25
   return k==[0,0] ? D*9047.5636 : tanh(norm([T1 T2]*k*D))/norm(k[1]*T1+k[2]*T2)*9047.5636
end


function Construct_DensityMatrix(loop_dic_Fock::Vector{Vector{Int}},loop_dic_Fock_val::Array{Float64},
                               loop_dic_Hartree::Vector{Vector{Int}},loop_dic_Hartree_val::Vector{Float64},
                              allowedq::Vector{Vector{Int}},
                               T1::Vector{Float64},T2::Vector{Float64},Nq::Int64,wave::Vector{Vector{Int64}},
                               input_DensityMatrix::Vector{Matrix{ComplexF64}},single_Ham::Vector{Matrix{ComplexF64}},
                               single_MoirePo::Vector{Matrix{ComplexF64}},constq::Float64,ϵr::Float64,overlapmatrix::Array{ComplexF64,4},
                               energy_input::Float64,filling::Int,Area::Float64)
  
  
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
        lfv=loop_dic_Fock_val[jk,jk1,:]
        opf=overlapmatrix[jk1,:,jk,:]
        opm=overlapmatrix[jk,:,jk1,:]

        for (ja,vvs) in pairs(loop_dic_Fock)
          Fk[vvs[1],vvs[4]]+=dmk[vvs[3],vvs[2]]*opm[vvs[1],vvs[3]]*opf[vvs[2],vvs[4]]*lfv[ja]
        end
     end
  end






   Hartree_Density=sum([input_DensityMatrix[ja] .* transpose(overlapmatrix[ja,:,ja,:]) for ja in 1:Nq^2])



  Threads.@threads for jk in 1:Nq^2
    Hk=HartreeMatrix[jk]
    oph=(overlapmatrix[jk,:,jk,:])   
      for (ja,vvs) in pairs(loop_dic_Hartree)
        Hk[vvs[1],vvs[3]]+=Hartree_Density[vvs[4],vvs[2]]*oph[vvs[1],vvs[3]]*loop_dic_Hartree_val[ja]
      end
  end


 for ja in 1:Nq^2
    HartreeMatrix[ja]=(HartreeMatrix[ja]+HartreeMatrix[ja]'-real(Diagonal(HartreeMatrix[ja])))/Area
    FockMatrix[ja]=(FockMatrix[ja]+FockMatrix[ja]'-real(Diagonal(FockMatrix[ja])))/Area
  end
 

 for ja in 1:Nq^2
   FFF=eigen(single_MoirePo[ja]+single_Ham[ja]+HartreeMatrix[ja]-FockMatrix[ja])
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
        NewDensityMatrix[ja]=0.5*(NewDensityMatrix[ja]'+ NewDensityMatrix[ja])
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
       ss=single_MoirePo[ja]+single_Ham[ja]+0.5*HartreeMatrix[ja]-0.5*FockMatrix[ja]
       energy+=real(tr(ss*output_DensityMatrix[ja]))
   end

   energy_change=real(energy-energy_input)



 return  eout,energy_change,output_DensityMatrix,DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,HartreeMatrix,FockMatrix
end





function iteration_loop(initial_DensityMatrix::Vector{Matrix{ComplexF64}},
                       allowedq::Vector{Vector{Int64}},T1::Vector{Float64},T2::Vector{Float64},
                       Nq::Int,wave::Vector{Vector{Int64}},single_Ham::Vector{Matrix{ComplexF64}},
                       single_MoirePo::Vector{Matrix{ComplexF64}},constq::Float64,ϵr::Float64,overlapmatrix::Array{ComplexF64,4},filling::Int,Area::Float64)
    eout=1.0
    itcount=0
    dimension=length(wave)
    #loop_dic=construct_loop_dic(wave)
     loop_dic_Fock,loop_dic_Fock_val,loop_dic_Hartree,loop_dic_Hartree_val=construct_loop_dic(wave,Nq,T1,T2,ϵr,constq)
   
    
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
  

    while (eout>1*10^(-20)) || (bad_count<4) || (energy_change>1*10^(-10))
      if eout<1*10^(-20)
       bad_count+=1
      end
      
      tic=time()

      if (itcount>100 && abs(eout)>10^(-2)) || (itcount>30 && abs(eout)<10^(-9))
      
        dmk=implement_DIIS(DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,Nq)
        if dmk==0
            itcount=0
            dmk=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
            for ja in 1:Nq^2
              A=randn(dimension,dimension)+im*randn(dimension,dimension)
              dmk[ja]+=(A+A')*0.01
            end
        end
      

       eout,energy_change,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue,HF_eigenvector,energy,HartreeMatrix,FockMatrix=Construct_DensityMatrix(loop_dic_Fock,loop_dic_Fock_val,loop_dic_Hartree,loop_dic_Hartree_val,allowedq,T1,T2,Nq,wave,dmk,single_Ham,single_MoirePo,constq,ϵr,overlapmatrix,energy,filling,Area)
       
        DIIS_input_DensityMatrix[mod(itcount,3)+1]=dmk
        input_DensityMatrix=output_DensityMatrix
        println("using DIIS")
       
      else
  
      

        eout,energy_change,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue,HF_eigenvector,energy,HartreeMatrix,FockMatrix=Construct_DensityMatrix(loop_dic_Fock,loop_dic_Fock_val,loop_dic_Hartree,loop_dic_Hartree_val,allowedq,T1,T2,Nq,wave,input_DensityMatrix,single_Ham,single_MoirePo,constq,ϵr,overlapmatrix,energy,filling,Area)
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





function metric(wavelist::Vector{Vector{Int64}},k::Vector{Int},
          q::Vector{Int},spinor_set::Matrix{Vector{ComplexF64}},Nq::Int)::Matrix{ComplexF64}
    Amatrix=zeros(ComplexF64,length(wavelist),length(wavelist))
 


    for ja in 1:length(wavelist)
 
    k1=findfirst(item->item==mod.(k,Nq),allowedq)
    g1=findfirst(item->item==wavelist[ja]+k-mod.(k,Nq),wavelist)
    k2=findfirst(item->item==mod.(k+q,Nq),allowedq)
    g2=findfirst(item->item==wavelist[ja]+k+q-mod.(k+q,Nq),wavelist)
      if g1≠nothing && g2≠nothing
         Amatrix[ja,ja]=spinor_set[k1,g1]'*spinor_set[k2,g2]
      end
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



function triangle_chern(Nq::Int,wave::Vector{Vector{Int}},
      allowedq::Vector{Vector{Int}},
       HF_eigenvector::Vector{Matrix{ComplexF64}},
       single_eigenvector::Vector{Matrix{ComplexF64}},spinor_set::Matrix{Vector{ComplexF64}})

  
    

    dimension=length(wave)
    
    
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
   
        Amatrix=metric(wave,[ja-1,jb-1],[1,0],spinor_set,Nq)
       Uonelink[ja,jb]=dot(eigenvector_bc[:,ja,jb],Amatrix*eigenvector_bc[:,ja+1,jb])/abs(dot(eigenvector_bc[:,ja,jb],Amatrix*eigenvector_bc[:,ja+1,jb]))
    end

    
    
    for ja in 1:Nq+1, jb in 1:Nq
        Amatrix=metric(wave,[ja-1,jb-1],[0,1],spinor_set,Nq)
     Utwolink[ja,jb]=dot(eigenvector_bc[:,ja,jb],Amatrix*eigenvector_bc[:,ja,jb+1])/abs(dot(eigenvector_bc[:,ja,jb],Amatrix*eigenvector_bc[:,ja,jb+1]))
    end
    
    dG=norm(T2)
    for ja in 1:Nq, jb in 1:Nq
      
        Amatrix=metric(wave,[ja-1,jb-1],[1,0],spinor_set,Nq)
        Bmatrix=metric(wave,[ja-1,jb-1],[0,1],spinor_set,Nq)
        Cmatrix=metric(wave,[ja-1,jb-1],[1,1],spinor_set,Nq)
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
        #Amatrix=metric(wave,NL,(ja-1)*T1+(jb-1)*T2,T1,T1,T2,uD)
        Amatrix=metric(wave,[ja-1,jb-1],[1,0],spinor_set,Nq)
       Uonelink[ja,jb]=dot(eigenvector_bc_single[:,ja,jb],Amatrix*eigenvector_bc_single[:,ja+1,jb])/abs(dot(eigenvector_bc_single[:,ja,jb],Amatrix*eigenvector_bc_single[:,ja+1,jb]))
    end
    
    for ja in 1:Nq+1, jb in 1:Nq
        #Amatrix=metric(wave,NL,(ja-1)*T1+(jb-1)*T2,T2,T1,T2,uD)
        Amatrix=metric(wave,[ja-1,jb-1],[0,1],spinor_set,Nq)
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
          return pinv(A, 0.1)  # Use pseudoinverse as an alternative
      else
          rethrow(e)  # If another error occurs, propagate it
      end
  end
end