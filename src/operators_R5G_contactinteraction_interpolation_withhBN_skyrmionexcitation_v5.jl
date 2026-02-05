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





function sendtomesh(Minv::Matrix{Int64},q1::Vector{Int})::Vector{Int}
  Qvec=q1'*inv(Minv)
  return q1 .-vec((Int.(floor.(round.(Qvec,digits=5)))*Minv)')
end




function Geometry(geonum::Int64)
  if geonum==1
      Nx=3;
      Ny=3;
      l1=[3,0]
      l2=[0,3]
  end

  if geonum==2
      Nx=3;
      Ny=9;
      l1=[1,1]*3
      l2=[-1,2]*3
  end

   if geonum==3
      Nx=6;
      Ny=6;
      l1=[6,0]
      l2=[0,6]
  end
  
  

  return Nx,Ny,l1,l2
end



function triangle_initial_Densitymatrix(NL::Int,θ::Float64,geonum::Int64,gcutoff::Float64,uD::Float64,λ::Float64,
                                        V0_hBN::Float64,V1_hBN::Float64,ψ_hBN::Float64,V2_scalar::Float64,ϕ::Float64)
   
    aGr=0.246
    ϵ=0.2504/aGr-1 #This is the normal one
  
    G1=2π/aGr*[1,-1/√3]
    G2=2π/aGr*[0,2/√3]
    Nx,Ny,l1,l2=Geometry(geonum)

    Rθ=[cos(θ) -sin(θ);sin(θ) cos(θ)]
    b1=(G1-(1+ϵ)^(-1)*Rθ*G1)
    b2=(G2-(1+ϵ)^(-1)*Rθ*G2)
    a1m=inv([b1';b2'])*[2π,0]
    a2m=inv([b1';b2'])*[0,2π]
    
    
      
    L1=l1[1]*a1m+l1[2]*a2m;
    L2=l2[1]*a1m+l2[2]*a2m;

    Area=abs(L1[1]*L2[2]-L2[1]*L1[2]);
    
    Rotminus90=[0 1;-1 0]
    T1=2*π/Area*Rotminus90*L2
    T2=-2*π/Area*Rotminus90*L1

    
    am=norm(a1m);
    

    b1T=Int.(round.(inv([T1 T2])*b1))
    b2T=Int.(round.(inv([T1 T2])*b2))




   
 
    
 
    
    

    
    

    
    
    
   allowedq=Vector{Int64}[]
  for ja in 0:Nx-1,jb in 0:Ny-1
      push!(allowedq,[ja,jb])
  end
  
  for ja in 1:Nx*Ny
    allowedq[ja]=sendtomesh([b1T';b2T'],allowedq[ja])
  end
    
  allowedq=unique(allowedq)
  if length(allowedq)<Nx*Ny
   throw("there is an error in allowedq")
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
   
    
    
    spinor_set=Matrix{Vector{ComplexF64}}(undef,length(allowedq),length(wave))
    for ja in eachindex(allowedq), jb in eachindex(wave)
      hh=(1-λ)*get_Ham_Holomorphic([T1 T2]*(allowedq[ja]+wave[jb]),uD,1,1,NL)+(λ)*get_Ham([T1 T2]*(allowedq[ja]+wave[jb]),uD,1,1,NL)
      v1=eigen(hh).vectors[:,NL+1]
      if uD>0.0
        v1angle=angle(v1[2*NL])
      else
        throw("there is an error")
      end
      spinor_set[ja,jb]=v1*exp(-im*v1angle)/norm(v1)
    end



     overlapmatrix=zeros(ComplexF64,length(allowedq),length(wave),length(allowedq),length(wave))
    for ja in eachindex(allowedq)
      for jb in eachindex(wave), jc in eachindex(allowedq), jd in eachindex(wave)
       overlapmatrix[ja,jb,jc,jd]=spinor_set[ja,jb]'*spinor_set[jc,jd]
      end
    end



    
     
    single_Ham=[zeros(ComplexF64,dimension,dimension) for _ in eachindex(allowedq)]
    single_MoirePo=[zeros(ComplexF64,dimension,dimension) for _ in eachindex(allowedq)]
    single_eigenvalue=[zeros(Float64,dimension) for _ in eachindex(allowedq)]
    single_eigenvector=[zeros(ComplexF64,dimension,dimension) for _ in eachindex(allowedq)]

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
    
    Threads.@threads for ja in eachindex(allowedq)
      
      
      k=allowedq[ja][1]*T1+allowedq[ja][2]*T2
    
      for jb in eachindex(wave)
          hh=(1-λ)*get_Ham_Holomorphic([T1 T2]*(allowedq[ja]+wave[jb]),uD,1,1,NL)+(λ)*get_Ham([T1 T2]*(allowedq[ja]+wave[jb]),uD,1,1,NL)
   
         
          single_Ham[ja][jb,jb]=real(eigen(hh).values[NL+1])

         
      end


     
     for jc in eachindex(wave)
    
        pos=findfirst(item->item==wave[jc]-b1T,wave)
        if pos≠nothing
       
          single_MoirePo[ja][jc,pos]+=V1_hBN*exp(-im*ψ_hBN)*(spinor_set[ja,jc]'*op_1*spinor_set[ja,pos])
          single_MoirePo[ja][jc,pos]+=V2_scalar*exp(-im*ϕ)*(spinor_set[ja,jc]'*op_5*spinor_set[ja,pos])
        end
        
        pos=findfirst(item->item==wave[jc]-b2T,wave)
        if pos≠nothing

            single_MoirePo[ja][jc,pos]+=V1_hBN*exp(-im*ψ_hBN)*(spinor_set[ja,jc]'*op_2*spinor_set[ja,pos])
            single_MoirePo[ja][jc,pos]+=V2_scalar*exp(-im*ϕ)*(spinor_set[ja,jc]'*op_5*spinor_set[ja,pos])
        end


        pos=findfirst(item->item==wave[jc]+(b2T+b1T),wave)
        if pos≠nothing

            single_MoirePo[ja][jc,pos]+=V1_hBN*exp(-im*ψ_hBN)*(spinor_set[ja,jc]'*op_3*spinor_set[ja,pos])
            single_MoirePo[ja][jc,pos]+=V2_scalar*exp(-im*ϕ)*(spinor_set[ja,jc]'*op_5*spinor_set[ja,pos])
        end
    
        single_MoirePo[ja][jc,jc]+=V0_hBN/2*(spinor_set[ja,jc]'*op_4*spinor_set[ja,jc])
     end


    
     single_MoirePo[ja]=single_MoirePo[ja]+single_MoirePo[ja]'
     FFF=eigen(single_MoirePo[ja]+single_Ham[ja])
    
      single_eigenvalue[ja]=real(FFF.values)
      single_eigenvector[ja]=FFF.vectors
    
    end


    input_DensityMatrix=[zeros(ComplexF64,dimension,dimension) for _ in eachindex(allowedq)]
 
    
     for ja in eachindex(allowedq)
      A=randn(ComplexF64,length(wave),length(wave))
      input_DensityMatrix[ja]=(A+A')*1.0
     end

    

    
    return overlapmatrix, wave, input_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,single_eigenvector,allowedq, T1, T2, a1m, a2m, b1,b2,spinor_set,Area,b1T,b2T
      

       
end










function construct_loop_dic(wave::Vector{Vector{Int}},allowedq::Vector{Vector{Int}},T1,T2,ϵr,constq)


    loop_dic_Fock=Vector{Int}[]

    for ja in eachindex(wave), jb in eachindex(wave), jc in eachindex(wave), jd in 1:ja
      if wave[ja]+wave[jb]==wave[jc]+wave[jd]
      
        push!(loop_dic_Fock,[ja,jb,jc,jd])
        
      end
    end

    loop_dic_Fock_val=zeros(Float64,length(allowedq),length(allowedq),length(loop_dic_Fock))

    Threads.@threads for k1 in eachindex(allowedq)
        for k2 in eachindex(allowedq)
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
                              wave::Vector{Vector{Int64}},
                               input_DensityMatrix::Vector{Matrix{ComplexF64}},single_Ham::Vector{Matrix{ComplexF64}},
                               single_MoirePo::Vector{Matrix{ComplexF64}},overlapmatrix::Array{ComplexF64,4},
                               energy_input::Float64,filling::Int,Area::Float64)
  
  
   dimension=length(wave)
  HartreeMatrix=[zeros(ComplexF64,dimension,dimension) for _ in eachindex(allowedq)]
  FockMatrix=[zeros(ComplexF64,dimension,dimension) for _ in eachindex(allowedq)]
  output_DensityMatrix=Vector{Matrix{ComplexF64}}(undef,length(allowedq))
  DeltaMatrix=Vector{Matrix{ComplexF64}}(undef,length(allowedq))
  NewDensityMatrix=[zeros(ComplexF64,dimension,dimension) for _ in eachindex(allowedq)]
  HF_eigenvalue=Vector{Vector{Float64}}(undef,length(allowedq))
  HF_eigenvector=Vector{Matrix{ComplexF64}}(undef,length(allowedq))
 

  Threads.@threads for jk in eachindex(allowedq)
    Fk = FockMatrix[jk]
    for jk1 in eachindex(allowedq)
        dmk = input_DensityMatrix[jk1]
        lfv=loop_dic_Fock_val[jk,jk1,:]
        opf=overlapmatrix[jk1,:,jk,:]
        opm=overlapmatrix[jk,:,jk1,:]

        for (ja,vvs) in pairs(loop_dic_Fock)
          Fk[vvs[1],vvs[4]]+=dmk[vvs[3],vvs[2]]*opm[vvs[1],vvs[3]]*opf[vvs[2],vvs[4]]*lfv[ja]
        end
     end
  end






   Hartree_Density=sum([input_DensityMatrix[ja] .* transpose(overlapmatrix[ja,:,ja,:]) for ja in eachindex(allowedq)])



  Threads.@threads for jk in eachindex(allowedq)
    Hk=HartreeMatrix[jk]
    oph=(overlapmatrix[jk,:,jk,:])   
      for (ja,vvs) in pairs(loop_dic_Hartree)
        Hk[vvs[1],vvs[3]]+=Hartree_Density[vvs[4],vvs[2]]*oph[vvs[1],vvs[3]]*loop_dic_Hartree_val[ja]
      end
  end


 for ja in eachindex(allowedq)
    HartreeMatrix[ja]=(HartreeMatrix[ja]+HartreeMatrix[ja]'-real(Diagonal(HartreeMatrix[ja])))/Area
    FockMatrix[ja]=(FockMatrix[ja]+FockMatrix[ja]'-real(Diagonal(FockMatrix[ja])))/Area
  end
 

 for ja in eachindex(allowedq)
   FFF=eigen(single_MoirePo[ja]+single_Ham[ja]+HartreeMatrix[ja]-FockMatrix[ja])
   HF_eigenvalue[ja]=real(FFF.values)
   HF_eigenvector[ja]=FFF.vectors
 end

 bound=(sort(reduce(vcat,HF_eigenvalue))[filling*length(allowedq)+1]+sort(reduce(vcat,HF_eigenvalue))[filling*length(allowedq)])/2

  for ja in eachindex(allowedq)
    
       for jd in eachindex(HF_eigenvalue[ja])
          if HF_eigenvalue[ja][jd]<bound
             NewDensityMatrix[ja]+=HF_eigenvector[ja][:,jd]*(HF_eigenvector[ja][:,jd])'
          end
       end
        NewDensityMatrix[ja]=0.5*(NewDensityMatrix[ja]'+ NewDensityMatrix[ja])
       DeltaMatrix[ja]=NewDensityMatrix[ja]-input_DensityMatrix[ja]
       mix_ratio=rand()
       output_DensityMatrix[ja]=mix_ratio*input_DensityMatrix[ja]+(1-mix_ratio)*NewDensityMatrix[ja]
  end


  
  e1=0.0
  for ja in eachindex(allowedq)
    e1+=sum(abs2,DeltaMatrix[ja])
  end
  eout=e1/length(allowedq)
  
  
  energy=0
   for ja in eachindex(allowedq)
       ss=single_MoirePo[ja]+single_Ham[ja]+0.5*HartreeMatrix[ja]-0.5*FockMatrix[ja]
    
        energy+=real(sum(ss .* transpose(input_DensityMatrix[ja])))
   end

   energy_change=real(energy-energy_input)



 return  eout,energy_change,output_DensityMatrix,DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,HartreeMatrix,FockMatrix
end





function iteration_loop(initial_DensityMatrix::Vector{Matrix{ComplexF64}},
                       allowedq::Vector{Vector{Int64}},T1::Vector{Float64},T2::Vector{Float64},
                       wave::Vector{Vector{Int64}},single_Ham::Vector{Matrix{ComplexF64}},
                       single_MoirePo::Vector{Matrix{ComplexF64}},constq::Float64,ϵr::Float64,overlapmatrix::Array{ComplexF64,4},filling::Int,Area::Float64)
    eout=1.0
    itcount=0
    dimension=length(wave)

     loop_dic_Fock,loop_dic_Fock_val,loop_dic_Hartree,loop_dic_Hartree_val=construct_loop_dic(wave,allowedq,T1,T2,ϵr,constq)
   
    
    HF_eigenvalue=Vector{Vector{Float64}}(undef,length(allowedq))
    HF_eigenvector=Vector{Matrix{ComplexF64}}(undef,length(allowedq))
    HartreeMatrix=[zeros(ComplexF64,dimension,dimension) for _ in eachindex(allowedq)]
    FockMatrix=[zeros(ComplexF64,dimension,dimension) for _ in eachindex(allowedq)]
  
    DIIS_size=3
    DIIS_input_DensityMatrix=Vector{Vector{Matrix{ComplexF64}}}(undef,DIIS_size)
    DIIS_input_DeltaMatrix=Vector{Vector{Matrix{ComplexF64}}}(undef,DIIS_size)
    input_DensityMatrix=initial_DensityMatrix
    bad_count=0
    energy=0.0
    energy_change=0.0
  

    while (eout>1*10^(-22)) || (bad_count<4) || (energy_change>1*10^(-10))

      if eout<1*10^(-22)
       bad_count+=1
      end
      
      tic=time()

      if (itcount>100 && abs(eout)>10^(-2)) || (itcount>30 && abs(eout)<10^(-10))
      
        dmk=implement_DIIS(DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,allowedq,DIIS_size)
        if dmk==0
            itcount=0
            dmk=[zeros(ComplexF64,dimension,dimension) for _ in eachindex(allowedq)]
            for ja in eachindex(allowedq)
              A=randn(dimension,dimension)+im*randn(dimension,dimension)
              dmk[ja]+=(A+A')*0.01
            end
        end
      

       eout,energy_change,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,DIIS_size)+1],HF_eigenvalue,HF_eigenvector,energy,HartreeMatrix,FockMatrix=Construct_DensityMatrix(loop_dic_Fock,loop_dic_Fock_val,loop_dic_Hartree,loop_dic_Hartree_val,allowedq,wave,dmk,single_Ham,single_MoirePo,overlapmatrix,energy,filling,Area)
       
        DIIS_input_DensityMatrix[mod(itcount,DIIS_size)+1]=dmk
        input_DensityMatrix=output_DensityMatrix
        println("using DIIS")
       
      else
  
      

        eout,energy_change,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,DIIS_size)+1],HF_eigenvalue,HF_eigenvector,energy,HartreeMatrix,FockMatrix=Construct_DensityMatrix(loop_dic_Fock,loop_dic_Fock_val,loop_dic_Hartree,loop_dic_Hartree_val,allowedq,wave,input_DensityMatrix,single_Ham,single_MoirePo,overlapmatrix,energy,filling,Area)
        DIIS_input_DensityMatrix[mod(itcount,DIIS_size)+1]=input_DensityMatrix
        input_DensityMatrix=output_DensityMatrix
        
       
         
     

      end

    



      itcount+=1
     
      toc=time()
      println(toc-tic,"eout=$eout","energy_change=$energy_change","itcount=$itcount")
      flush(stdout)
     
    
  end







  return DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,eout,HartreeMatrix,FockMatrix

end











function shift_vector(vector_toshift::Array{ComplexF64,3},shiftamount::Vector{Int},wave::Vector{Vector{Int}},NL::Int,filling::Int)
    
  
  shifted_vector=zeros(ComplexF64,2*NL,length(wave),filling)
    for ja in eachindex(wave)
        pos=findfirst(item->item==wave[ja]+shiftamount,wave)
        if pos≠nothing
        shifted_vector[:,ja,:]=vector_toshift[:,pos,:]
        end
    
    end

  shifted_vector=reshape(shifted_vector,2*NL*length(wave),filling)
  for ja in 1:filling
   shifted_vector[:,ja]=shifted_vector[:,ja]/norm(shifted_vector[:,ja])
  end

    return shifted_vector

end



function triangle_chern(wave::Vector{Vector{Int}},
       allowedq::Vector{Vector{Int}},
       HF_eigenvector::Vector{Matrix{ComplexF64}},
       single_eigenvector::Vector{Matrix{ComplexF64}},spinor_set::Matrix{Vector{ComplexF64}},
       geonum::Int,b1T::Vector{Int},b2T::Vector{Int},NL::Int,filling::Int)

   Nx,Ny,_,_=Geometry(geonum)
    

    
    
    chern_allowedq=Vector{Int64}[]
    for ja in 0:Nx,jb in 0:Ny
        push!(chern_allowedq,[ja,jb])
    end

    PW_HF_eigenvector=[zeros(ComplexF64,2*NL,length(wave),filling) for _ in eachindex(allowedq)]
    PW_single_eigenvector=[zeros(ComplexF64,2*NL,length(wave),filling) for _ in eachindex(allowedq)]

    for ja in eachindex(allowedq), jb in eachindex(wave), jc in 1:filling
        PW_HF_eigenvector[ja][:,jb,jc]=HF_eigenvector[ja][jb,jc]*spinor_set[ja,jb]
        PW_single_eigenvector[ja][:,jb,jc]=single_eigenvector[ja][jb,jc]*spinor_set[ja,jb]
    end



  

    eigenvector_intermediate_bc=Vector{Array{ComplexF64,2}}(undef,(Nx+1)*(Ny+1))
    eigenvector_intermediate_single=Vector{Array{ComplexF64,2}}(undef,(Nx+1)*(Ny+1))

    Threads.@threads for ja in eachindex(chern_allowedq)
      
      
      kbraket=sendtomesh([b1T';b2T'],chern_allowedq[ja])
      kbraket_pos=findfirst(item->item==kbraket,allowedq)
      if kbraket==chern_allowedq[ja]
        eigenvector_intermediate_bc[ja]=reshape(PW_HF_eigenvector[kbraket_pos],2*NL*length(wave),filling)
        eigenvector_intermediate_single[ja]=reshape(PW_single_eigenvector[kbraket_pos],2*NL*length(wave),filling)
      else
        eigenvector_intermediate_bc[ja]=shift_vector(PW_HF_eigenvector[kbraket_pos],chern_allowedq[ja]-kbraket,wave,NL,filling)
        eigenvector_intermediate_single[ja]=shift_vector(PW_single_eigenvector[kbraket_pos],chern_allowedq[ja]-kbraket,wave,NL,filling)

      end  
    end

    eigenvector_bc=zeros(ComplexF64,length(wave)*2*NL,filling,Nx+1,Ny+1)
    eigenvector_bc_single=zeros(ComplexF64,length(wave)*2*NL,filling,Nx+1,Ny+1)
    for ja in eachindex(chern_allowedq)
       eigenvector_bc[:,:,chern_allowedq[ja][1]+1,chern_allowedq[ja][2]+1]=eigenvector_intermediate_bc[ja]
       eigenvector_bc_single[:,:,chern_allowedq[ja][1]+1,chern_allowedq[ja][2]+1]=eigenvector_intermediate_single[ja]
    end
    
    
    Uonelink=zeros(ComplexF64,Nx,Ny+1)
    Utwolink=zeros(ComplexF64,Nx+1,Ny)
    
    
    for ja in 1:Nx, jb in 1:Ny+1
       Uonelink[ja,jb]=det(eigenvector_bc[:,:,ja,jb]'*eigenvector_bc[:,:,ja+1,jb])/abs(det(eigenvector_bc[:,:,ja,jb]'*eigenvector_bc[:,:,ja+1,jb]))
    end

    
    
    for ja in 1:Nx+1, jb in 1:Ny
      Utwolink[ja,jb]=det(eigenvector_bc[:,:,ja,jb]'*eigenvector_bc[:,:,ja,jb+1])/abs(det(eigenvector_bc[:,:,ja,jb]'*eigenvector_bc[:,:,ja,jb+1]))
    end
    
    
   


    Flink=zeros(ComplexF64,Nx,Ny)
    for ja in 1:Nx, jb in 1:Ny
     Flink[ja,jb]=log(Uonelink[ja,jb]*Utwolink[ja+1,jb]/(Uonelink[ja,jb+1]*Utwolink[ja,jb]))
    end
    chern=sum(Flink)/(2*π*im)
    aveF=sum(Flink)/(Nx*Ny)
    uniform=0
    for ja in 1:Nx, jb in 1:Ny
        uniform+=(imag(Flink[ja,jb])-imag(aveF))^2*(Nx*Ny)/(2π)^2
    end


    Uonelink=zeros(ComplexF64,Nx,Ny+1)
    Utwolink=zeros(ComplexF64,Nx+1,Ny)
  
    
  for ja in 1:Nx, jb in 1:Ny+1
       Uonelink[ja,jb]=det(eigenvector_bc_single[:,:,ja,jb]'*eigenvector_bc_single[:,:,ja+1,jb])/abs(det(eigenvector_bc_single[:,:,ja,jb]'*eigenvector_bc_single[:,:,ja+1,jb]))
    end

    
    
    for ja in 1:Nx+1, jb in 1:Ny
      Utwolink[ja,jb]=det(eigenvector_bc_single[:,:,ja,jb]'*eigenvector_bc_single[:,:,ja,jb+1])/abs(det(eigenvector_bc_single[:,:,ja,jb]'*eigenvector_bc_single[:,:,ja,jb+1]))
    end
    

   

    Flink_single=zeros(ComplexF64,Nx,Ny)
    for ja in 1:Nx, jb in 1:Ny
     Flink_single[ja,jb]=log(Uonelink[ja,jb]*Utwolink[ja+1,jb]/(Uonelink[ja,jb+1]*Utwolink[ja,jb]))
    end
    chern_single=sum(Flink_single)/(2*π*im)
    aveF=sum(Flink_single)/(Nx*Ny)
    uniform_single=0
    for ja in 1:Nx, jb in 1:Ny
        uniform_single+=(imag(Flink_single[ja,jb])-imag(aveF))^2*(Nx*Ny)/(2π)^2
    end
    
    

    
 
      


    
      
   return chern,Flink,chern_single,Flink_single,uniform,uniform_single

end



function implement_DIIS(DIIS_input_projector::Vector{Vector{Matrix{ComplexF64}}},DIIS_input_DeltaMatrix::Vector{Vector{Matrix{ComplexF64}}},allowedq::Vector{Vector{Int}},DIIS_size::Int)



      Bmatrix=zeros(ComplexF64,DIIS_size+1,DIIS_size+1)
      for ja in 1:DIIS_size
       Bmatrix[ja,DIIS_size+1]=1
       Bmatrix[DIIS_size+1,ja]=1
      end
  
      for ja in 1:DIIS_size,jb in 1:DIIS_size
          for jc in eachindex(allowedq)
             #Bmatrix[ja,jb]+=real(tr((DIIS_input_DeltaMatrix[ja][jc])'*(DIIS_input_DeltaMatrix[jb][jc])))
             Bmatrix[ja,jb]+=real(dot(DIIS_input_DeltaMatrix[ja][jc],DIIS_input_DeltaMatrix[jb][jc]))
          end
      end

      inB=safe_inverse(Bmatrix)
      if inB≠0
         onh=zeros(Float64,DIIS_size+1)
         onh[end]=1
         coeff=inB* onh
         #AA=coeff[1]*(DIIS_input_projector[1]+DIIS_input_DeltaMatrix[1])+coeff[2]*(DIIS_input_projector[2]+DIIS_input_DeltaMatrix[2])+coeff[3]*(DIIS_input_projector[3]+DIIS_input_DeltaMatrix[3])
         
        AA =coeff[1]*(DIIS_input_projector[1]+DIIS_input_DeltaMatrix[1])           # alloc once
        @inbounds for ja in 2:DIIS_size
          AA+=coeff[ja]*(DIIS_input_projector[ja]+DIIS_input_DeltaMatrix[ja])  # dmk += coeff[ja] * projector[ja]
        end
         return AA
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