using LinearAlgebra
BLAS.set_num_threads(1)
using Arpack
using Combinatorics
using Random
using SparseArrays


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








function triangle_initial_Densitymatrix(NL::Int,am::Float64,Nq::Int64,gcutoff::Float64,uD::Float64,λ::Float64,enlarge_factor::Int64,
                                        V0_hBN::Float64,V1_hBN::Float64,ψ_hBN::Float64,V2_scalar::Float64,ϕ::Float64,flux1::Float64,flux2::Float64)
   


    a1m=am*[1,0]
    a2m=am*[1/2,√3/2]

    b1=inv([a1m';a2m'])*[2π,0]
    b2=inv([a1m';a2m'])*[0,2π]
  

    T1=b1/Nq
    T2=b2/Nq

 
    shift_vector=flux1*T1+flux2*T2
    

    b1T=Int.(round.(inv([T1 T2])*b1))
    b2T=Int.(round.(inv([T1 T2])*b2))


    
    
    allowedq=Vector{Int64}[]
    for ja in 0:Nq-1,jb in 0:Nq-1
        push!(allowedq,[ja,jb])
    end

    
    
    
    wave=Vector{Int64}[]
    cutoff=Int(round(gcutoff))*5
    cutoffstandard=gcutoff*norm(b1)
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        gtest=ja*b1+jb*b2;
        if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
            push!(wave,ja*b1T+jb*b2T)
        end
    end
    dimension=length(wave)
   

        wave_diff=Vector{Int64}[]

    cutoffstandard_diff=gcutoff*norm(b1)*2
    cutoff=Int(round(gcutoff))*10
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
            gtest=ja*b1+jb*b2;
            if (gtest[1]^2+gtest[2]^2)<cutoffstandard_diff^2
                push!(wave_diff,ja*b1T+jb*b2T)
            end
    end
    println(length(wave_diff))
    
    
    spinor_set=Matrix{Vector{ComplexF64}}(undef,Nq^2,length(wave))
    for ja in 1:Nq^2, jb in eachindex(wave)
      hh=(1-λ)*get_Ham_Holomorphic([T1 T2]*(allowedq[ja]+wave[jb])-shift_vector,uD,1,1,NL)+(λ)*get_Ham([T1 T2]*(allowedq[ja]+wave[jb])-shift_vector,uD,1,1,NL)
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

       op_5=Matrix{Float64}(I,2*NL,2*NL)
    
    Threads.@threads for ja in 1:Nq^2
      
      
      k=allowedq[ja][1]*T1+allowedq[ja][2]*T2
    
      for jb in eachindex(wave)
          hh=(1-λ)*get_Ham_Holomorphic([T1 T2]*(allowedq[ja]+wave[jb])-shift_vector,uD,1,1,NL)+(λ)*get_Ham([T1 T2]*(allowedq[ja]+wave[jb])-shift_vector,uD,1,1,NL)
   
         
          single_Ham[ja][jb,jb]=real(eigen(hh).values[NL+1])

         
      end


     
     for jc in eachindex(wave)
    
        pos=findfirst(item->item==wave[jc]-b1T*enlarge_factor,wave)
        if pos≠nothing
       
          single_MoirePo[ja][jc,pos]+=V1_hBN*exp(-im*ψ_hBN)*(spinor_set[ja,jc]'*op_1*spinor_set[ja,pos])
          single_MoirePo[ja][jc,pos]+=V2_scalar*exp(-im*ϕ)*(spinor_set[ja,jc]'*op_5*spinor_set[ja,pos])
        end
        
        pos=findfirst(item->item==wave[jc]-b2T*enlarge_factor,wave)
        if pos≠nothing

            single_MoirePo[ja][jc,pos]+=V1_hBN*exp(-im*ψ_hBN)*(spinor_set[ja,jc]'*op_2*spinor_set[ja,pos])
            single_MoirePo[ja][jc,pos]+=V2_scalar*exp(-im*ϕ)*(spinor_set[ja,jc]'*op_5*spinor_set[ja,pos])
        end


        pos=findfirst(item->item==wave[jc]+(b2T+b1T)*enlarge_factor,wave)
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


    input_DensityMatrix=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
 
    
     for ja in 1:Nq^2
      A=randn(ComplexF64,length(wave),length(wave))
      input_DensityMatrix[ja]=(A+A')*1.0
     end

    

    
    return overlapmatrix, wave, wave_diff, input_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,single_eigenvector,allowedq, T1, T2, a1m, a2m, b1,b2,spinor_set
      

       
end


function Cdag(op_index::Int64,state_index::Int64,sign::Int64)::Tuple{Int64,Int64}
    if (sign==0) | ((2^(op_index-1)& state_index)==2^(op_index-1))
    return 0,0
    end
 
    new_sign=sign*(-1)^count_ones((2^(op_index-1)-1) & (state_index))
    new_state_index=state_index+2^(op_index-1)
   
    return new_sign,new_state_index
end





function Cann(op_index::Int64,state_index::Int64,sign::Int64)::Tuple{Int64,Int64}
    if (sign==0) | ((2^(op_index-1)& state_index)≠2^(op_index-1))
    return 0,0
    end
 
    new_sign=sign*(-1)^count_ones((2^(op_index-1)-1) & (state_index))
    new_state_index=state_index-2^(op_index-1)
   
    return new_sign,new_state_index
end




function Construct_Manybodymatrix(reduced_Vcol::Vector{ComplexF64},reduced_Vcoor::Vector{Vector{Int64}},reduced_sinval::Vector{ComplexF64},reduced_sincoor::Vector{Vector{Int64}},state_can,state_integer)

    

  
   
  Rows = [Vector{Int64}() for _ in 1:Threads.maxthreadid()]
  Cols = [Vector{Int64}() for _ in 1:Threads.maxthreadid()]
  Vals = [Vector{ComplexF64}() for _ in 1:Threads.maxthreadid()]

  Threads.@threads for jc in eachindex(reduced_Vcol)
    for jb in eachindex(state_can)
   
      (sign,state)=Cann(reduced_Vcoor[jc][3],state_integer[jb],1)
      (sign,state)=Cann(reduced_Vcoor[jc][4],state,sign)
      (sign,state)=Cdag(reduced_Vcoor[jc][2],state,sign)
      (sign,state)=Cdag(reduced_Vcoor[jc][1],state,sign)
 
       state_index=searchsortedfirst(state_integer,state)
       
      if (state*sign)≠0 &&  (state_integer[state_index]==state) 
       push!(Rows[Threads.threadid()],state_index)
       push!(Cols[Threads.threadid()],jb)
       push!(Vals[Threads.threadid()],sign*reduced_Vcol[jc]/2)
      end
   end

  end


  Threads.@threads for jc in eachindex(reduced_sinval)
    for jb in eachindex(state_can)
   
      (sign,state)=Cann(reduced_sincoor[jc][2],state_integer[jb],1)
      (sign,state)=Cdag(reduced_sincoor[jc][1],state,sign)
 
       state_index=searchsortedfirst(state_integer,state)
       
      if (state*sign)≠0 &&  (state_integer[state_index]==state) 
       push!(Rows[Threads.threadid()],state_index)
       push!(Cols[Threads.threadid()],jb)
       push!(Vals[Threads.threadid()],sign*reduced_sinval[jc])
      end
   end

  end




  matrix_index1=reduce(vcat,Rows)
  matrix_index2=reduce(vcat,Cols)
  matrix_value=reduce(vcat,Vals)
  
  Rows=nothing;
  Cols=nothing;
  Vals=nothing;
  GC.gc()

 
 

    Hsp = sparse(
        matrix_index1,
        matrix_index2,
        matrix_value,
        length(state_can),
        length(state_can),
    )

    H = Hermitian(Matrix(Hsp))

    F = eigen(H)

    MB_spectrum = F.values
    ζ = F.vectors





 return MB_spectrum,ζ

end


function Construct_MBstate(Nq::Int64,Nparticle::Int64,allowedq::Vector{Vector{Int64}},sector::Int64)


  MB_state=collect(combinations(1:Nq^2,Nparticle))
  MB_state_can_threads =
      [Vector{Int64}[] for _ in 1:Threads.maxthreadid()]

  MB_state_integer_threads =
      [Int64[] for _ in 1:Threads.maxthreadid()]
  
  
   Threads.@threads for ja in eachindex(MB_state)
       QN=sum(allowedq[MB_state[ja]])
       QN_M=[mod(QN[1],Nq),mod(QN[2],Nq)]
       if QN_M==allowedq[sector]
           push!(MB_state_can_threads[Threads.threadid()],MB_state[ja])
           push!(MB_state_integer_threads[Threads.threadid()],sum(2 .^ (MB_state[ja].-1)))
         end
   end

   MB_state_can=reduce(vcat,MB_state_can_threads)
   MB_state_integer=reduce(vcat,MB_state_integer_threads)
   MB_state_can_threads=nothing
   MB_state_integer_threads=nothing

   
  
     sortindex=sortperm(MB_state_integer)
     MB_state_integer=MB_state_integer[sortindex]
     MB_state_can=MB_state_can[sortindex]
   

  return MB_state_can, MB_state_integer


end


function do_ED(allowedq::Vector{Vector{Int}},wave::Vector{Vector{Int}},NL::Int,wave_diff::Vector{Vector{Int}},
               spinor_set::Matrix{Vector{ComplexF64}},HF_eigenvector::Vector{Matrix{ComplexF64}},Nq::Int,
               gateD::Float64,ϵr::Float64,Area::Float64,constq::Float64,Nparticle::Int)
  
      HF_eigenvector_pw=zeros(ComplexF64,length(allowedq),length(wave),2*NL)
      for ja in eachindex(allowedq), jb in eachindex(wave)
        HF_eigenvector_pw[ja,jb,:]=spinor_set[ja,jb]*HF_eigenvector[ja][jb,1]
      end


      Fmatrix=zeros(ComplexF64,Nq^2,Nq^2,length(wave_diff))

      diff_eigenvector=zeros(
          ComplexF64,
          Nq^2,
          length(wave),
          2*NL,
          length(wave_diff)
      )


      for jc in 1:length(wave_diff), jd in 1:length(wave)
          pos=findfirst(item->item==wave[jd]+wave_diff[jc],wave)
          if pos≠nothing 
            diff_eigenvector[:,pos,:,jc]=HF_eigenvector_pw[:,jd,:]
          
          end 
      end

      for ja in 1:Nq^2, jb in 1:Nq^2, jc in 1:length(wave_diff)
          Fmatrix[ja,jb,jc]=sum(conj(vec(HF_eigenvector_pw[ja,:,:])).*vec(diff_eigenvector[jb,:,:,jc]))
      end

      diff_eigenvector=nothing




      
      Vmatrix=zeros(ComplexF64,Nq^2,Nq^2,Nq^2,Nq^2)


      for jc in 1:Nq^2, jd in 1:(jc-1)
        for qmesh in 1:Nq^2
            k1mesh=[mod(allowedq[jc][1]+allowedq[qmesh][1],Nq),mod(allowedq[jc][2]+allowedq[qmesh][2],Nq)]
            k2mesh=[mod(allowedq[jd][1]-allowedq[qmesh][1],Nq),mod(allowedq[jd][2]-allowedq[qmesh][2],Nq)]
            k1mesh_pos=findfirst(item->item==k1mesh,allowedq)
            k2mesh_pos=findfirst(item->item==k2mesh,allowedq)

            for qg in 1:length(wave_diff)
        
            q_vec=allowedq[qmesh][1]*T1+allowedq[qmesh][2]*T2+wave_diff[qg][1]*b1+wave_diff[qg][2]*b2

            
            gk3pq_int=allowedq[jc]-k1mesh+(allowedq[qmesh]+wave_diff[qg])
            gk4mq_int=allowedq[jd]-k2mesh-(allowedq[qmesh]+wave_diff[qg])

            gk3pq_pos=findfirst(item->item==gk3pq_int,wave_diff)
            gk4mq_pos=findfirst(item->item==gk4mq_int,wave_diff)
            
            if k1mesh_pos≠nothing && k2mesh_pos≠nothing && gk3pq_pos≠nothing && gk4mq_pos≠nothing 
          
            Vmatrix[k1mesh_pos,k2mesh_pos,jc,jd]+=(Coulomb(allowedq[qmesh]+wave_diff[qg],T1,T2,gateD)/ϵr+constq)/Area*Fmatrix[k1mesh_pos,jc,gk3pq_pos]*Fmatrix[k2mesh_pos,jd,gk4mq_pos]
        
            end 

        end
        end
      end



      reduced_Vcol=ComplexF64[]
      reduced_Vcoor=Vector{Int64}[]

      for i in 1:Nq^2, j in 1:i-1, k in 1:Nq^2, p in 1:k-1
        if abs(Vmatrix[i,j,k,p])>10^(-10)
          push!(reduced_Vcol,2*Vmatrix[i,j,k,p]-2*Vmatrix[j,i,k,p])
          push!(reduced_Vcoor,[i,j,k,p])
        end
      end
      print(length(reduced_Vcol))


      reduced_sinval=ComplexF64[]
      reduced_sincoor=Vector{Int64}[]

      for i in 1:Nq^2
          transformed_singleHam=HF_eigenvector[i]'*(single_Ham[i]+single_MoirePo[i])*HF_eigenvector[i]
        push!(reduced_sinval,transformed_singleHam[1,1])
        push!(reduced_sincoor,[i,i])

      end
      print(length(reduced_sincoor))



      all_MB_state_can=Vector{Vector{Vector{Int64}}}(undef,length(allowedq))
    all_MB_state_integer=Vector{Vector{Int64}}(undef,length(allowedq))
 
    for ja in eachindex(allowedq)
        all_MB_state_can[ja],all_MB_state_integer[ja]=Construct_MBstate(Nq,Nparticle,allowedq,ja)

    end


    values_record=Vector{Vector{ComplexF64}}(undef,Nq^2)
    eig_vec_record=Vector{Matrix{ComplexF64}}(undef,Nq^2)

    for sector in 1:Nq^2
        if length(all_MB_state_can[sector])>0
            values_record[sector],eig_vec_record[sector]=Construct_Manybodymatrix(reduced_Vcol,reduced_Vcoor,reduced_sinval,reduced_sincoor,all_MB_state_can[sector],all_MB_state_integer[sector])
        end
    end


    G_dic_record=Matrix{Any}(undef,Nq^2,3)
    rho_mat_record=Matrix{Any}(undef,Nq^2,3)

   for sector in 1:Nq^2,  which_eig in 1:3
        G_dic_record[sector,which_eig]=get_G(all_MB_state_integer[sector],
                        eig_vec_record[sector][:,which_eig],
                        allowedq,
                        wave_diff,
                        Fmatrix,
                        Nq)

        rho_mat_record[sector,which_eig]=get_rho(all_MB_state_integer[sector],
                          eig_vec_record[sector][:,which_eig],
                          Nq)
   end

    return values_record,eig_vec_record,all_MB_state_can, all_MB_state_integer,Fmatrix,G_dic_record,rho_mat_record

end





function construct_loop_dic(wave::Vector{Vector{Int}},Nq,T1,T2,ϵr,constq,gateD)


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
      cc=Coulomb(allowedq[k2]+wave[loop_dic_Fock[waveset][3]]-allowedq[k1]-wave[loop_dic_Fock[waveset][1]],T1,T2,gateD)/ϵr+constq
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
      cc=Coulomb(wave[loop_dic_Hartree[waveset][3]]-wave[loop_dic_Hartree[waveset][1]],T1,T2,gateD)/ϵr+constq
      loop_dic_Hartree_val[waveset]=cc
      end
 


    
   return loop_dic_Fock,loop_dic_Fock_val,loop_dic_Hartree,loop_dic_Hartree_val

end





function Coulomb(k::Vector{Int},T1::Vector{Float64},T2::Vector{Float64},gateD::Float64)::Float64
  
   return k==[0,0] ? gateD*9047.5636 : tanh(norm([T1 T2]*k*gateD))/norm(k[1]*T1+k[2]*T2)*9047.5636
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
       mix_ratio=0.5
       output_DensityMatrix[ja]=mix_ratio*input_DensityMatrix[ja]+(1-mix_ratio)*NewDensityMatrix[ja]
  end


  
  eout=0.0
  for ja in 1:Nq^2
    eout+=sum(abs2,DeltaMatrix[ja])/Nq^2
  end
 
  
  
  energy=0
   for ja in 1:Nq^2
       ss=single_MoirePo[ja]+single_Ham[ja]+0.5*HartreeMatrix[ja]-0.5*FockMatrix[ja]
       energy+=real(sum(ss .* transpose(input_DensityMatrix[ja])))
   end

   energy_change=real(energy-energy_input)



 return  eout,energy_change,output_DensityMatrix,DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,HartreeMatrix,FockMatrix
end





function iteration_loop(initial_DensityMatrix::Vector{Matrix{ComplexF64}},
                       allowedq::Vector{Vector{Int64}},T1::Vector{Float64},T2::Vector{Float64},
                       Nq::Int,wave::Vector{Vector{Int64}},single_Ham::Vector{Matrix{ComplexF64}},
                       single_MoirePo::Vector{Matrix{ComplexF64}},constq::Float64,ϵr::Float64,overlapmatrix::Array{ComplexF64,4},
                       filling::Int,Area::Float64,gateD::Float64)
    eout=1.0
    itcount=0
    dimension=length(wave)
    #loop_dic=construct_loop_dic(wave)
     loop_dic_Fock,loop_dic_Fock_val,loop_dic_Hartree,loop_dic_Hartree_val=construct_loop_dic(wave,Nq,T1,T2,ϵr,constq,gateD)
   
    
    HF_eigenvalue=Vector{Vector{Float64}}(undef,Nq^2)
    HF_eigenvector=Vector{Matrix{ComplexF64}}(undef,Nq^2)
    HartreeMatrix=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
    FockMatrix=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
  
    DIIS_size=5
    DIIS_input_DensityMatrix=Vector{Vector{Matrix{ComplexF64}}}(undef,DIIS_size)
    DIIS_input_DeltaMatrix=Vector{Vector{Matrix{ComplexF64}}}(undef,DIIS_size)
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
  

    while (eout>1*10^(-18)) || (bad_count<4) || (energy_change>1*10^(-9))
      if eout<1*10^(-18)
       bad_count+=1
      else
        bad_count=0
      end
      
      tic=time()

      if (itcount>100 && abs(eout)>10^(-2)) || (itcount>30 && abs(eout)<10^(-8)) || diis_fire_once
      
        dmk=implement_DIIS(DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,Nq,DIIS_size)
        if dmk==0
            itcount=0
            dmk=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
            for ja in 1:Nq^2
              A=randn(dimension,dimension)+im*randn(dimension,dimension)
              dmk[ja]+=(A+A')*0.01
            end
        end
      

       eout,energy_change,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,DIIS_size)+1],HF_eigenvalue,HF_eigenvector,energy,HartreeMatrix,FockMatrix=Construct_DensityMatrix(loop_dic_Fock,loop_dic_Fock_val,loop_dic_Hartree,loop_dic_Hartree_val,allowedq,T1,T2,Nq,wave,dmk,single_Ham,single_MoirePo,constq,ϵr,overlapmatrix,energy,filling,Area)
       
        DIIS_input_DensityMatrix[mod(itcount,DIIS_size)+1]=dmk
        input_DensityMatrix=output_DensityMatrix
        println("using DIIS")
          if diis_fire_once
            diis_fire_once = false
            empty!(eout_hist)          # <-- yes: clear history after firing
            diis_cooldown = DIIS_COOLDOWN
           end
       
      else
  
      

        eout,energy_change,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,DIIS_size)+1],HF_eigenvalue,HF_eigenvector,energy,HartreeMatrix,FockMatrix=Construct_DensityMatrix(loop_dic_Fock,loop_dic_Fock_val,loop_dic_Hartree,loop_dic_Hartree_val,allowedq,T1,T2,Nq,wave,input_DensityMatrix,single_Ham,single_MoirePo,constq,ϵr,overlapmatrix,energy,filling,Area)
        DIIS_input_DensityMatrix[mod(itcount,DIIS_size)+1]=input_DensityMatrix
        input_DensityMatrix=output_DensityMatrix
        
       
         
     

      end

    



      itcount+=1
     
      toc=time()
      println(toc-tic,"eout=$eout","energy_change=$energy_change","itcount=$itcount")
      flush(stdout)

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



function implement_DIIS(DIIS_input_projector::Vector{Vector{Matrix{ComplexF64}}},DIIS_input_DeltaMatrix::Vector{Vector{Matrix{ComplexF64}}},Nq::Int,DIIS_size::Int)



     
      Bmatrix=zeros(Float64,DIIS_size+1,DIIS_size+1)
      for ja in 1:DIIS_size
       Bmatrix[ja,DIIS_size+1]=1
       Bmatrix[DIIS_size+1,ja]=1
      end
  
  
   
  
      for ja in 1:DIIS_size,jb in 1:DIIS_size
          for jc in 1:Nq^2
             Bmatrix[ja,jb]+=real(dot(DIIS_input_DeltaMatrix[ja][jc],DIIS_input_DeltaMatrix[jb][jc]))
          end
      end



      inB=safe_inverse(Bmatrix)
      if inB≠0
         onh=zeros(Float64,DIIS_size+1)
         onh[end]=1
         coeff=inB* onh
      
        dmk = zero.(DIIS_input_projector[1])          # alloc once
        @inbounds for ja in 1:DIIS_size
           dmk.+= coeff[ja].*(DIIS_input_projector[ja]+DIIS_input_DeltaMatrix[ja])  # dmk += coeff[ja] * projector[ja]
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
          return pinv(A, 0.001)  # Use pseudoinverse as an alternative
      else
          rethrow(e)  # If another error occurs, propagate it
      end
  end
end

















function get_G(state_integer::Vector{Int64},
               state_vector::Vector{ComplexF64},
               allowedq::Vector{Vector{Int64}},
               wave_diff::Vector{Vector{Int64}},
               Fmatrix::Array{ComplexF64,3},
               Nq::Int64)



    allowedq_dict=Dict{Tuple{Int64,Int64},Int64}()
    for ja in eachindex(allowedq)
        allowedq_dict[Tuple(allowedq[ja])]=ja
    end


    wave_diff_dict=Dict{Tuple{Int64,Int64},Int64}()
    for ja in eachindex(wave_diff)
        wave_diff_dict[Tuple(wave_diff[ja])]=ja
    end



    G_dic=Dict{Tuple{Int64,Int64},ComplexF64}()



    for qmesh in 1:Nq^2

        for jc in 1:Nq^2, jd in 1:Nq^2


            k1mesh=[
                mod(allowedq[jc][1]+allowedq[qmesh][1],Nq),
                mod(allowedq[jc][2]+allowedq[qmesh][2],Nq)
            ]

            k2mesh=[
                mod(allowedq[jd][1]-allowedq[qmesh][1],Nq),
                mod(allowedq[jd][2]-allowedq[qmesh][2],Nq)
            ]


            k1mesh_pos=get(allowedq_dict,Tuple(k1mesh),0)
            k2mesh_pos=get(allowedq_dict,Tuple(k2mesh),0)


            if k1mesh_pos==0 || k2mesh_pos==0
                continue
            end



            fourpoint=0.0+0.0im


            for je in eachindex(state_integer)

                (sign,state)=Cann(jc,state_integer[je],1)
                (sign,state)=Cann(jd,state,sign)
                (sign,state)=Cdag(k2mesh_pos,state,sign)
                (sign,state)=Cdag(k1mesh_pos,state,sign)


                if state*sign≠0

                    state_index=searchsortedfirst(state_integer,state)

                    if state_index<=length(state_integer) &&
                       state_integer[state_index]==state

                        fourpoint+=sign*
                                   conj(state_vector[state_index])*
                                   state_vector[je]

                    end

                end

            end


            abs(fourpoint)<1e-14 && continue



            for qg in eachindex(wave_diff)

                qtotal=allowedq[qmesh]+wave_diff[qg]
                qtotal_key=Tuple(qtotal)


                gk3pq_int=allowedq[jc]-k1mesh+qtotal
                gk4mq_int=allowedq[jd]-k2mesh-qtotal


                gk3pq_pos=get(
                    wave_diff_dict,
                    Tuple(gk3pq_int),
                    0
                )

                gk4mq_pos=get(
                    wave_diff_dict,
                    Tuple(gk4mq_int),
                    0
                )


                if gk3pq_pos≠0 && gk4mq_pos≠0

                    G_dic[qtotal_key]=
                        get(G_dic,qtotal_key,0.0+0.0im)+
                        Fmatrix[k1mesh_pos,jc,gk3pq_pos]*
                        Fmatrix[k2mesh_pos,jd,gk4mq_pos]*
                        fourpoint

                end

            end

        end

    end


    return G_dic

end


function get_rho(state_integer::Vector{Int64},
                 state_vector::Vector{ComplexF64},
                 Nq::Int64)


    rho_matrix=zeros(ComplexF64,Nq^2,Nq^2)


    for a in 1:Nq^2, b in 1:Nq^2

        rho_ab=0.0+0.0im


        for je in eachindex(state_integer)

            (sign,state)=Cann(a,state_integer[je],1)
            (sign,state)=Cdag(b,state,sign)


            if state*sign≠0

                state_index=searchsortedfirst(state_integer,state)

                if state_index<=length(state_integer) &&
                   state_integer[state_index]==state

                    rho_ab+=sign*
                            conj(state_vector[state_index])*
                            state_vector[je]

                end

            end

        end


        rho_matrix[a,b]=rho_ab

    end


    return rho_matrix

end