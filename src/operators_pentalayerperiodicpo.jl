function single_particle_periodicpo(θ::Float64,Nq::Int64,uD::Float64,Nband::Int64,Vperiod::Float64,period::Float64)

    ac=0.246
    R1=ac*[1,0]
    R2=ac*[1/2,√3/2]
    G1=2π/ac*[1,-1/√3]
    G2=2π/ac*[0,2/√3]
    ϵ=0.2504/ac-1
    
    NL=5
 

    Rθ=[cos(θ) -sin(θ);sin(θ) cos(θ)]
    #g1=(G1-(1+ϵ)^(-1)*Rθ*G1)
    #g2=(G2-(1+ϵ)^(-1)*Rθ*G2)
    g1=G1/norm(G1)*4π/(period*√3)
    g2=G2/norm(G1)*4π/(period*√3)
    T1=g1/Nq
    T2=g2/Nq
    g1T=Int.(round.(inv([T1 T2])*g1))
    g2T=Int.(round.(inv([T1 T2])*g2))
    

    phaseangle=0.0
    
    
    
    g3=-g1-g2
    a1=inv([g1';g2'])*[2π,0]
    a2=inv([g1';g2'])*[0,2π]
    KGr=4π/(3*ac)*[1,0]
    am=norm(a1)

    
    allowedq=Vector{Int}[]
    for ja in 0:Nq-1, jb in 0:Nq-1
        push!(allowedq,[ja,jb])
    end

    allowedq_dic=Dict{Vector{Int},Int}()

  for ja in eachindex(allowedq)
     allowedq_dic[allowedq[ja]]=ja
   end
    
    
    wave=Vector{Int64}[]
    cutoff=18
    cutoffstandard=5.1*norm(g1)
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        gtest=ja*g1+jb*g2;
        if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
             push!(wave,ja*g1T+jb*g2T)
         end
    end
    
    

    eigenvalue=[zeros(Float64,Nband) for _ in 1:Nq^2]
    eigenvector=[zeros(ComplexF64,length(wave)*2*NL,Nband) for _ in 1:Nq^2]
    
    eigenvalue_forrecord=zeros(Float64,2*NL*length(wave),Nq^2)
    

   for ja in eachindex(allowedq)
       Ham=get_MoireHam_periodicpo([T1 T2]*allowedq[ja]+KGr,wave,NL,uD,T1,T2,Vperiod,phaseangle,g1T,g2T)
       FFF=eigen(Ham)
       eigenvalue[ja]=FFF.values[length(wave)*NL+1:length(wave)*NL+Nband]
       eigenvector[ja]=FFF.vectors[:,length(wave)*NL+1:length(wave)*NL+Nband]
       eigenvalue_forrecord[:,ja]=FFF.values
   end
   
  BW=real(sort(eigenvalue_forrecord[length(wave)*NL+1,:])[Nq^2]-sort(eigenvalue_forrecord[length(wave)*NL+1,:])[1])
  direct_gapup=sort(real.(eigenvalue_forrecord[length(wave)*NL+2,:]-eigenvalue_forrecord[length(wave)*NL+1,:]))[1]
  direct_gapdown=sort(real.(eigenvalue_forrecord[length(wave)*NL+1,:]-eigenvalue_forrecord[length(wave)*NL,:]))[1]
  indirect_gapup=sort(real.(eigenvalue_forrecord[length(wave)*NL+2,:]))[1]-sort(real.(eigenvalue_forrecord[length(wave)*NL+1,:]))[Nq^2]
  indirect_gapdown=sort(real.(eigenvalue_forrecord[length(wave)*NL+1,:]))[1]-sort(real.(eigenvalue_forrecord[length(wave)*NL,:]))[Nq^2]
 

   return eigenvector,eigenvalue,wave,allowedq,allowedq_dic,T1,T2,BW,direct_gapup,direct_gapdown,indirect_gapup, indirect_gapdown

end


function reshuffle(state::Vector{ComplexF64},NL,wave,shuffle_vec)
    permu_state=zeros(ComplexF64,2*NL*length(wave))
    for ja in eachindex(wave)
      pos=findfirst(item->item==wave[ja]-shuffle_vec,wave)
      if pos≠nothing
         permu_state[2*NL*(pos-1)+1:2*NL*pos]=state[2*NL*(ja-1)+1:2*NL*ja]
      end
    end

    return permu_state
end


function get_f(k::Vector{Float64})
    delta1=1/√3*0.246*[0,1]
    delta2=1/√3*0.246*[√3/2,-1/2]
    delta3=1/√3*0.246*[-√3/2,-1/2]
 
    return exp(im*dot(k,delta1))+exp(im*dot(k,delta2))+exp(im*dot(k,delta3))
end


function get_RNGham(k::Vector{Float64},NL::Int,uD::Float64)
    Ham=zeros(ComplexF64,2*NL,2*NL)
    t0=3100
    t1=380
    t2=-21
    t3=290
    t4=141
    for layer in 1:NL-1
       Ham[2*layer-1:2*layer,2*layer+1:2*layer+2]=[t4*get_f(k) t3*conj(get_f(k));t1 t4*get_f(k)]
    end

    for layer in 1:NL-2
        Ham[2*layer-1:2*layer,2*layer+3:2*layer+4]=[0.0 t2/2;0.0 0.0]
    end

    Ham=Ham+Ham'

    for layer in 1:NL
        Ham[2*layer-1:2*layer,2*layer-1:2*layer]=[uD*(layer+1-(NL-1)/2) -t0*get_f(k);-t0*conj(get_f(k)) uD*(layer+1-(NL-1)/2)]
    end
   return Ham
end



function get_MoireHam_periodicpo(k::Vector{Float64},wave::Vector{Vector{Int}},NL::Int,uD::Float64,T1::Vector{Float64},T2::Vector{Float64},Vperiod::Float64,angle::Float64,g1T::Vector{Int64},g2T::Vector{Int64})
    

    iden=Matrix{Float64}(I, 2*NL, 2*NL) 
      Hamiltonian=zeros(ComplexF64,length(wave)*2*NL,length(wave)*2*NL)
      for jb in eachindex(wave)
          kvec=k+[T1 T2]*wave[jb]
          Hamiltonian[2*NL*(jb-1)+1:2*NL*jb,2*NL*(jb-1)+1:2*NL*jb]=get_RNGham(kvec,NL,uD)
      end
     
      Moire=zeros(ComplexF64,length(wave)*2*NL,length(wave)*2*NL)
  
      for jc in eachindex(wave)
         pos=findfirst(item->item==wave[jc]+g1T,wave)
         if pos≠nothing
            Moire[2*NL*(pos-1)+1:2*NL*(pos),2*NL*(jc-1)+1:2*NL*(jc)]+=Vperiod*exp(-im*angle)*iden
         end
  
         pos=findfirst(item->item==wave[jc]+g2T,wave)
         if pos≠nothing
            Moire[2*NL*(pos-1)+1:2*NL*(pos),2*NL*(jc-1)+1:2*NL*(jc)]+=Vperiod*exp(-im*angle)*iden
  
         end
  
         pos=findfirst(item->item==wave[jc]-g1T-g2T,wave)
         if pos≠nothing
          Moire[2*NL*(pos-1)+1:2*NL*(pos),2*NL*(jc-1)+1:2*NL*(jc)]+=Vperiod*exp(-im*angle)*iden
         end
     
      end
  
  
 
      return Moire+Moire'+Hamiltonian
  
  
  
 
 end



 function get_chern(Nq::Int64,wave::Vector{Vector{Int64}},allowedq::Vector{Vector{Int}},allowedq_dic::Dict{Vector{Int},Int},T1::Vector{Float64},T2::Vector{Float64},eigenvector::Vector{Matrix{ComplexF64}})
    NL=5
    chern_eigenvector=zeros(ComplexF64,2*NL*length(wave),Nq+1,Nq+1)
    for ja in 1:Nq+1, jb in 1:Nq+1
        
        wavevec=mod.([ja-1,jb-1],Nq)
        if (ja<Nq+1) & (jb<Nq+1)
          chern_eigenvector[:,ja,jb]=eigenvector[allowedq_dic[wavevec]][:,1]
        end
    
        if (ja>Nq) || (jb>Nq)
            shuffle_vec=[ja-1,jb-1]-wavevec
    
           chern_eigenvector[:,ja,jb]=reshuffle(eigenvector[allowedq_dic[wavevec]][:,1],NL,wave,shuffle_vec)
        end
    end


    Uonelink=zeros(ComplexF64,Nq,Nq+1)
    Utwolink=zeros(ComplexF64,Nq+1,Nq)
    tra=zeros(Float64,Nq,Nq)
        
    for ja in 1:Nq, jb in 1:Nq+1
        Uonelink[ja,jb]=dot(chern_eigenvector[:,ja,jb],chern_eigenvector[:,ja+1,jb])/abs(dot(chern_eigenvector[:,ja,jb],chern_eigenvector[:,ja+1,jb]))
    end
    
        
        
    for ja in 1:Nq+1, jb in 1:Nq
        Utwolink[ja,jb]=dot(chern_eigenvector[:,ja,jb],chern_eigenvector[:,ja,jb+1])/abs(dot(chern_eigenvector[:,ja,jb],chern_eigenvector[:,ja,jb+1]))
    end
        
    Flink=zeros(ComplexF64,Nq,Nq)
    for ja in 1:Nq, jb in 1:Nq
     Flink[ja,jb]=log(Uonelink[ja,jb]*Utwolink[ja+1,jb]/(Uonelink[ja,jb+1]*Utwolink[ja,jb]))
    end
    chern=sum(Flink)/(2*π*im)
    
    aveF=sum(Flink)/Nq^2
    uniform=0.0
    for ja in 1:Nq, jb in 1:Nq
        uniform+=(imag(Flink[ja,jb])-imag(aveF))^2*Nq^2/(2π)^2
    end
    
    dG=norm(T2)
    
    for ja in 1:Nq, jb in 1:Nq
       
        A1=dot(chern_eigenvector[:,ja,jb],chern_eigenvector[:,ja,jb+1])
        B1=dot(chern_eigenvector[:,ja,jb],chern_eigenvector[:,ja+1,jb])
        C1=dot(chern_eigenvector[:,ja,jb],chern_eigenvector[:,ja+1,jb+1])
       gyy=(1-abs(A1)^2)/dG^2
       gxx=(2-1/2*gyy*dG^2-abs(B1)^2-abs(C1)^2)/(1.5*dG^2)
       tra[ja,jb]+=(gxx+gyy)*3^(1/2)/2*dG^2
    
      
    end
     
    trace_condition=sum(tra)-sum(abs.(Flink))

    return trace_condition,tra,Flink,chern,uniform

end