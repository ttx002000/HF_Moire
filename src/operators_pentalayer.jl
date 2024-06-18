



function single_particle(ϵr::Float64,θ::Float64,Nq::Int64,uD::Float64,Nband::Int64)

    ac=0.246
    R1=ac*[1,0]
    R2=ac*[1/2,√3/2]
    G1=2π/ac*[1,-1/√3]
    G2=2π/ac*[0,2/√3]
    ϵ=0.2504/ac-1
    
    #ϵ=0.650313445592362/(norm(G1)-0.650313445592362)
    #ϵ=0.66/(norm(G1)-0.66)
    NL=5
    V0=28.9
    V1=21.0
    #V0=0.0
    #V1=0.0

    ψ=-0.29
    Rθ=[cos(θ) -sin(θ);sin(θ) cos(θ)]
    g1=(G1-(1+ϵ)^(-1)*Rθ*G1)*1/2
    g2=(G2-(1+ϵ)^(-1)*Rθ*G2)*1/2
    T1=g1/Nq
    T2=g2/Nq
    g1T=Int.(round.(inv([T1 T2])*g1))
    g2T=Int.(round.(inv([T1 T2])*g2))
    
    
    
    
    
    g3=-g1-g2
    a1=inv([g1';g2'])*[2π,0]
    a2=inv([g1';g2'])*[0,2π]
    KGr=4π/(3*ac)*[1,0]
    am=norm(a1)
    
    constq=1/(√3/2*am^2*ϵr*Nq^2)*9047.5636
    
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
    
    
    
    wave_diff=Vector{Int64}[]
    cutoff=18
    cutoffstandard_diff=5.1*norm(g1)
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        gtest=ja*g1+jb*g2;
        if (gtest[1]^2+gtest[2]^2)<cutoffstandard_diff^2
             push!(wave_diff,ja*g1T+jb*g2T)
         end
    end


    eigenvalue=[zeros(Float64,Nband) for _ in 1:Nq^2]
    eigenvector=[zeros(ComplexF64,length(wave)*2*NL,Nband) for _ in 1:Nq^2]



   for ja in eachindex(allowedq)
       Ham=get_MoireHam([T1 T2]*allowedq[ja]+KGr,wave,NL,uD,T1,T2,V1,V0,ψ,g1T,g2T)
       FFF=eigen(Ham)
       eigenvalue[ja]=FFF.values[length(wave)*NL+1:length(wave)*NL+Nband]
       eigenvector[ja]=FFF.vectors[:,length(wave)*NL+1:length(wave)*NL+Nband]
   end
   



   form_factors=Array{Matrix{ComplexF64}}(undef,Nq^2,Nq^2,length(wave_diff))
   #diff_eigenvector=zeros(ComplexF64,length(wave)*2*NL,Nband,Nq^2,length(wave_diff))
   diff_eigenvector_threads=[zeros(ComplexF64,length(wave)*2*NL,Nband,Nq^2) for _ in 1:length(wave_diff)]
   form_factors_threads=[Array{Matrix{ComplexF64}}(undef,Nq^2,length(wave_diff)) for _ in 1:Nq^2]

   #=
   for ja in eachindex(wave), jd in eachindex(wave_diff)
       pos=findfirst(item->item==wave[ja]-wave_diff[jd],wave)
       if pos≠nothing
          for jc in 1:Nq^2
                 diff_eigenvector[(ja-1)*2*NL+1:ja*2*NL,:,jc,jd]=eigenvector[jc][(pos-1)*2*NL+1:pos*2*NL,:]
          end
       end
   end
   
   =#
 Threads.@threads for jd in eachindex(wave_diff)
   for ja in eachindex(wave)
    pos=findfirst(item->item==wave[ja]-wave_diff[jd],wave)
    if pos≠nothing
       for jc in 1:Nq^2
              diff_eigenvector_threads[jd][(ja-1)*2*NL+1:ja*2*NL,:,jc]=eigenvector[jc][(pos-1)*2*NL+1:pos*2*NL,:]
       end
    end
  end
  end



   #=
   for ja in 1:Nq^2, jb in 1:Nq^2
     
     meshk1plusq=mod.(allowedq[ja]+allowedq[jb],Nq)
     k2pos=findfirst(item->item==meshk1plusq,allowedq)
     for  jc in eachindex(wave_diff)
        gk1plusq=wave_diff[jc]+allowedq[ja]+allowedq[jb]-meshk1plusq
        gk1plusq_pos=findfirst(item->item==gk1plusq,wave_diff)
   
        form_factors[ja,jb,jc]=zeros(ComplexF64,Nband,Nband)
        if gk1plusq_pos≠nothing
          for α in 1:Nband, β in 1:Nband
            form_factors[ja,jb,jc][α,β]=diff_eigenvector[:,α,ja,gk1plusq_pos]'*eigenvector[k2pos][:,β]
          end
        end
      end
   
   end
   =#

   Threads.@threads for ja in 1:Nq^2
    for jb in 1:Nq^2
     
        meshk1plusq=mod.(allowedq[ja]+allowedq[jb],Nq)
        k2pos=allowedq_dic[meshk1plusq]
       for  jc in eachindex(wave_diff)
           gk1plusq=wave_diff[jc]+allowedq[ja]+allowedq[jb]-meshk1plusq
           gk1plusq_pos=findfirst(item->item==gk1plusq,wave_diff)
           form_factors_threads[ja][jb,jc]=zeros(ComplexF64,Nband,Nband)
           if gk1plusq_pos≠nothing
               form_factors_threads[ja][jb,jc]=(diff_eigenvector_threads[gk1plusq_pos][:,:,ja])'*eigenvector[k2pos][:,:]
           end
       end
     end
   end

   for ja in 1:Nq^2,jb in 1:Nq^2, jc in eachindex(wave_diff)
      form_factors[ja,jb,jc]=form_factors_threads[ja][jb,jc]
   end

   return eigenvector,eigenvalue,wave,wave_diff,allowedq,allowedq_dic,T1,T2,form_factors,constq,ϵ

end



function get_initial_projector(Nq::Int64,Nband::Int64,eigenvalue::Vector{Vector{Float64}})

    band_Ham=[zeros(ComplexF64,Nband,Nband) for _ in 1:Nq^2]
    band_eigenvector=[zeros(ComplexF64,Nband,Nband) for _ in 1:Nq^2]
    for ja in 1:Nq^2
      band_Ham[ja]=diagm(eigenvalue[ja])
      FFF=eigen(band_Ham[ja])
      band_eigenvector=FFF.vectors
    end
    
    initial_projector=[zeros(ComplexF64,Nband,Nband) for _ in 1:Nq^2]
    for ja in 1:Nq^2
      initial_projector[ja]=band_eigenvector[:,1]*band_eigenvector[:,1]'
    end
    
    for ja in 1:Nq^2
      A=randn(Nband,Nband)+im*randn(Nband,Nband)
      initial_projector[ja]+=(A+A')*1.0
    end
   
  return initial_projector,band_Ham
end





function iteration(Nq::Int64,Nband::Int64,initial_projector::Vector{Matrix{ComplexF64}},form_factors::Array{Matrix{ComplexF64}},constq::Float64,wave_diff::Vector{Vector{Int}},allowedq::Vector{Vector{Int}},allowedq_dic::Dict{Vector{Int},Int},T1::Vector{Float64},T2::Vector{Float64},band_Ham::Vector{Matrix{ComplexF64}})
    eout=1.0
    itcount=0
    bad_count=0
    energy=0.0
    energy_change=0.0
    bound=0.0
    HF_eigenvalue=Vector{Any}(undef,Nq^2)
    HF_eigenvector=[zeros(ComplexF64,Nband,Nband) for _ in 1:Nq^2]
    DIIS_input_projector=Vector{Vector{Matrix{ComplexF64}}}(undef,3)
    DIIS_input_DeltaMatrix=Vector{Vector{Matrix{ComplexF64}}}(undef,3)
    input_projector=initial_projector
    
    while (eout>1*10^-15) || (bad_count<4) || (abs(energy_change)>1*10^-10)
      if eout<1*10^-15
       bad_count+=1
      end
      tic=time()
      eout,energy_change,output_projector,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue,bound,HF_eigenvector,energy=Construct_projector(form_factors,input_projector,constq,Nq,Nband,wave_diff,allowedq,allowedq_dic,T1,T2,band_Ham,energy)
      DIIS_input_projector[mod(itcount,3)+1]=input_projector
      input_projector=output_projector
      itcount+=1
     
      toc=time()
      println(toc-tic,"eout=$eout","energy_change=$energy_change")
      flush(stdout)
     
    
    end
    

   #=
    println("startDIIS",itcount)
    
    bad_count=0
    while (eout>1*10^-15) || (bad_count<4) || (abs(energy_change)>1*10^-10)
        if  eout<1*10^-15 
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
       
        dmk=coeff[1]*(DIIS_input_projector[1]+DIIS_input_DeltaMatrix[1])+coeff[2]*(DIIS_input_projector[2]+DIIS_input_DeltaMatrix[2])+coeff[3]*(DIIS_input_projector[3]+DIIS_input_DeltaMatrix[3])
        eout,energy_change,output_projector,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue,bound,HF_eigenvector,energy=Construct_projector(form_factors,dmk,constq,Nq,Nband,wave_diff,allowedq,allowedq_dic,T1,T2,band_Ham,energy)
        DIIS_input_projector[mod(itcount,3)+1]=dmk
        itcount+=1

        toc=time()
        println(toc-tic,"eout=$eout","energy_change=$energy_change")
        flush(stdout)
    end
  =#
   





    return DIIS_input_projector,energy,HF_eigenvalue,HF_eigenvector,bound
end








function get_f(k::Vector{Float64})
    delta1=1/√3*0.246*[0,1]
    delta2=1/√3*0.246*[√3/2,-1/2]
    delta3=1/√3*0.246*[-√3/2,-1/2]
 
    return exp(im*dot(k,delta1))+exp(im*dot(k,delta2))+exp(im*dot(k,delta3))
end


function Coulomb(k::Vector{Int},T1::Vector{Float64},T2::Vector{Float64})::Float64
   D=25
   return k==[0,0] ? D : tanh(norm([T1 T2]*k*D))/norm(k[1]*T1+k[2]*T2)
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


function get_MoireHam(k::Vector{Float64},wave::Vector{Vector{Int}},NL::Int,uD::Float64,T1::Vector{Float64},T2::Vector{Float64},V1::Float64,V0::Float64,ψ::Float64,g1T::Vector{Int},g2T::Vector{Int})
    
      ω=exp(im*2π/3)
        Hamiltonian=zeros(ComplexF64,length(wave)*2*NL,length(wave)*2*NL)
        for jb in eachindex(wave)
            kvec=k+[T1 T2]*wave[jb]
            Hamiltonian[2*NL*(jb-1)+1:2*NL*jb,2*NL*(jb-1)+1:2*NL*jb]=get_RNGham(kvec,NL,uD)
        end
       
        Moire=zeros(ComplexF64,length(wave)*2*NL,length(wave)*2*NL)
    
        for jc in eachindex(wave)
           pos=findfirst(item->item==wave[jc]+2*g1T,wave)
           if pos≠nothing
              Moire[2*NL*(pos-1)+1:2*NL*(pos-1)+2,2*NL*(jc-1)+1:2*NL*(jc-1)+2]+=V1*exp(-im*ψ)*[1 1;ω ω]
           end
    
           pos=findfirst(item->item==wave[jc]+2*g2T,wave)
           if pos≠nothing
              Moire[2*NL*(pos-1)+1:2*NL*(pos-1)+2,2*NL*(jc-1)+1:2*NL*(jc-1)+2]+=V1*exp(-im*ψ)*[1 ω^2;ω^2 ω]
    
           end
    
           pos=findfirst(item->item==wave[jc]-2*g1T-2*g2T,wave)
           if pos≠nothing
            Moire[2*NL*(pos-1)+1:2*NL*(pos-1)+2,2*NL*(jc-1)+1:2*NL*(jc-1)+2]+=V1*exp(-im*ψ)*[1 ω;1 ω]
           end
       
           Moire[2*NL*(jc-1)+1:2*NL*(jc-1)+2,2*NL*(jc-1)+1:2*NL*(jc-1)+2]+=V0/2*[1 0; 0 1]
        end
    
    

        return Moire+Moire'+Hamiltonian
    
    
    

end



function get_MoireHam_periodicpo(k::Vector{Float64},wave::Vector{Vector{Int}},NL::Int,uD::Float64,T1::Vector{Float64},T2::Vector{Float64},Vperiod::Float64,angle::Float64)
    

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



function Construct_projector(form_factor::Array{Matrix{ComplexF64}},projector::Vector{Matrix{ComplexF64}},constq::Float64,Nq::Int64,Nband::Int64,wave_diff::Vector{Vector{Int64}},allowedq::Vector{Vector{Int}},allowedq_dic::Dict{Vector{Int},Int},T1::Vector{Float64},T2::Vector{Float64},band_Ham::Vector{Matrix{ComplexF64}},energy_input::Float64)
   Hartree=[zeros(ComplexF64,Nband,Nband) for _ in 1:Nq^2]
   Fock=[zeros(ComplexF64,Nband,Nband) for _ in 1:Nq^2]
   New_projector=[zeros(ComplexF64,Nband,Nband) for _ in 1:Nq^2]
   
   Threads.@threads for k2 in eachindex(allowedq)
    for jqmesh in eachindex(allowedq)
       k1pos=allowedq_dic[mod.(allowedq[k2]+allowedq[jqmesh],Nq)]
      for jqg in eachindex(wave_diff)
         Fock[k2]+=Coulomb(allowedq[jqmesh]+wave_diff[jqg],T1,T2)*form_factor[k2,jqmesh,jqg]*projector[k1pos]*(form_factor[k2,jqmesh,jqg])'
       end
   end 
   end
  
  HartreeDensity=zeros(ComplexF64,length(wave_diff))
  zeropos=allowedq_dic[[0,0]]
  for jqg in eachindex(wave_diff), jb in 1:Nq^2
     HartreeDensity[jqg]+=tr(projector[jb]*form_factor[jb,zeropos,jqg])
  end

  Threads.@threads for ja in eachindex(allowedq)
    for jqg in eachindex(wave_diff)
       Hartree[ja]+=Coulomb(wave_diff[jqg],T1,T2)*form_factor[ja,zeropos,jqg]*HartreeDensity[jqg]
    end
  end

  HF_eigenvalue=[zeros(Float64,Nband) for _ in 1:Nq^2]
  HF_eigenvector=[zeros(ComplexF64,Nband,Nband) for _ in 1:Nq^2]
  for ja in eachindex(allowedq)
     FFF=eigen(constq*Hartree[ja]+band_Ham[ja]-constq*Fock[ja])
     HF_eigenvalue[ja]=real.(FFF.values)
     HF_eigenvector[ja]=FFF.vectors
  end



  sorted=sort(reduce(vcat,HF_eigenvalue))
  bound=(sorted[Nq^2+1]+sorted[Nq^2])/2
 
  Threads.@threads for ja in 1:Nq^2
        for jd in eachindex(HF_eigenvalue[ja])
           if HF_eigenvalue[ja][jd]<bound
              New_projector[ja]+=HF_eigenvector[ja][:,jd]*(HF_eigenvector[ja][:,jd])'
           end
        end
   end

   DeltaMatrix=New_projector-projector
   output_projector=0.4*projector+0.6*New_projector

   e1=0.0
   for ja in 1:Nq^2
     e1+=tr(DeltaMatrix[ja]'*DeltaMatrix[ja])
   end
   eout=real(e1)/Nq^2
    
  energy=0.0
  for ja in 1:Nq^2
     energy+=tr((constq/2*Hartree[ja]+band_Ham[ja]-constq/2*Fock[ja])*output_projector[ja])
  end


  energy_change=real(energy-energy_input)

   return  eout,energy_change,output_projector,DeltaMatrix,HF_eigenvalue,bound, HF_eigenvector,real(energy)
  

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


function get_chern(Nq::Int64,wave::Vector{Vector{Int64}},allowedq::Vector{Vector{Int}},allowedq_dic::Dict{Vector{Int},Int},T1::Vector{Float64},T2::Vector{Float64},eigenvector::Vector{Matrix{ComplexF64}},HF_eigenvector::Vector{Matrix{ComplexF64}})
    NL=5
    chern_eigenvector=zeros(ComplexF64,2*NL*length(wave),Nq+1,Nq+1)
    for ja in 1:Nq+1, jb in 1:Nq+1
        
        wavevec=mod.([ja-1,jb-1],Nq)
        if (ja<Nq+1) & (jb<Nq+1)
          chern_eigenvector[:,ja,jb]=eigenvector[allowedq_dic[wavevec]]*HF_eigenvector[allowedq_dic[wavevec]][:,1]
        end
    
        if (ja>Nq) || (jb>Nq)
            shuffle_vec=[ja-1,jb-1]-wavevec
    
           chern_eigenvector[:,ja,jb]=reshuffle(eigenvector[allowedq_dic[wavevec]]*HF_eigenvector[allowedq_dic[wavevec]][:,1],NL,wave,shuffle_vec)
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