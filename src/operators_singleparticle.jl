using LinearAlgebra
using Arpack
using Combinatorics

using Random



function overlap(k::Vector{Float64},q::Vector{Float64},β::Float64)::ComplexF64
    v=q[1]^2+q[2]^2+2*im*(k[1]*q[2]-k[2]*q[1])
    #v=2*im*(k[1]*q[2]-k[2]*q[1])
    return exp(-β/4*v)
end
 
 
 
 



function square_initial_Densitymatrix(flux::Float64,V0::Float64,ϕ::Float64,scale::Float64,Nq::Int64,maxg::Float64)
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
    cutoff=3*floor(maxg)
    cutoffstandard=maxg*scale
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        gtest=ja*b1+jb*b2;
        if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
            push!(wave,ja*b1T+jb*b2T)
        end
    end
    dimension=length(wave)
   
    

    
     
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



     

    

    
    return  wave, single_MoirePo, single_Ham, single_eigenvalue,allowedq, T1, T2, a1m, a2m
      

       
end





function metric(wavelist::Vector{Vector{Int64}},β::Float64,k::Vector{Float64},q::Vector{Float64},T1::Vector{Float64},T2::Vector{Float64})::Matrix{ComplexF64}
    Amatrix=zeros(ComplexF64,length(wavelist),length(wavelist))
    for ja in 1:length(wavelist)
    Amatrix[ja,ja]=overlap(k+wavelist[ja][1]*T1+wavelist[ja][2]*T2,q,β)
    end
    return Amatrix
end









function square_chern(Nq::Int,wave::Vector{Vector{Int}},scale::Float64,ϕ::Float64,flux::Float64)

  
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
      eigenvector_intermediate_single[ja]=eigvecs(chern_Ham+chern_MoirePo)[:,1] 
    
    end


    eigenvector_bc_single=zeros(ComplexF64,dimension,Nq+1,Nq+1)


    for ja in eachindex(chern_allowedq)
      eigenvector_bc_single[:,chern_allowedq[ja][1]+1,chern_allowedq[ja][2]+1]=eigenvector_intermediate_single[ja]
    end
    
    

    tra_single=zeros(ComplexF64,Nq,Nq)

    

    
    dG=norm(T2)
    for ja in 1:Nq, jb in 1:Nq
        Amatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T1,T1,T2)
        Bmatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T2,T1,T2)
   
        A1=dot(eigenvector_bc_single[:,ja,jb],Amatrix*eigenvector_bc_single[:,ja+1,jb])
        B1=dot(eigenvector_bc_single[:,ja,jb],Bmatrix*eigenvector_bc_single[:,ja,jb+1])
      
       gyy=(1-abs(A1)^2)/dG^2
       gxx=(1-abs(B1)^2)/(dG^2)
       tra_single[ja,jb]+=(gxx+gyy)*dG^2
      
    end

   


    Uonelink=zeros(ComplexF64,Nq,Nq+1)
    Utwolink=zeros(ComplexF64,Nq+1,Nq)
  
    
    for ja in 1:Nq, jb in 1:Nq+1
        Amatrix=(metric(wave,β,(ja-1)*T1+(jb-1)*T2,T1,T1,T2))
        Uonelink[ja,jb]=dot(eigenvector_bc_single[:,ja,jb],Amatrix*eigenvector_bc_single[:,ja+1,jb])/abs(dot(eigenvector_bc_single[:,ja,jb],Amatrix*eigenvector_bc_single[:,ja+1,jb]))
    end
    
    for ja in 1:Nq+1, jb in 1:Nq
        Amatrix=(metric(wave,β,(ja-1)*T1+(jb-1)*T2,T2,T1,T2))
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
    




  trace_condition_single=sum(tra_single)-sum(abs.(Flink_single))
 

    
      
   return chern_single,Flink_single,trace_condition_single,uniform_single

end









