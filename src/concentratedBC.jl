
function overlap(k::Vector{Float64},q::Vector{Float64},spin::Float64,M::Float64)::ComplexF64
   
    v=(M^2+norm(k)^2+k[1]*q[1]+k[2]*q[2]-im*(k[1]*q[2]-k[2]*q[1]))^(Int(2*spin))/((M^2+norm(k)^2)^spin*(M^2+norm(k+q)^2)^spin)
    
    return v
end


function metric(wavelist::Vector{Vector{Int64}},spin::Float64,M::Float64,k::Vector{Float64},q::Vector{Float64},T1::Vector{Float64},T2::Vector{Float64})::Matrix{ComplexF64}
    Amatrix=zeros(ComplexF64,length(wavelist),length(wavelist))
    for ja in 1:length(wavelist)
    Amatrix[ja,ja]=overlap(k+wavelist[ja][1]*T1+wavelist[ja][2]*T2,q,spin,M)
    end
    return Amatrix
end



function get_MoireHam(k::Vector{Float64},mass::Float64,M::Float64,spin::Float64,wave::Vector{Vector{Int64}},T1::Vector{Float64},T2::Vector{Float64},V0::Float64,ϕ::Float64,g1T::Vector{Int64},g2T::Vector{Int64})
    dimension=length(wave)
    diagHam=zeros(ComplexF64,dimension,dimension)
    for ja in eachindex(wave)
       kvec=k+[T1 T2]*wave[ja]
       diagHam[ja,ja]=norm(kvec)^2/(2*mass)
    end
   
    MoirePo=zeros(ComplexF64,dimension,dimension)
    for ja in eachindex(wave)
    
       pos=findfirst(item->item==wave[ja]-g2T,wave)
       if pos≠nothing
           MoirePo[ja,pos]=V0*exp(im*ϕ)*overlap(k+[T1 T2]*wave[ja],[T1 T2]*g2T,spin,M)
       end
       
     
       pos=findfirst(item->item==wave[ja]-g1T,wave)
       if pos≠nothing
           MoirePo[ja,pos]=V0*exp(im*ϕ)*overlap(k+[T1 T2]*wave[ja],[T1 T2]*g1T,spin,M)
       end
   
       pos=findfirst(item->item==wave[ja]+g1T+g2T,wave)
       if pos≠nothing
           MoirePo[ja,pos]=V0*exp(im*ϕ)*overlap(k+[T1 T2]*wave[ja],-[T1 T2]*(g2T+g1T),spin,M)
       end
    end

    return diagHam+MoirePo+MoirePo'
 
end

function get_eigenvector(Nq::Int64,spin::Float64,scale::Float64,ϕ::Float64,V0::Float64,mass::Float64,M::Float64)


    g1=scale*[0,1]
    g2=scale*[√3/2,-1/2]
    T1=g1/Nq
    T2=g2/Nq
    
    chern_allowedq=Vector{Int}[]
    wave=Vector{Int}[]
    for ja in 0:Nq, jb in 0:Nq
     push!(chern_allowedq,[ja,jb])
    end
    
    g1T=Int.(round.(inv([T1 T2])*g1))
    g2T=Int.(round.(inv([T1 T2])*g2))
    
    
    
    cutoff=18
    cutoffstandard=5.01*scale
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        gtest=ja*g1+jb*g2;
        if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
            push!(wave,ja*g1T+jb*g2T)
        end
    end
    dimension=length(wave)
    
    
    
    chern_eigenvector=zeros(ComplexF64,dimension,Nq+1,Nq+1)
    eigenvalue=zeros(ComplexF64,dimension,(Nq+1)^2)
    chern_eigenvector_threads=Vector{Vector{ComplexF64}}(undef,(Nq+1)^2)
    eigenvalue_threads=Vector{Vector{ComplexF64}}(undef,(Nq+1)^2)
    Threads.@threads for ja in eachindex(chern_allowedq)
       k=[T1 T2]*chern_allowedq[ja]
       Hamiltonian=get_MoireHam(k,mass,M,spin,wave,T1,T2,V0,ϕ,g1T,g2T)
       FFF=eigen(Hamiltonian)
       #chern_eigenvector[:,chern_allowedq[ja][1]+1,chern_allowedq[ja][2]+1]=FFF.vectors[:,1]
       #eigenvalue[:,ja]=FFF.values
       chern_eigenvector_threads[ja]=FFF.vectors[:,1]
       eigenvalue_threads[ja]=FFF.values
    end


    for ja in eachindex(chern_allowedq)
        chern_eigenvector[:,chern_allowedq[ja][1]+1,chern_allowedq[ja][2]+1]=chern_eigenvector_threads[ja]
        eigenvalue[:,ja]=eigenvalue_threads[ja]
    end

   directgap=sort(real.(eigenvalue[2,:]-eigenvalue[1,:]))[1]
    BW=sort(real.(eigenvalue[1,:]))[(Nq+1)^2]-sort(real.(eigenvalue[1,:]))[1]
   indirectgap=sort(real.(eigenvalue[2,:]))[1]-sort(real.(eigenvalue[1,:]))[(Nq+1)^2]


   return directgap,BW,indirectgap,chern_eigenvector,wave,T1,T2
end


function get_chern(Nq::Int64,wave::Vector{Vector{Int64}},spin::Float64,M::Float64,T1::Vector{Float64},T2::Vector{Float64})

    Uonelink=zeros(ComplexF64,Nq,Nq+1)
    Utwolink=zeros(ComplexF64,Nq+1,Nq)
        
    for ja in 1:Nq, jb in 1:Nq+1
        Amatrix=metric(wave,spin,M,(ja-1)*T1+(jb-1)*T2,T1,T1,T2)
        Uonelink[ja,jb]=dot(chern_eigenvector[:,ja,jb],Amatrix*chern_eigenvector[:,ja+1,jb])/abs(dot(chern_eigenvector[:,ja,jb],Amatrix*chern_eigenvector[:,ja+1,jb]))
    end
    
        
        
    for ja in 1:Nq+1, jb in 1:Nq
        Amatrix=metric(wave,spin,M,(ja-1)*T1+(jb-1)*T2,T2,T1,T2)
        Utwolink[ja,jb]=dot(chern_eigenvector[:,ja,jb],Amatrix*chern_eigenvector[:,ja,jb+1])/abs(dot(chern_eigenvector[:,ja,jb],Amatrix*chern_eigenvector[:,ja,jb+1]))
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

    tra=zeros(ComplexF64,Nq,Nq)
    dG=norm(T2)
    for ja in 1:Nq, jb in 1:Nq
        Amatrix=metric(wave,spin,M,(ja-1)*T1+(jb-1)*T2,T1,T1,T2)
        Bmatrix=metric(wave,spin,M,(ja-1)*T1+(jb-1)*T2,T2,T1,T2)
        Cmatrix=metric(wave,spin,M,(ja-1)*T1+(jb-1)*T2,T2+T1,T1,T2)
        A1=dot(chern_eigenvector[:,ja,jb],Amatrix*chern_eigenvector[:,ja+1,jb])
        B1=dot(chern_eigenvector[:,ja,jb],Bmatrix*chern_eigenvector[:,ja,jb+1])
        C1=dot(chern_eigenvector[:,ja,jb],Cmatrix*chern_eigenvector[:,ja+1,jb+1])
       gyy=(1-abs(A1)^2)/dG^2
       gxx=(2-1/2*gyy*dG^2-abs(B1)^2-abs(C1)^2)/(1.5*dG^2)
       tra[ja,jb]+=(gxx+gyy)*3^(1/2)/2*dG^2
      
    end
        
  trace_condition=sum(tra)-sum(abs.(Flink))



    return chern,uniform,trace_condition
end