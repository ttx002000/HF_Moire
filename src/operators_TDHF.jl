function single_part(flux::Float64,V0::Float64,ϕ::Float64,scale::Float64,Nq::Int64)
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
    cutoffstandard=6.01*scale
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

    
    return overlapmatrix, wave, single_MoirePo, single_Ham,allowedq, T1, T2
      

       
end