function get_Fmatrix(allowedq,HF_eigenvector,wave,wave_diff,form_overlapmatrix)
    dimension=length(wave)
    Fmatrix=zeros(ComplexF64,length(allowedq),length(allowedq),length(wave_diff))

    


    diff_eigenvector=zeros(ComplexF64,dimension,length(allowedq),length(wave_diff))
    for jc in eachindex(wave_diff), jd in eachindex(wave)
        pos=findfirst(item->item==wave[jd]+wave_diff[jc],wave)
        if pos≠nothing 
            diff_eigenvector[jd,:,jc]=HF_eigenvector[pos,:]
        end 
    end

    for ja in eachindex(allowedq), jb in eachindex(allowedq), jc in eachindex(wave_diff)
        Fmatrix[ja,jb,jc]=sum(conj(diff_eigenvector[:,ja,jc]).*HF_eigenvector[:,jb].*form_overlapmatrix[ja,jb,:,jc])
    end
  
    pos=findfirst(item->item==[0,0],wave_diff)
    for ja in eachindex(allowedq)
        Fmatrix[ja,ja,pos]=1.0+0.0*im
    end
    
  return Fmatrix
end


function get_Velement(qindex::Int64,v1::Int64,v2::Int64,v3::Int64,v4::Int64,Fmatrix::Array{ComplexF64},gkpqmap::Array{Int64},gkmqmap::Array{Int64},T1::Vector{Float64},T2::Vector{Float64},wave_diff::Vector{Vector{Int64}},constq::Float64)::ComplexF64
    Velement=0 
   for ja in eachindex(wave_diff)
        qvec=allowedq[qindex]+wave_diff[ja]
       if  (gkpqmap[v3,qindex,ja]≠0) && (gkmqmap[v4,qindex,ja]≠0) &&(qvec≠[0,0])
          Velement+=constq/norm(qvec[1]*T1+qvec[2]*T2)*Fmatrix[v1,v3,gkpqmap[v3,qindex,ja]]*Fmatrix[v2,v4,gkmqmap[v4,qindex,ja]]
       end
   end

   return Velement
end


function overlap(k::Vector{Float64},q::Vector{Float64},β::Float64)::ComplexF64
    v=q[1]^2+q[2]^2+2*im*(k[1]*q[2]-k[2]*q[1])
    #v=2*im*(k[1]*q[2]-k[2]*q[1])
    return exp(-β/4*v)
end


function single_Hamiltonian(flux::Float64,V0::Float64,ϕ::Float64,scale::Float64,Nq::Int64,wave::Vector{Vector{Int64}})
    am=4*π/(√3*scale);
    β=4*flux/(√3*scale^2)
    mass=0.5;
    
    b1=4*π/(√3*am)*[0,1]
    b2=4*π/(√3*am)*[√3/2,-1/2]


    
    T1=b1/(Nq)
    T2=b2/(Nq)
    
    
    
    b1T=Int.(round.(inv([T1 T2])*b1))
    b2T=Int.(round.(inv([T1 T2])*b2))
    
    
    
    allowedq=Vector{Int64}[]
    for ja in 0:Nq-1,jb in 0:Nq-1
        push!(allowedq,[ja,jb])
    end

    dimension=length(wave)
   
    
    
     
    single_Ham=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
    single_MoirePo=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
 
    
    
   for ja in 1:Nq^2
      
      
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
     
    
    end



    

    
    return single_MoirePo+single_Ham
      

       
end




function get_wavefunction(scale::Float64,flux::Float64,Nq::Int64,V0::Float64,ϕ::Float64,chispace::Vector{Float64},constq::Float64)
 
   
    am=4*π/(√3*scale);

    
    b1=4*π/(√3*am)*[0,1]
    b2=4*π/(√3*am)*[3^(1/2)/2,-1/2]

    T1=b1/(Nq)
    T2=b2/(Nq)
        
        
        
    b1T=Int.(round.(inv([T1 T2])*b1))
    b2T=Int.(round.(inv([T1 T2])*b2))
        
        
        
    allowedq=Vector{Int64}[]
    for ja in 0:Nq-1,jb in 0:Nq-1
            push!(allowedq,[ja,jb])
    end
    
        
        
        
    wave=Vector{Int64}[]
    cutoff=30
    cutoffstandard=8.01*scale
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        gtest=ja*b1+jb*b2;
        if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
                push!(wave,ja*b1T+jb*b2T)
        end
    end
    dimension=length(wave)
       
    
    wave_diff=Vector{Int64}[]
    cutoff=60
    cutoffstandard_diff=2*cutoffstandard
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        gtest=ja*b1+jb*b2;
        if (gtest[1]^2+gtest[2]^2)<cutoffstandard_diff^2
                push!(wave_diff,ja*b1T+jb*b2T)
        end
    end

   
    eigenvector=[zeros(ComplexF64,dimension,Nq^2) for _ in 1:length(chispace)];
    Abz=scale^2*√3/2

   for jchi in eachindex(chispace)
   for ja in 1:(Nq)^2,jc in 1:dimension
      n1=wave[jc][1]/Nq
      n2=wave[jc][2]/Nq
      kvec=[T1 T2]*allowedq[ja]
      gvec=[T1 T2]*wave[jc]
    
      phase=exp(im*π*((kvec[1]*gvec[2]-kvec[2]*gvec[1])/Abz+(n1-1)*(n2-1))) #I change flux to pi
      eigenvector[jchi][jc,ja]=exp(-1/(4*chispace[jchi]^2)*norm([T1 T2]*(allowedq[ja]+wave[jc]))^2)*conj(phase)
   end

   for ja in 1:Nq^2
       eigenvector[jchi][:,ja]=eigenvector[jchi][:,ja]/norm(eigenvector[jchi][:,ja])
   end
   end


   single_Ham=single_Hamiltonian(flux,V0,ϕ,scale,Nq,wave)


   return single_Ham, eigenvector, wave, wave_diff, allowedq, T1, T2





end


function get_formoverlap(Nq::Int64,allowedq::Vector{Vector{Int}},wave::Vector{Vector{Int}},wave_diff::Vector{Vector{Int}},flux::Float64, T1::Vector{Float64},T2::Vector{Float64})::Array{ComplexF64,4}
    form_overlapmatrix=zeros(ComplexF64,Nq^2,Nq^2,length(wave),length(wave_diff))
    β=4*flux/(√3*scale^2)
    form_overlapmatrix_parallel=[zeros(ComplexF64,Nq^2,Nq^2,length(wave)) for _ in 1:length(wave_diff)]

   Threads.@threads for jd in eachindex(wave_diff)
        for ja in 1:Nq^2, jb in 1:Nq^2, jc in eachindex(wave)
        k1=allowedq[ja][1]*T1+allowedq[ja][2]*T2
        k2=allowedq[jb][1]*T1+allowedq[jb][2]*T2
        gprime=wave[jc][1]*T1+wave[jc][2]*T2
        gdiff=wave_diff[jd][1]*T1+wave_diff[jd][2]*T2
        form_overlapmatrix_parallel[jd][ja,jb,jc]=overlap(k1+gprime+gdiff,k2-k1-gdiff,β)
        end
    end
    
   for jd in eachindex(wave_diff)
     form_overlapmatrix[:,:,:,jd]=form_overlapmatrix_parallel[jd]
   end
   form_overlapmatrix_parallel=nothing

    return form_overlapmatrix
end


function get_gkpgmap(Nq::Int64,allowedq::Vector{Vector{Int}},wave_diff::Vector{Vector{Int}})

    gkpqmap=zeros(Int64,length(allowedq),length(allowedq),length(wave_diff))
    gkmqmap=zeros(Int64,length(allowedq),length(allowedq),length(wave_diff))
    
    for ja in eachindex(allowedq),jb in eachindex(allowedq),jc in eachindex(wave_diff)
        kpq=allowedq[ja]+allowedq[jb]+wave_diff[jc]
        gkpq=kpq-[mod(kpq[1],Nq),mod(kpq[2],Nq)]
    
        kmq=allowedq[ja]-allowedq[jb]-wave_diff[jc]
        gkmq=kmq-[mod(kmq[1],Nq),mod(kmq[2],Nq)]
        
        gkpq_pos=findfirst(item->item==gkpq,wave_diff)
        gkmq_pos=findfirst(item->item==gkmq,wave_diff)
        if gkpq_pos≠nothing
            gkpqmap[ja,jb,jc]=gkpq_pos
        end
    
        if gkmq_pos≠nothing
            gkmqmap[ja,jb,jc]=gkmq_pos
        end
    end

   return gkpqmap,gkmqmap

end



function get_energy(chispace::Vector{Float64},constq::Float64,allowedq::Vector{Vector{Int}},form_overlapmatrix::Array{ComplexF64},single_Ham::Vector{Matrix{ComplexF64}},eigenvector::Vector{Matrix{ComplexF64}},wave::Vector{Vector{Int}},wave_diff::Vector{Vector{Int}},T1::Vector{Float64},T2::Vector{Float64},gkpqmap::Array{Int},gkmqmap::Array{Int})
   
   
    Energy=zeros(ComplexF64,length(chispace))
    kinetic=zeros(ComplexF64,length(chispace))
    Fock=zeros(ComplexF64,length(chispace))

    Threads.@threads for jchi in eachindex(chispace)
     Fmatrix=get_Fmatrix(allowedq,eigenvector[jchi],wave,wave_diff,form_overlapmatrix)
      for ja in 1:Nq^2
       Energy[jchi]+=eigenvector[jchi][:,ja]'*single_Ham[ja]*eigenvector[jchi][:,ja]
       kinetic[jchi]+=eigenvector[jchi][:,ja]'*single_Ham[ja]*eigenvector[jchi][:,ja]
      end
    
     for ja in 1:Nq^2,jb in 1:Nq^2
       if jb≠ja
         deltaq1=[0,0]
          q1index=findfirst(item->item==[mod(deltaq1[1],Nq),mod(deltaq1[2],Nq)],allowedq)
          deltaq2=allowedq[ja]-allowedq[jb]
          q2index=findfirst(item->item==[mod(deltaq2[1],Nq),mod(deltaq2[2],Nq)],allowedq)
    
          Energy[jchi]+=1/2*get_Velement(q1index,ja,jb,ja,jb,Fmatrix,gkpqmap,gkmqmap,T1,T2,wave_diff,constq)-1/2*get_Velement(q2index,ja,jb,jb,ja,Fmatrix,gkpqmap,gkmqmap,T1,T2,wave_diff,constq)
          Fock[jchi]+=-1/2*get_Velement(q2index,ja,jb,jb,ja,Fmatrix,gkpqmap,gkmqmap,T1,T2,wave_diff,constq)
       end
     end
    end

  return Energy,kinetic,Fock

end


function get_Fmatrix_noformoverlap(allowedq,HF_eigenvector,wave,wave_diff,β)
    dimension=length(wave)
    Fmatrix=zeros(ComplexF64,length(allowedq),length(allowedq),length(wave_diff))
    Fmatrix_threaded=[zeros(ComplexF64,length(allowedq),length(allowedq)) for _ in eachindex(wave_diff)]
    
    diff_eigenvector=zeros(ComplexF64,dimension,length(allowedq),length(wave_diff))
    for jc in eachindex(wave_diff), jd in eachindex(wave)
        pos=findfirst(item->item==wave[jd]+wave_diff[jc],wave)
        if pos≠nothing 
            diff_eigenvector[jd,:,jc]=HF_eigenvector[pos,:]
        end 
    end

   Threads.@threads for jc in eachindex(wave_diff)
      for ja in eachindex(allowedq), jb in eachindex(allowedq)

        k1=[T1 T2]*allowedq[ja]
        k2=[T1 T2]*allowedq[jb]
        gdiff=[T1 T2]*wave_diff[jc]
        
        form=zeros(ComplexF64,dimension)
        for jd in eachindex(wave)
            gprime=[T1 T2]*wave[jd]
            form[jd]=overlap(k1+gprime+gdiff,k2-k1-gdiff,β)
        end
        Fmatrix_threaded[jc][ja,jb]=sum(conj(diff_eigenvector[:,ja,jc]).*HF_eigenvector[:,jb].*form)
     end
   end

   for jc in eachindex(wave_diff)
     Fmatrix[:,:,jc]=Fmatrix_threaded[jc]
   end
  
    pos=findfirst(item->item==[0,0],wave_diff)
    for ja in eachindex(allowedq)
        Fmatrix[ja,ja,pos]=1.0+0.0*im
    end
    
  return Fmatrix
end


function get_energy_noformoverlap(chispace::Vector{Float64},constq::Float64,allowedq::Vector{Vector{Int}},single_Ham::Vector{Matrix{ComplexF64}},eigenvector::Vector{Matrix{ComplexF64}},wave::Vector{Vector{Int}},wave_diff::Vector{Vector{Int}},T1::Vector{Float64},T2::Vector{Float64},gkpqmap::Array{Int},gkmqmap::Array{Int},flux::Float64,scale::Float64)
   
    β=4*flux/(√3*scale^2)
    Energy=zeros(ComplexF64,length(chispace))
    kinetic=zeros(ComplexF64,length(chispace))
    Fock=zeros(ComplexF64,length(chispace))

    Threads.@threads for jchi in eachindex(chispace)
     Fmatrix=get_Fmatrix_noformoverlap(allowedq,eigenvector[jchi],wave,wave_diff,β)
      for ja in 1:Nq^2
       Energy[jchi]+=eigenvector[jchi][:,ja]'*single_Ham[ja]*eigenvector[jchi][:,ja]
       kinetic[jchi]+=eigenvector[jchi][:,ja]'*single_Ham[ja]*eigenvector[jchi][:,ja]
      end
    
     for ja in 1:Nq^2,jb in 1:Nq^2
       if jb≠ja
         deltaq1=[0,0]
          q1index=findfirst(item->item==[mod(deltaq1[1],Nq),mod(deltaq1[2],Nq)],allowedq)
          deltaq2=allowedq[ja]-allowedq[jb]
          q2index=findfirst(item->item==[mod(deltaq2[1],Nq),mod(deltaq2[2],Nq)],allowedq)
    
          Energy[jchi]+=1/2*get_Velement(q1index,ja,jb,ja,jb,Fmatrix,gkpqmap,gkmqmap,T1,T2,wave_diff,constq)-1/2*get_Velement(q2index,ja,jb,jb,ja,Fmatrix,gkpqmap,gkmqmap,T1,T2,wave_diff,constq)
          Fock[jchi]+=-1/2*get_Velement(q2index,ja,jb,jb,ja,Fmatrix,gkpqmap,gkmqmap,T1,T2,wave_diff,constq)
       end
     end
    end

  return Energy,kinetic,Fock

end