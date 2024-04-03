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


    wave_diff=Vector{Int64}[]
    cutoff=36
    cutoffstandard=12.01*scale
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        gtest=ja*b1+jb*b2;
        if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
            push!(wave_diff,ja*b1T+jb*b2T)
        end
    end
   
    
    
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

    
    return overlapmatrix, wave,wave_diff, single_MoirePo, single_Ham,allowedq, T1, T2
      

       
end


function overlap(k::Vector{Float64},q::Vector{Float64},β::Float64)::ComplexF64
    v=q[1]^2+q[2]^2+2*im*(k[1]*q[2]-k[2]*q[1])
    #v=2*im*(k[1]*q[2]-k[2]*q[1])
    return exp(-β/4*v)
end


function  get_foverlap(wave_diff,wave,allowedq,Nq,T1,T2)
    form_overlapmatrix=zeros(ComplexF64,Nq^2,Nq^2,length(wave),length(wave_diff))


    for ja in 1:Nq^2, jb in 1:Nq^2, jc in eachindex(wave), jd in eachindex(wave_diff)
        k1=allowedq[ja][1]*T1+allowedq[ja][2]*T2
        k2=allowedq[jb][1]*T1+allowedq[jb][2]*T2
        gprime=wave[jc][1]*T1+wave[jc][2]*T2
        gdiff=wave_diff[jd][1]*T1+wave_diff[jd][2]*T2
        form_overlapmatrix[ja,jb,jc,jd]=overlap(k1+gprime+gdiff,k2-k1-gdiff,β)
    end

   return form_overlapmatrix
end



function get_HFeigenvectors(loop_dic::Dict{Vector{Int},Any},allowedq::Vector{Vector{Int}},T1::Vector{Float64},T2::Vector{Float64},Nq::Int64,wave::Vector{Vector{Int64}},input_DensityMatrix::Vector{Matrix{ComplexF64}},single_Ham::Vector{Matrix{ComplexF64}},single_MoirePo::Vector{Matrix{ComplexF64}},constq::Float64,overlapmatrix::Array{ComplexF64,4})::Tuple{Float64,Vector{Matrix{ComplexF64}},Vector{Vector{Float64}}}
  
 
    dimension=length(wave)
    HartreeMatrix=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
    FockMatrix=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
   HF_eigenvalue=Vector{Any}(undef,Nq^2)
   HF_eigenvector=Vector{Any}(undef,Nq^2)
 
  
   for jk in 1:Nq^2
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
 
   
 
 
   for jk in 1:Nq^2
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
 
  Flevel=(sort(reduce(vcat,HF_eigenvalue))[Nq^2+1]+sort(reduce(vcat,HF_eigenvalue))[Nq^2])/2
   
  
 
  return  Flevel,HF_eigenvector,HF_eigenvalue
end


function get_indexset(Flevel,HF_eigenvalue,wave,wave_diff,num_bandbelow,num_bandup,allowedq)
        
 dimension=length(wave)    
 FLindex=0
 for ja in 1:dimension
   if HF_eigenvalue[1][ja]<Flevel
    FLindex+=1
   end
 end
 Bandvector=Vector{Int64}[]

 for ja in 1:Nq^2
   for jb in FLindex-num_bandbelow+1:FLindex
       push!(Bandvector,[ja,jb])
   end
   kpQ=[mod(allowedq[ja][1]+bigQ[1],Nq),mod(allowedq[ja][2]+bigQ[2],Nq)]
   kmQ=[mod(allowedq[ja][1]-bigQ[1],Nq),mod(allowedq[ja][2]-bigQ[2],Nq)]
   kpQpos=findfirst(item->item==kpQ,allowedq)
   kmQpos=findfirst(item->item==kmQ,allowedq)
    for jb in FLindex+1:FLindex+num_bandup
    push!(Bandvector,[kpQpos,jb])
    push!(Bandvector,[kmQpos,jb])
    end
  end
 Bandvector=sort(unique(Bandvector))


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

 Aindexset=Vector{Vector{Int64}}[]
 for kindex in 1:Nq^2
     kpQ=[mod(allowedq[kindex][1]+bigQ[1],Nq),mod(allowedq[kindex][2]+bigQ[2],Nq)]
     kpQpos=findfirst(item->item==kpQ,allowedq)
     for hband in FLindex-num_bandbelow+1:FLindex, pband in FLindex+1:FLindex+num_bandup
         pos1=findfirst(item->item==[kpQpos,pband],Bandvector)
         pos2=findfirst(item->item==[kindex,hband],Bandvector)
         push!(Aindexset,[[pos1,kpQpos,pband],[pos2,kindex,hband]])
     end
 end
 
 AmQindexset=Vector{Vector{Int64}}[]
 for kindex in 1:Nq^2
     kmQ=[mod(allowedq[kindex][1]-bigQ[1],Nq),mod(allowedq[kindex][2]-bigQ[2],Nq)]
     kmQpos=findfirst(item->item==kmQ,allowedq)
     for hband in FLindex-num_bandbelow+1:FLindex, pband in FLindex+1:FLindex+num_bandup
         pos1=findfirst(item->item==[kmQpos,pband],Bandvector)
         pos2=findfirst(item->item==[kindex,hband],Bandvector)
         push!(AmQindexset,[[pos1,kmQpos,pband],[pos2,kindex,hband]])
     end
 end


 B2indexset=Vector{Vector{Int64}}[]
 for kindex in 1:Nq^2
    kmQ=[mod(allowedq[kindex][1]-bigQ[1],Nq),mod(allowedq[kindex][2]-bigQ[2],Nq)]
    kmQpos=findfirst(item->item==kmQ,allowedq)
    for hband in FLindex-num_bandbelow+1:FLindex, pband in FLindex+1:FLindex+num_bandup
        pos1=findfirst(item->item==[kmQpos,pband],Bandvector)
        pos2=findfirst(item->item==[kindex,hband],Bandvector)
        push!(B2indexset,[[pos1,kmQpos,pband],[pos2,kindex,hband]])
    end
 end
end