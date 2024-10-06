function single_part(flux::Float64,V0::Float64,ϕ::Float64,scale::Float64,Nq::Int64)
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

    
    
    
    wave=Vector{Int64}[]
    cutoff=18
    cutoffstandard=4.01*scale
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        gtest=ja*b1+jb*b2;
        if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
            push!(wave,ja*b1T+jb*b2T)
        end
    end
    dimension=length(wave)


    wave_diff=Vector{Int64}[]
    cutoff=36
    cutoffstandard=8.01*scale
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


function  get_foverlap(wave_diff::Vector{Vector{Int64}},wave::Vector{Vector{Int64}},allowedq::Vector{Vector{Int64}},Nq::Int64,T1::Vector{Float64},T2::Vector{Float64},flux::Float64)::Array{ComplexF64}
    form_overlapmatrix=zeros(ComplexF64,Nq^2,Nq^2,length(wave),length(wave_diff))
    scale=1.0
    β=4*flux/(√3*scale^2)
    for ja in 1:Nq^2, jb in 1:Nq^2, jc in eachindex(wave), jd in eachindex(wave_diff)
        k1=allowedq[ja][1]*T1+allowedq[ja][2]*T2
        k2=allowedq[jb][1]*T1+allowedq[jb][2]*T2
        gprime=wave[jc][1]*T1+wave[jc][2]*T2
        gdiff=wave_diff[jd][1]*T1+wave_diff[jd][2]*T2
        form_overlapmatrix[ja,jb,jc,jd]=overlap(k1+gprime+gdiff,k2-k1-gdiff,β)
    end

   return form_overlapmatrix
end


function Coulomb(k::Vector{Int64},T1::Vector{Float64},T2::Vector{Float64})::Float64
 

    return k==[0,0] ? 0.0 : 1/norm(k[1]*T1+k[2]*T2)
    
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


function get_indexset(Flevel::Float64,HF_eigenvalue::Vector{Vector{Float64}},wave::Vector{Vector{Int64}},wave_diff::Vector{Vector{Int64}},num_bandbelow::Int64,num_bandup::Int64,allowedq::Vector{Vector{Int64}},Nq::Int64,bigQ::Vector{Int64})::Tuple{Vector{Vector{Int64}},Vector{Vector{Vector{Int64}}},Vector{Vector{Vector{Int64}}},Vector{Vector{Vector{Int64}}},Array{Int64},Array{Int64}}
        
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
       push!(Bandvector,[ja,jb]) #I think this is the hole band
   end
      kpQ=[mod(allowedq[ja][1]+bigQ[1],Nq),mod(allowedq[ja][2]+bigQ[2],Nq)]
      kmQ=[mod(allowedq[ja][1]-bigQ[1],Nq),mod(allowedq[ja][2]-bigQ[2],Nq)]
      kpQpos=findfirst(item->item==kpQ,allowedq)
      kmQpos=findfirst(item->item==kmQ,allowedq)
    for jb in FLindex+1:FLindex+num_bandup
      push!(Bandvector,[kpQpos,jb])#I think this is the electron band
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

 return Bandvector,Aindexset,AmQindexset,B2indexset,gkpqmap,gkmqmap
end


function get_Fmatrix(Bandvector::Vector{Vector{Int64}},HF_eigenvector::Vector{Matrix{ComplexF64}},wave::Vector{Vector{Int}},wave_diff::Vector{Vector{Int}},form_overlapmatrix::Array{ComplexF64})::Array{ComplexF64}
    dimension=length(wave)
    Fmatrix=zeros(ComplexF64,length(Bandvector),length(Bandvector),length(wave_diff))
    eigenvector_single=zeros(ComplexF64,length(wave),length(Bandvector))

    for ja in eachindex(Bandvector)
        eigenvector_single[:,ja]=HF_eigenvector[Bandvector[ja][1]][:,Bandvector[ja][2]]
    end       


    diff_eigenvector=zeros(ComplexF64,dimension,length(Bandvector),length(wave_diff))
    for jc in eachindex(wave_diff), jd in eachindex(wave)
        pos=findfirst(item->item==wave[jd]+wave_diff[jc],wave)
        if pos≠nothing 
            diff_eigenvector[jd,:,jc]=eigenvector_single[pos,:]
        end 
    end

    for ja in eachindex(Bandvector), jb in eachindex(Bandvector), jc in eachindex(wave_diff)
        Fmatrix[ja,jb,jc]=sum(conj(diff_eigenvector[:,ja,jc]).*eigenvector_single[:,jb].*form_overlapmatrix[Bandvector[ja][1],Bandvector[jb][1],:,jc])
    end
  
    pos=findfirst(item->item==[0,0],wave_diff)
    for ja in eachindex(Bandvector)
        Fmatrix[ja,ja,pos]=1.0+0.0*im
    end
    
  return Fmatrix
end


function get_Velement(qindex::Int64,v1::Vector{Int64},v2::Vector{Int64},v3::Vector{Int64},v4::Vector{Int64},Fmatrix::Array{ComplexF64},gkpqmap::Array{Int64},gkmqmap::Array{Int64},T1::Vector{Float64},T2::Vector{Float64},wave_diff::Vector{Vector{Int64}},allowedq::Vector{Vector{Int64}},constq::Float64)::ComplexF64
    Velement=0 #v1v2v3v4 contains two informatin, which band and which momentum
   for ja in eachindex(wave_diff)
        qvec=allowedq[qindex]+wave_diff[ja]
       if  (gkpqmap[v3[2],qindex,ja]≠0) && (gkmqmap[v4[2],qindex,ja]≠0) &&(qvec≠[0,0])
          Velement+=constq/norm(qvec[1]*T1+qvec[2]*T2)*Fmatrix[v1[1],v3[1],gkpqmap[v3[2],qindex,ja]]*Fmatrix[v2[1],v4[1],gkmqmap[v4[2],qindex,ja]]
       end
   end

   return Velement
end



function test_function()
    Amatrixvector=[zeros(Float64,5) for _ in 1:5]
    Threads.@threads for ja in 1:5
        for jb in 1:5
       Amatrixvector[ja][jb]=ja-jb
    end
    end
    Amatrix=reduce(hcat,Amatrixvector)
    return Amatrix,Amatrixvector
end


function Construct_Amatrix(Aindexset::Vector{Vector{Vector{Int64}}},AmQindexset::Vector{Vector{Vector{Int64}}},B2indexset::Vector{Vector{Vector{Int64}}},allowedq::Vector{Vector{Int64}},Fmatrix::Array{ComplexF64},gkpqmap::Array{Int64},gkmqmap::Array{Int64},T1::Vector{Float64},T2::Vector{Float64},wave_diff::Vector{Vector{Int64}},HF_eigenvalue::Vector{Vector{Float64}},Nq::Int64,constq::Float64)::Tuple{Matrix{ComplexF64},Matrix{ComplexF64},Matrix{ComplexF64}}
    #=
    Amatrix=zeros(ComplexF64,length(Aindexset),length(Aindexset))
 

    for ja in eachindex(Aindexset),jb in eachindex(Aindexset)
   v1=Aindexset[ja][1]
   v2=Aindexset[jb][2]
   v3=Aindexset[jb][1]
   v4=Aindexset[ja][2]
   deltaq1=allowedq[v1[2]]-allowedq[v3[2]]
   q1index=findfirst(item->item==[mod(deltaq1[1],Nq),mod(deltaq1[2],Nq)],allowedq)

   deltaq2=allowedq[v1[2]]-allowedq[v4[2]]
   q2index=findfirst(item->item==[mod(deltaq2[1],Nq),mod(deltaq2[2],Nq)],allowedq)
   Amatrix[ja,jb]=-get_Velement(q1index,v1[1:2],v2[1:2],v3[1:2],v4[1:2],Fmatrix,gkpqmap,gkmqmap,T1,T2,wave_diff)+get_Velement(q2index,v1[1:2],v2[1:2],v4[1:2],v3[1:2],Fmatrix,gkpqmap,gkmqmap,T1,T2,wave_diff)
  
    end



 for ja in eachindex(Aindexset)
   v1=Aindexset[ja][1]
   v2=Aindexset[ja][2]
   Amatrix[ja,ja]+=HF_eigenvalue[v1[2]][v1[3]]-HF_eigenvalue[v2[2]][v2[3]]
 end


 
 AmQmatrix=zeros(ComplexF64,length(AmQindexset),length(AmQindexset))


 for ja in eachindex(AmQindexset),jb in eachindex(AmQindexset)
   v1=AmQindexset[ja][1]
   v2=AmQindexset[jb][2]
   v3=AmQindexset[jb][1]
   v4=AmQindexset[ja][2]
   deltaq1=allowedq[v1[2]]-allowedq[v3[2]]
   q1index=findfirst(item->item==[mod(deltaq1[1],Nq),mod(deltaq1[2],Nq)],allowedq)

   deltaq2=allowedq[v1[2]]-allowedq[v4[2]]
   q2index=findfirst(item->item==[mod(deltaq2[1],Nq),mod(deltaq2[2],Nq)],allowedq)
   AmQmatrix[ja,jb]=-get_Velement(q1index,v1[1:2],v2[1:2],v3[1:2],v4[1:2],Fmatrix,gkpqmap,gkmqmap,T1,T2,wave_diff)+get_Velement(q2index,v1[1:2],v2[1:2],v4[1:2],v3[1:2],Fmatrix,gkpqmap,gkmqmap,T1,T2,wave_diff)
 end
 





  for ja in eachindex(AmQindexset)
   v1=AmQindexset[ja][1]
   v2=AmQindexset[ja][2]
   AmQmatrix[ja,ja]+=HF_eigenvalue[v1[2]][v1[3]]-HF_eigenvalue[v2[2]][v2[3]]
 end




  Bmatrix=zeros(ComplexF64,length(Aindexset),length(B2indexset))

   for ja in eachindex(Aindexset),jb in eachindex(B2indexset)
     v1=Aindexset[ja][1]
     v2=B2indexset[jb][1]
     v3=B2indexset[jb][2]
     v4=Aindexset[ja][2]
     deltaq1=allowedq[v1[2]]-allowedq[v3[2]]
     q1index=findfirst(item->item==[mod(deltaq1[1],Nq),mod(deltaq1[2],Nq)],allowedq)
  
     deltaq2=allowedq[v1[2]]-allowedq[v4[2]]
     q2index=findfirst(item->item==[mod(deltaq2[1],Nq),mod(deltaq2[2],Nq)],allowedq)
     Bmatrix[ja,jb]=-get_Velement(q1index,v1[1:2],v2[1:2],v3[1:2],v4[1:2],Fmatrix,gkpqmap,gkmqmap,T1,T2,wave_diff)+get_Velement(q2index,v1[1:2],v2[1:2],v4[1:2],v3[1:2],Fmatrix,gkpqmap,gkmqmap,T1,T2,wave_diff)
   end


  =#


   
    Amatrix=zeros(ComplexF64,length(Aindexset),length(Aindexset))
    Amatrixvec=[zeros(ComplexF64,length(Aindexset)) for _ in eachindex(Aindexset)]

 Threads.@threads for ja in eachindex(Aindexset)
  for  jb in eachindex(Aindexset)
   v1=Aindexset[ja][1]
   v2=Aindexset[jb][2]
   v3=Aindexset[jb][1]
   v4=Aindexset[ja][2]
   deltaq1=allowedq[v1[2]]-allowedq[v3[2]]
   q1index=findfirst(item->item==[mod(deltaq1[1],Nq),mod(deltaq1[2],Nq)],allowedq)

   deltaq2=allowedq[v1[2]]-allowedq[v4[2]]
   q2index=findfirst(item->item==[mod(deltaq2[1],Nq),mod(deltaq2[2],Nq)],allowedq)
   Amatrixvec[ja][jb]=-get_Velement(q1index,v1[1:2],v2[1:2],v3[1:2],v4[1:2],Fmatrix,gkpqmap,gkmqmap,T1,T2,wave_diff,allowedq,constq)+get_Velement(q2index,v1[1:2],v2[1:2],v4[1:2],v3[1:2],Fmatrix,gkpqmap,gkmqmap,T1,T2,wave_diff,allowedq,constq)
  end
 end

 for ja in eachindex(Aindexset)
    Amatrix[ja,:]=Amatrixvec[ja]
 end

 for ja in eachindex(Aindexset)
   v1=Aindexset[ja][1]
   v2=Aindexset[ja][2]
   Amatrix[ja,ja]+=HF_eigenvalue[v1[2]][v1[3]]-HF_eigenvalue[v2[2]][v2[3]]
 end


 
 AmQmatrix=zeros(ComplexF64,length(AmQindexset),length(AmQindexset))
 AmQmatrixvec=[zeros(ComplexF64,length(AmQindexset)) for _ in eachindex(AmQindexset)]

 Threads.@threads for ja in eachindex(AmQindexset)
 for jb in eachindex(AmQindexset)
   v1=AmQindexset[ja][1]
   v2=AmQindexset[jb][2]
   v3=AmQindexset[jb][1]
   v4=AmQindexset[ja][2]
   deltaq1=allowedq[v1[2]]-allowedq[v3[2]]
   q1index=findfirst(item->item==[mod(deltaq1[1],Nq),mod(deltaq1[2],Nq)],allowedq)

   deltaq2=allowedq[v1[2]]-allowedq[v4[2]]
   q2index=findfirst(item->item==[mod(deltaq2[1],Nq),mod(deltaq2[2],Nq)],allowedq)
   AmQmatrixvec[ja][jb]=-get_Velement(q1index,v1[1:2],v2[1:2],v3[1:2],v4[1:2],Fmatrix,gkpqmap,gkmqmap,T1,T2,wave_diff,allowedq,constq)+get_Velement(q2index,v1[1:2],v2[1:2],v4[1:2],v3[1:2],Fmatrix,gkpqmap,gkmqmap,T1,T2,wave_diff,allowedq,constq)
 end
 end

 for ja in eachindex(AmQindexset)
    AmQmatrix[ja,:]=AmQmatrixvec[ja]
 end



  for ja in eachindex(AmQindexset)
   v1=AmQindexset[ja][1]
   v2=AmQindexset[ja][2]
   AmQmatrix[ja,ja]+=HF_eigenvalue[v1[2]][v1[3]]-HF_eigenvalue[v2[2]][v2[3]]
 end




  Bmatrix=zeros(ComplexF64,length(Aindexset),length(B2indexset))
  Bmatrixvec=[zeros(ComplexF64,length(B2indexset)) for _ in eachindex(Aindexset)]

  Threads.@threads for ja in eachindex(Aindexset)
   for jb in eachindex(B2indexset)
     v1=Aindexset[ja][1]
     v2=B2indexset[jb][1]
     v3=B2indexset[jb][2]
     v4=Aindexset[ja][2]
     deltaq1=allowedq[v1[2]]-allowedq[v3[2]]
     q1index=findfirst(item->item==[mod(deltaq1[1],Nq),mod(deltaq1[2],Nq)],allowedq)
  
     deltaq2=allowedq[v1[2]]-allowedq[v4[2]]
     q2index=findfirst(item->item==[mod(deltaq2[1],Nq),mod(deltaq2[2],Nq)],allowedq)
     Bmatrixvec[ja][jb]=-get_Velement(q1index,v1[1:2],v2[1:2],v3[1:2],v4[1:2],Fmatrix,gkpqmap,gkmqmap,T1,T2,wave_diff,allowedq,constq)+get_Velement(q2index,v1[1:2],v2[1:2],v4[1:2],v3[1:2],Fmatrix,gkpqmap,gkmqmap,T1,T2,wave_diff,allowedq,constq)
   end
   end

   for ja in eachindex(Aindexset)
    Bmatrix[ja,:]=Bmatrixvec[ja]
   end



 return  (Amatrix+Amatrix')/2,(AmQmatrix+AmQmatrix')/2,Bmatrix
end