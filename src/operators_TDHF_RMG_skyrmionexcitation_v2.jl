










function get_indexset(wave_diff::Vector{Vector{Int64}},num_bandbelow::Int64,num_bandup::Int64,allowedq::Vector{Vector{Int64}},Nq::Int64,bigQ::Vector{Int64})::Tuple{Vector{Vector{Int64}},Vector{Vector{Vector{Int64}}},Vector{Vector{Vector{Int64}}},Vector{Vector{Vector{Int64}}},Array{Int64},Array{Int64}}
        
  
 FLindex=1
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


function get_Fmatrix(Bandvector::Vector{Vector{Int64}},spinor_set::Matrix{Vector{ComplexF64}},
                          HF_eigenvector::Vector{Matrix{ComplexF64}},
                          wave::Vector{Vector{Int}},wave_diff::Vector{Vector{Int}},NL::Int64)::Array{ComplexF64}
    dimension=length(wave)*2*NL
    Fmatrix=zeros(ComplexF64,length(Bandvector),length(Bandvector),length(wave_diff))
    eigenvector_single=zeros(ComplexF64,dimension,length(Bandvector))

    for ja in eachindex(Bandvector), jb in eachindex(wave)
        eigenvector_single[2*NL*(jb-1)+1:2*NL*jb,ja]=spinor_set[Bandvector[ja][1],jb]*HF_eigenvector[Bandvector[ja][1]][jb,Bandvector[ja][2]] # Rewrite this so that eigenvector_single has the right dimension
    end       


   println("inside F, I am here")


    diff_eigenvector=zeros(ComplexF64,dimension,length(Bandvector),length(wave_diff))
    for jc in eachindex(wave_diff), jd in eachindex(wave)
        pos=findfirst(item->item==wave[jd]+wave_diff[jc],wave)
        if pos≠nothing 
            diff_eigenvector[2*NL*(jd-1)+1:2*NL*jd,:,jc]=eigenvector_single[2*NL*(pos-1)+1:2*NL*pos,:]
        end 
    end
 

  
    Threads.@threads for ja in eachindex(Bandvector) 
      for jb in eachindex(Bandvector), jc in eachindex(wave_diff)
        Fmatrix[ja,jb,jc]=diff_eigenvector[:,ja,jc]'*eigenvector_single[:,jb]
      end
    end

    println("inside F, I am here")
    flush(stdout)
  
    pos=findfirst(item->item==[0,0],wave_diff)
    for ja in eachindex(Bandvector)
        Fmatrix[ja,ja,pos]=1.0+0.0*im
    end
    
  return Fmatrix
end


function get_Velement(qindex::Int64,v1::Vector{Int64},v2::Vector{Int64},v3::Vector{Int64},v4::Vector{Int64}
                    ,Fmatrix::Array{ComplexF64},gkpqmap::Array{Int64},gkmqmap::Array{Int64}
                     ,T1::Vector{Float64},T2::Vector{Float64},wave_diff::Vector{Vector{Int64}}
                      ,allowedq::Vector{Vector{Int64}},Area::Float64,contact_strength::Float64,ϵr::Float64)::ComplexF64
    Velement=0 #v1v2v3v4 contains two informatin, which band and which momentum
   for ja in eachindex(wave_diff)
       qvec=allowedq[qindex]+wave_diff[ja]
       if  (gkpqmap[v3[2],qindex,ja]≠0) && (gkmqmap[v4[2],qindex,ja]≠0)
          Velement+=(Coulomb(qvec,T1,T2)/ϵr+contact_strength)*Fmatrix[v1[1],v3[1],gkpqmap[v3[2],qindex,ja]]*Fmatrix[v2[1],v4[1],gkmqmap[v4[2],qindex,ja]]
       end
   end

   return Velement/Area
end


function Coulomb(k::Vector{Int},T1::Vector{Float64},T2::Vector{Float64})::Float64
   D=25
   return k==[0,0] ? D*9047.5636 : tanh(norm([T1 T2]*k*D))/norm(k[1]*T1+k[2]*T2)*9047.5636
end



function Construct_Amatrix(Aindexset::Vector{Vector{Vector{Int64}}},AmQindexset::Vector{Vector{Vector{Int64}}},B2indexset::Vector{Vector{Vector{Int64}}},
                  allowedq::Vector{Vector{Int64}},Fmatrix::Array{ComplexF64},gkpqmap::Array{Int64},gkmqmap::Array{Int64},
                  T1::Vector{Float64},T2::Vector{Float64},wave_diff::Vector{Vector{Int64}},
                  HF_eigenvalue::Vector{Vector{Float64}},Nq::Int64,Area::Float64,contact_strength::Float64,ϵr::Float64)::Tuple{Matrix{ComplexF64},Matrix{ComplexF64},Matrix{ComplexF64}}
   


   
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
   Amatrixvec[ja][jb]=-get_Velement(q1index,v1[1:2],v2[1:2],v3[1:2],v4[1:2],Fmatrix,gkpqmap,gkmqmap,T1,T2,wave_diff,allowedq,Area,contact_strength,ϵr)+get_Velement(q2index,v1[1:2],v2[1:2],v4[1:2],v3[1:2],Fmatrix,gkpqmap,gkmqmap,T1,T2,wave_diff,allowedq,Area,contact_strength,ϵr)
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
   AmQmatrixvec[ja][jb]=-get_Velement(q1index,v1[1:2],v2[1:2],v3[1:2],v4[1:2],Fmatrix,gkpqmap,gkmqmap,T1,T2,wave_diff,allowedq,Area,contact_strength,ϵr)+get_Velement(q2index,v1[1:2],v2[1:2],v4[1:2],v3[1:2],Fmatrix,gkpqmap,gkmqmap,T1,T2,wave_diff,allowedq,Area,contact_strength,ϵr)
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
     Bmatrixvec[ja][jb]=-get_Velement(q1index,v1[1:2],v2[1:2],v3[1:2],v4[1:2],Fmatrix,gkpqmap,gkmqmap,T1,T2,wave_diff,allowedq,Area,contact_strength,ϵr)+get_Velement(q2index,v1[1:2],v2[1:2],v4[1:2],v3[1:2],Fmatrix,gkpqmap,gkmqmap,T1,T2,wave_diff,allowedq,Area,contact_strength,ϵr)
   end
   end

   for ja in eachindex(Aindexset)
    Bmatrix[ja,:]=Bmatrixvec[ja]
   end



 return  (Amatrix+Amatrix')/2,(AmQmatrix+AmQmatrix')/2,Bmatrix
end



function get_wavediff(gcutoff::Float64,T1::Vector{Float64},T2::Vector{Float64},b1::Vector{Float64},b2::Vector{Float64})

   

    
    
    
    b1T=Int.(round.(inv([T1 T2])*b1))
    b2T=Int.(round.(inv([T1 T2])*b2))
    
    
    
    wavediff=Vector{Int64}[]
    cutoff=25
    cutoffstandard=2*gcutoff*norm(b1)
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        gtest=ja*b1+jb*b2;
        if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
            push!(wavediff,ja*b1T+jb*b2T)
        end
    end

    return wavediff
end