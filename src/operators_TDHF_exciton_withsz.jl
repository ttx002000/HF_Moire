









function Coulomb(kvec::Vector{Float64},dis::Float64)::Float64
  
  if norm(kvec)==0.0
      return 9047.5636*(-dis)

    else
     return 9047.5636/norm(kvec)*exp(-norm(kvec)*dis)

  end
  
end
  

function Coulomb_matrix(kvec::Vector{Float64},zpos::Matrix{Float64},NL::Int)
  sublattice_num=2*NL
  valley_num=2
  layer_num=2
  
  Coulomb_matrix=zeros(Float64,sublattice_num,valley_num,layer_num,sublattice_num,valley_num,layer_num)
  for ja in 1:sublattice_num, jb in 1:layer_num, jc in 1:sublattice_num, jd in 1:layer_num
    Coulomb_matrix[ja,:,jb,jc,:,jd].=Coulomb(kvec,abs(zpos[jb,ja]-zpos[jd,jc]))
  end
  Coulomb_matrix=reshape(Coulomb_matrix,sublattice_num*valley_num*layer_num,sublattice_num*valley_num*layer_num)

  return Coulomb_matrix
end






function get_indexset(num_bandbelow::Int64,num_bandup::Int64,mode::Int64,bigQ::Vector{Int64},TDHF_k_set::Vector{Vector{Float64}},TDHF_k_index::Vector{Vector{Int}},TDHF_k_pos::Vector{Int},k_index::Vector{Vector{Int}},k_set::Vector{Vector{Float64}})::Tuple{Vector{Vector{Int64}},Vector{Vector{Vector{Int64}}},Vector{Vector{Vector{Int64}}},Vector{Vector{Vector{Int64}}}}
        

  valley_num=2
  layer_num=2
  spin_num=2

  #Let's input three variables, TDHF_k_set, TDHF_k_index and TDHF_k_pos. TDHF_k_pos[ja] is the index of the corresponding vector in the original set
 

 bandnum=valley_num*layer_num 

 FLindex=Int(round(bandnum/2)) # the higest occupied band index
  println(FLindex)
 Bandvector=Vector{Int64}[]

 for ja in eachindex(TDHF_k_set), sih in 1:spin_num, sip in 1:spin_num
  if (mode==0 && sih==sip)|| (mode==1 && sih≠sip) || (mode==2)
   for jb in FLindex-num_bandbelow+1:FLindex
       push!(Bandvector,[TDHF_k_pos[ja],jb,sih]) #I think this is the hole band
   end
  
      kpQ=TDHF_k_index[ja]+bigQ
      kmQ=TDHF_k_index[ja]-bigQ

      kpQpos=findfirst(item->item==kpQ,k_index)
      kmQpos=findfirst(item->item==kmQ,k_index)

      if kpQpos≠nothing 
        for jb in FLindex+1:FLindex+num_bandup
          push!(Bandvector,[kpQpos,jb,sip])#I think this is the electron band
        end
      end

      if kmQpos≠nothing 
        for jb in FLindex+1:FLindex+num_bandup
            push!(Bandvector,[kmQpos,jb,sip])
        end
      end
    end


  end
 Bandvector=sort(unique(Bandvector))




 Aindexset=Vector{Vector{Int64}}[]
 for ja in eachindex(TDHF_k_index), sih in 1:spin_num, sip in 1:spin_num
   if (mode==0 && sih==sip)|| (mode==1 && sih≠sip) || (mode==2)
     kpQ=TDHF_k_index[ja]+bigQ
     kpQpos=findfirst(item->item==kpQ,k_index)
    if kpQpos≠nothing
     for hband in FLindex-num_bandbelow+1:FLindex, pband in FLindex+1:FLindex+num_bandup
         pos1=findfirst(item->item==[kpQpos,pband,sip],Bandvector)
         pos2=findfirst(item->item==[TDHF_k_pos[ja],hband,sih],Bandvector)
         push!(Aindexset,[[pos1,kpQpos,pband,sip],[pos2,TDHF_k_pos[ja],hband,sih]])
     end
    end

  end
 end
 
 AmQindexset=Vector{Vector{Int64}}[]
 for ja in eachindex(TDHF_k_index), sih in 1:spin_num, sip in 1:spin_num
     if (mode==0 && sih==sip)|| (mode==1 && sih≠sip) || (mode==2)
     kmQ=TDHF_k_index[ja]-bigQ
     kmQpos=findfirst(item->item==kmQ,k_index)
    if kmQpos≠nothing
     for hband in FLindex-num_bandbelow+1:FLindex, pband in FLindex+1:FLindex+num_bandup
         pos1=findfirst(item->item==[kmQpos,pband,sip],Bandvector)
         pos2=findfirst(item->item==[TDHF_k_pos[ja],hband,sih],Bandvector)
         push!(AmQindexset,[[pos1,kmQpos,pband,sip],[pos2,TDHF_k_pos[ja],hband,sih]])
     end
    end
  end
  end


 B2indexset=Vector{Vector{Int64}}[]
 for ja in eachindex(TDHF_k_index), sih in 1:spin_num, sip in 1:spin_num
    if (mode==0 && sih==sip)|| (mode==1 && sih≠sip) || (mode==2)
    kmQ=TDHF_k_index[ja]-bigQ
    kmQpos=findfirst(item->item==kmQ,k_index)
    if kmQpos≠nothing
     for hband in FLindex-num_bandbelow+1:FLindex, pband in FLindex+1:FLindex+num_bandup
        pos1=findfirst(item->item==[kmQpos,pband,sip],Bandvector)
        pos2=findfirst(item->item==[TDHF_k_pos[ja],hband,sih],Bandvector)
        push!(B2indexset,[[pos1,kmQpos,pband,sip],[pos2,TDHF_k_pos[ja],hband,sih]])
     end
   end
  end
 end

 return Bandvector,Aindexset,AmQindexset,B2indexset
end





function get_Velement(v1::Vector{Int64},v2::Vector{Int64},v3::Vector{Int64},v4::Vector{Int64},HF_eigenvectors_sub::Array{ComplexF64},Vmatrix::Matrix{Float64})::ComplexF64
 
  # the first index of v1 is which index in the Bandvector, the second index is which momentum, the third index is which band, the 4th is the spin
  #Vmatrix should


  #HF_eigenvectors_sublatticebasis_TDHF=zeros(ComplexF64,sublattice_num*valley_num*layer_num,valley_num*layer_num,spin_num,length(k_set))

  if v1[4]==v3[4] && v2[4]==v4[4]
   p1=transpose(conj.(HF_eigenvectors_sub[:,v1[3],v1[4],v1[2]]) .*HF_eigenvectors_sub[:,v3[3],v3[4],v3[2]]) #We should put this into the sublattice basis
   p2=conj.(HF_eigenvectors_sub[:,v2[3],v2[4],v2[2]]) .*HF_eigenvectors_sub[:,v4[3],v4[4],v4[2]]
   Velement=p1*Vmatrix*p2
    return Velement
  else
    return 0.0

  end
   

  
end







function get_Hunds(v1::Vector{Int64},v2::Vector{Int64},v3::Vector{Int64},v4::Vector{Int64},HF_eigenvectors_sub::Array{ComplexF64},NL::Int,ppmatrix::Array{ComplexF64})
   
  valley_num=2
  sublattice_num=2*NL
  layer_num=2


  reshaped_size=(sublattice_num,valley_num,layer_num)
  
  u1=reshape(HF_eigenvectors_sub[:,v1[3],v1[4],v1[2]],reshaped_size)
  u2=reshape(HF_eigenvectors_sub[:,v2[3],v2[4],v2[2]],reshaped_size)
  u3=reshape(HF_eigenvectors_sub[:,v3[3],v3[4],v3[2]],reshaped_size)
  u4=reshape(HF_eigenvectors_sub[:,v4[3],v4[4],v4[2]],reshaped_size)



  Velement=0.0 
  for li in 1:layer_num, vone in 1:valley_num, vtwo in 1:valley_num
     if vone≠vtwo
       f1=u1[:,vone,li]'*u3[:,vone,li]
       f2=u2[:,vtwo,li]'*u4[:,vtwo,li]
       Velement+=f1*f2
     end   
  end
  
  return Velement*ppmatrix[v1[4],v2[4],v3[4],v4[4]]

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











function Construct_Amatrix(Aindexset::Vector{Vector{Vector{Int64}}},AmQindexset::Vector{Vector{Vector{Int64}}},B2indexset::Vector{Vector{Vector{Int64}}},HF_eigenvalues::Array{Float64},HF_eigenvectors_sub::Array{ComplexF64},TDHF_Area::Float64,ϵr::Float64,NL::Int,k_set::Vector{Vector{Float64}},JH::Float64,z_pos::Array{Float64})::Tuple{Matrix{ComplexF64},Matrix{ComplexF64},Matrix{ComplexF64}}

    


  layer_num=2
  valley_num=2
  sublattice_num=2*NL
  spin_num=2

   paulimatrix=[[0 1;1 0],[0 -im;im 0],[1 0;0 -1]]
   ppmatrix=zeros(ComplexF64,spin_num,spin_num,spin_num,spin_num)
   for ppi in 1:3
     ppmatrix+=[paulimatrix[ppi][sone,sthree]*paulimatrix[ppi][stwo,sfour] for sone in 1:spin_num, stwo in 1:spin_num, sthree in 1:spin_num, sfour in 1:spin_num] 
   end


    Vmatrix_set=Matrix{Matrix{Float64}}(undef,length(k_set),length(k_set))


   for ja in eachindex(k_set), jb in eachindex(k_set)
     Vmatrix_set[ja,jb]=Coulomb_matrix(k_set[ja]-k_set[jb],z_pos,NL) 
   end

  Amatrix=zeros(ComplexF64,length(Aindexset),length(Aindexset))
  Amatrixvec=[zeros(ComplexF64,length(Aindexset)) for _ in eachindex(Aindexset)]


  for ja in eachindex(Aindexset)
  for  jb in eachindex(Aindexset)
   v1=Aindexset[ja][1]
   v2=Aindexset[jb][2]
   v3=Aindexset[jb][1]
   v4=Aindexset[ja][2]
  


   #V1=Coulomb_matrix(k_set[v1[2]]-k_set[v3[2]])
    V1=Vmatrix_set[v1[2],v3[2]]
   #V2=Coulomb_matrix(k_set[v1[2]]-k_set[v4[2]])
    V2=Vmatrix_set[v1[2],v4[2]]

   Amatrixvec[ja][jb]+=(-get_Velement(v1,v2,v3,v4,HF_eigenvectors_sub,V1)+get_Velement(v1,v2,v4,v3,HF_eigenvectors_sub,V2))*1/(TDHF_Area*ϵr)
    if !(abs(JH)==0.0)
        Amatrixvec[ja][jb]+=(-get_Hunds(v1,v2,v3,v4,HF_eigenvectors_sub,NL,ppmatrix)+get_Hunds(v1,v2,v4,v3,HF_eigenvectors_sub,NL,ppmatrix))*JH/(TDHF_Area)
    end
  
  end
 end
   
  for ja in eachindex(Aindexset)
    Amatrix[ja,:]+=Amatrixvec[ja]
 end

 Amatrixvec=nothing
 for ja in eachindex(Aindexset)
   v1=Aindexset[ja][1]
   v2=Aindexset[ja][2]
   Amatrix[ja,ja]+=HF_eigenvalues[v1[3],v1[4],v1[2]]-HF_eigenvalues[v2[3],v2[4],v2[2]]
 end

 println("finishA")
 flush(stdout)
 
 AmQmatrix=zeros(ComplexF64,length(AmQindexset),length(AmQindexset))
 AmQmatrixvec=[zeros(ComplexF64,length(AmQindexset)) for _ in eachindex(AmQindexset)]

 for ja in eachindex(AmQindexset)
  
 for jb in eachindex(AmQindexset)
   v1=AmQindexset[ja][1]
   v2=AmQindexset[jb][2]
   v3=AmQindexset[jb][1]
   v4=AmQindexset[ja][2]


   V1=Vmatrix_set[v1[2],v3[2]]
   V2=Vmatrix_set[v1[2],v4[2]]

   AmQmatrixvec[ja][jb]+=(-get_Velement(v1,v2,v3,v4,HF_eigenvectors_sub,V1)+get_Velement(v1,v2,v4,v3,HF_eigenvectors_sub,V2))*1/(TDHF_Area*ϵr)
     if !(abs(JH)==0.0)
   AmQmatrixvec[ja][jb]+=(-get_Hunds(v1,v2,v3,v4,HF_eigenvectors_sub,NL,ppmatrix)+get_Hunds(v1,v2,v4,v3,HF_eigenvectors_sub,NL,ppmatrix))*JH/(TDHF_Area)
     end
  end
 end

 for ja in eachindex(AmQindexset)
    AmQmatrix[ja,:]+=AmQmatrixvec[ja]
 end

 AmQmatrixvec=nothing
 println("finishA")
 flush(stdout)


  for ja in eachindex(AmQindexset)
   v1=AmQindexset[ja][1]
   v2=AmQindexset[ja][2]
   AmQmatrix[ja,ja]+=HF_eigenvalues[v1[3],v1[4],v1[2]]-HF_eigenvalues[v2[3],v2[4],v2[2]]
 end




  Bmatrix=zeros(ComplexF64,length(Aindexset),length(B2indexset))
  Bmatrixvec=[zeros(ComplexF64,length(B2indexset)) for _ in eachindex(Aindexset)]

  for ja in eachindex(Aindexset)
   for jb in eachindex(B2indexset)
     v1=Aindexset[ja][1]
     v2=B2indexset[jb][1]
     v3=B2indexset[jb][2]
     v4=Aindexset[ja][2]
     
   
     #V1=Coulomb_matrix(k_set[v1[2]]-k_set[v3[2]],zpos)
 
     #V2=Coulomb_matrix(k_set[v1[2]]-k_set[v4[2]],zpos)
     
     
      V1=Vmatrix_set[v1[2],v3[2]]
      V2=Vmatrix_set[v1[2],v4[2]]

     Bmatrixvec[ja][jb]+=(-get_Velement(v1,v2,v3,v4,HF_eigenvectors_sub,V1)+get_Velement(v1,v2,v4,v3,HF_eigenvectors_sub,V2))*1/(TDHF_Area*ϵr)
     if !(abs(JH)==0.0)

      Bmatrixvec[ja][jb]+=(-get_Hunds(v1,v2,v3,v4,HF_eigenvectors_sub,NL,ppmatrix)+get_Hunds(v1,v2,v4,v3,HF_eigenvectors_sub,NL,ppmatrix))*JH/(TDHF_Area)
     end
    end
   end

   for ja in eachindex(Aindexset)
    Bmatrix[ja,:]+=Bmatrixvec[ja]
   end

   Bmatrixvec=nothing
   println("finishB")
   flush(stdout)

 return  (Amatrix+Amatrix')/2,(AmQmatrix+AmQmatrix')/2,Bmatrix
end



