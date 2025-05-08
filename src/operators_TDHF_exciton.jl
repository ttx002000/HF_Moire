





function Coulomb(dis::Float64,kvec::Vector{Float64})::Float64
  
  if norm(kvec)==0.0
      return 9047.5636*(-dis)
  
    else
     return 9047.5636/norm(kvec)*exp(-norm(kvec)*dis)
    
  end
  



end
  




function get_indexset(num_bandbelow::Int64,num_bandup::Int64,bigQ::Vector{Int64},TDHF_k_set::Vector{Vector{Float64}},TDHF_k_index::Vector{Vector{Int}},TDHF_k_pos::Vector{Int},k_index::Vector{Vector{Int}},k_set::Vector{Vector{Float64}})::Tuple{Vector{Vector{Int64}},Vector{Vector{Vector{Int64}}},Vector{Vector{Vector{Int64}}},Vector{Vector{Vector{Int64}}}}
        
  spin_num=2
  valley_num=2
  sublattice_num=Int(2*NL)

  #Let's input three variables, TDHF_k_set, TDHF_k_index and TDHF_k_pos. TDHF_k_pos[ja] is the index of the corresponding vector in the original set


 dimension=spin_num*valley_num*sublattice_num   
 FLindex=Int(round(dimension/2)) # the higest occupied band index

 Bandvector=Vector{Int64}[]

 for ja in eachindex(TDHF_k_set)
   for jb in FLindex-num_bandbelow+1:FLindex
       push!(Bandvector,[TDHF_k_pos[ja],jb]) #I think this is the hole band
   end
  
      kpQ=TDHF_k_index[ja]+bigQ
      kmQ=TDHF_k_index[ja]-bigQ

      kpQpos=findfirst(item->item==kpQ,k_index)
      kmQpos=findfirst(item->item==kmQ,k_index)

      if kpQpos≠nothing 
        for jb in FLindex+1:FLindex+num_bandup
          push!(Bandvector,[kpQpos,jb])#I think this is the electron band
        end
      end

      if kmQpos≠nothing 
        for jb in FLindex+1:FLindex+num_bandup
            push!(Bandvector,[kmQpos,jb])
        end
      end
  
  end
 Bandvector=sort(unique(Bandvector))




 Aindexset=Vector{Vector{Int64}}[]
 for ja in eachindex(TDHF_k_index)
     kpQ=TDHF_k_index[ja]+bigQ
     kpQpos=findfirst(item->item==kpQ,k_index)
    if kpQpos≠nothing
     for hband in FLindex-num_bandbelow+1:FLindex, pband in FLindex+1:FLindex+num_bandup
         pos1=findfirst(item->item==[kpQpos,pband],Bandvector)
         pos2=findfirst(item->item==[TDHF_k_pos[ja],hband],Bandvector)
         push!(Aindexset,[[pos1,kpQpos,pband],[pos2,TDHF_k_pos[ja],hband]])
     end
    end


 end
 
 AmQindexset=Vector{Vector{Int64}}[]
 for ja in eachindex(TDHF_k_index)
     kmQ=TDHF_k_index[ja]-bigQ
     kmQpos=findfirst(item->item==kmQ,k_index)
    if kmQpos≠nothing
     for hband in FLindex-num_bandbelow+1:FLindex, pband in FLindex+1:FLindex+num_bandup
         pos1=findfirst(item->item==[kmQpos,pband],Bandvector)
         pos2=findfirst(item->item==[TDHF_k_pos[ja],hband],Bandvector)
         push!(AmQindexset,[[pos1,kmQpos,pband],[pos2,TDHF_k_pos[ja],hband]])
     end
    end

 end


 B2indexset=Vector{Vector{Int64}}[]
 for ja in eachindex(TDHF_k_index)
    kmQ=TDHF_k_index[ja]-bigQ
    kmQpos=findfirst(item->item==kmQ,k_index)
    if kmQpos≠nothing
     for hband in FLindex-num_bandbelow+1:FLindex, pband in FLindex+1:FLindex+num_bandup
        pos1=findfirst(item->item==[kmQpos,pband],Bandvector)
        pos2=findfirst(item->item==[TDHF_k_pos[ja],hband],Bandvector)
        push!(B2indexset,[[pos1,kmQpos,pband],[pos2,TDHF_k_pos[ja],hband]])
     end
   end
 end

 return Bandvector,Aindexset,AmQindexset,B2indexset
end






function get_formfactors(NL::Int,k_set::Vector{Vector{Float64}},ildis::Float64)
  
  valley_num=2
  spin_num=2
  layer_num=2
  sublattice_num=Int(2*NL)
  dimension=valley_num*spin_num*sublattice_num*layer_num


  z_pos=zeros(Float64,layer_num,sublattice_num)
  z_pos[1,:]=0.335*[i for i in 0:NL-1 for _ in 1:2]
  z_pos[2,:]=0.335*[i for i in -NL+1:0 for _ in 1:2].-ildis
  z_pos=reshape(z_pos,layer_num*sublattice_num)
  sublayer_i=zeros(Int,valley_num,spin_num,sublattice_num*layer_num)
  for ja in 1:valley_num,jb in 1:spin_num, jc in 1:sublattice_num*layer_num
    sublayer_i[ja,jb,jc]=jc
  end
  sublayer_i=reshape(sublayer_i,dimension)

  fcmatrix=Matrix{Matrix{Float64}}(undef,length(k_set),length(k_set))


  tic=time()

  Threads.@threads for ja in 1:length(k_set)
   for jb in 1:ja
    fcmatrix[ja,jb]=zeros(Float64,dimension,dimension)
        for jc in 1:dimension, jd in 1:jc
          fcmatrix[ja,jb][jc,jd]+=Coulomb(abs(z_pos[sublayer_i[jc]]-z_pos[sublayer_i[jd]]),k_set[ja]-k_set[jb])
        end

        for jc in 1:dimension, jd in jc+1:dimension
          fcmatrix[ja,jb][jc,jd]+=fcmatrix[ja,jb][jd,jc]
        end
   end
 end



 toc=time()
 println("formfactorstime",toc-tic) 


  return fcmatrix
  


 
    

end







function get_Velement(Vmatrix::Matrix{Float64},v1::Vector{Int64},v2::Vector{Int64},v3::Vector{Int64},v4::Vector{Int64},HF_eigenvectors::Array{ComplexF64})::ComplexF64
 
  # the first index of v1 is which index in the Bandvector, the second index is which momentum, the third index is which band
  #Vmatrix should

   p1=transpose((conj.(HF_eigenvectors[:,v1[3],v1[2]])) .* HF_eigenvectors[:,v3[3],v3[2]])
   p2=(conj.(HF_eigenvectors[:,v2[3],v2[2]])) .* HF_eigenvectors[:,v4[3],v4[2]]
   Velement=p1*(Vmatrix)*p2
   

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




function construct_Vmatrix(k_vec::Vector{Float64},z_pos::Vector{Float64},NL::Int)

  valley_num=2
  spin_num=2
  layer_num=2
  sublattice_num=2*NL
  dimension=valley_num*spin_num*layer_num*sublattice_num


  s1=zeros(Float64,layer_num*sublattice_num,layer_num*sublattice_num)
  for ja in 1:layer_num*sublattice_num, jb in 1:layer_num*sublattice_num
     s1[ja,jb]=Coulomb(z_pos[ja]-z_pos[jb],k_vec)
    
  end

  Vmatrix=zeros(valley_num*spin_num,layer_num*sublattice_num,valley_num*spin_num,layer_num*sublattice_num)
  for ja in 1:valley_num*spin_num, jb in 1:valley_num*spin_num
    Vmatrix[ja,:,jb,:]=s1
  end

  return(reshape(Vmatrix,dimension,dimension))

end









function Construct_Amatrix(Aindexset::Vector{Vector{Vector{Int64}}},AmQindexset::Vector{Vector{Vector{Int64}}},B2indexset::Vector{Vector{Vector{Int64}}},HF_eigenvalues::Matrix{Float64},HF_eigenvectors::Array{ComplexF64},TDHF_Area::Float64,ϵr::Float64,ildis::Float64,NL::Int,k_set::Vector{Vector{Float64}})::Tuple{Matrix{ComplexF64},Matrix{ComplexF64},Matrix{ComplexF64}}

    

  valley_num=2
  spin_num=2
  layer_num=2
  sublattice_num=Int(2*NL)
  dimension=valley_num*spin_num*sublattice_num*layer_num


  z_pos=zeros(Float64,layer_num,sublattice_num)
  z_pos[1,:]=0.335*[i for i in 0:NL-1 for _ in 1:2]
  z_pos[2,:]=0.335*[i for i in -NL+1:0 for _ in 1:2].-ildis
  z_pos=reshape(z_pos,layer_num*sublattice_num)
   
    Amatrix=zeros(ComplexF64,length(Aindexset),length(Aindexset))
    Amatrixvec=[zeros(ComplexF64,length(Aindexset)) for _ in eachindex(Aindexset)]

 Threads.@threads for ja in eachindex(Aindexset)
  for  jb in eachindex(Aindexset)
   v1=Aindexset[ja][1]
   v2=Aindexset[jb][2]
   v3=Aindexset[jb][1]
   v4=Aindexset[ja][2]
  


   V1matrix=construct_Vmatrix(k_set[v1[2]]-k_set[v3[2]],z_pos,NL)

   V2matrix=construct_Vmatrix(k_set[v1[2]]-k_set[v4[2]],z_pos,NL)

   
   Amatrixvec[ja][jb]+=-get_Velement(V1matrix,v1,v2,v3,v4,HF_eigenvectors)+get_Velement(V2matrix,v1,v2,v4,v3,HF_eigenvectors)
  end
 end
   

 for ja in eachindex(Aindexset)
    Amatrix[ja,:]=Amatrixvec[ja]*1/(TDHF_Area*ϵr)
 end

 for ja in eachindex(Aindexset)
   v1=Aindexset[ja][1]
   v2=Aindexset[ja][2]
   Amatrix[ja,ja]+=HF_eigenvalues[v1[3],v1[2]]-HF_eigenvalues[v2[3],v2[2]]
 end

 println("finishA")
 flush(stdout)
 
 AmQmatrix=zeros(ComplexF64,length(AmQindexset),length(AmQindexset))
 AmQmatrixvec=[zeros(ComplexF64,length(AmQindexset)) for _ in eachindex(AmQindexset)]

 Threads.@threads for ja in eachindex(AmQindexset)
  
 for jb in eachindex(AmQindexset)
   v1=AmQindexset[ja][1]
   v2=AmQindexset[jb][2]
   v3=AmQindexset[jb][1]
   v4=AmQindexset[ja][2]
   
   V1matrix=construct_Vmatrix(k_set[v1[2]]-k_set[v3[2]],z_pos,NL)

   V2matrix=construct_Vmatrix(k_set[v1[2]]-k_set[v4[2]],z_pos,NL)

   AmQmatrixvec[ja][jb]+=-get_Velement(V1matrix,v1,v2,v3,v4,HF_eigenvectors)+get_Velement(V2matrix,v1,v2,v4,v3,HF_eigenvectors)
 end
 end

 for ja in eachindex(AmQindexset)
    AmQmatrix[ja,:]=AmQmatrixvec[ja]*1/(TDHF_Area*ϵr)
 end
 println("finishA")
 flush(stdout)


  for ja in eachindex(AmQindexset)
   v1=AmQindexset[ja][1]
   v2=AmQindexset[ja][2]
   AmQmatrix[ja,ja]+=HF_eigenvalues[v1[3],v1[2]]-HF_eigenvalues[v2[3],v2[2]]
 end




  Bmatrix=zeros(ComplexF64,length(Aindexset),length(B2indexset))
  Bmatrixvec=[zeros(ComplexF64,length(B2indexset)) for _ in eachindex(Aindexset)]

  Threads.@threads for ja in eachindex(Aindexset)
   for jb in eachindex(B2indexset)
     v1=Aindexset[ja][1]
     v2=B2indexset[jb][1]
     v3=B2indexset[jb][2]
     v4=Aindexset[ja][2]
     
      V1matrix=construct_Vmatrix(k_set[v1[2]]-k_set[v3[2]],z_pos,NL)

     V2matrix=construct_Vmatrix(k_set[v1[2]]-k_set[v4[2]],z_pos,NL)

     Bmatrixvec[ja][jb]+=-get_Velement(V1matrix,v1,v2,v3,v4,HF_eigenvectors)+get_Velement(V2matrix,v1,v2,v4,v3,HF_eigenvectors)
   end
   end

   for ja in eachindex(Aindexset)
    Bmatrix[ja,:]+=Bmatrixvec[ja]*1/(TDHF_Area*ϵr)
   end

   println("finishB")
   flush(stdout)

 return  (Amatrix+Amatrix')/2,(AmQmatrix+AmQmatrix')/2,Bmatrix
end







#=
function Construct_Amatrix(Aindexset::Vector{Vector{Vector{Int64}}},AmQindexset::Vector{Vector{Vector{Int64}}},B2indexset::Vector{Vector{Vector{Int64}}},HF_eigenvalues::Matrix{Float64},HF_eigenvectors::Array{ComplexF64},fcmatrix::Matrix{Matrix{Float64}},TDHF_Area::Float64,ϵr::Float64)::Tuple{Matrix{ComplexF64},Matrix{ComplexF64},Matrix{ComplexF64}}


   
    Amatrix=zeros(ComplexF64,length(Aindexset),length(Aindexset))
    Amatrixvec=[zeros(ComplexF64,length(Aindexset)) for _ in eachindex(Aindexset)]

 Threads.@threads for ja in eachindex(Aindexset)
  for  jb in eachindex(Aindexset)
   v1=Aindexset[ja][1]
   v2=Aindexset[jb][2]
   v3=Aindexset[jb][1]
   v4=Aindexset[ja][2]
  
   if v1[2]>=v3[2]
    V1matrix=fcmatrix[v1[2],v3[2]]
   else
    V1matrix=fcmatrix[v3[2],v1[2]]
   end

   
   if v1[2]>=v4[2]
    V2matrix=fcmatrix[v1[2],v4[2]]
   else
    V2matrix=fcmatrix[v4[2],v1[2]]
   end
   
   
   Amatrixvec[ja][jb]=-get_Velement(V1matrix,v1,v2,v3,v4,HF_eigenvectors)+get_Velement(V2matrix,v1,v2,v4,v3,HF_eigenvectors)
  end
 end
   

 for ja in eachindex(Aindexset)
    Amatrix[ja,:]=Amatrixvec[ja]*1/(TDHF_Area*ϵr)
 end

 for ja in eachindex(Aindexset)
   v1=Aindexset[ja][1]
   v2=Aindexset[ja][2]
   Amatrix[ja,ja]+=HF_eigenvalues[v1[3],v1[2]]-HF_eigenvalues[v2[3],v2[2]]
 end

 println("finishA")
 flush(stdout)
 
 AmQmatrix=zeros(ComplexF64,length(AmQindexset),length(AmQindexset))
 AmQmatrixvec=[zeros(ComplexF64,length(AmQindexset)) for _ in eachindex(AmQindexset)]

 Threads.@threads for ja in eachindex(AmQindexset)
 for jb in eachindex(AmQindexset)
   v1=AmQindexset[ja][1]
   v2=AmQindexset[jb][2]
   v3=AmQindexset[jb][1]
   v4=AmQindexset[ja][2]
   if v1[2]>=v3[2]
    V1matrix=fcmatrix[v1[2],v3[2]]
   else
    V1matrix=fcmatrix[v3[2],v1[2]]
   end

   
   if v1[2]>=v4[2]
    V2matrix=fcmatrix[v1[2],v4[2]]
   else
    V2matrix=fcmatrix[v4[2],v1[2]]
   end
   AmQmatrixvec[ja][jb]=-get_Velement(V1matrix,v1,v2,v3,v4,HF_eigenvectors)+get_Velement(V2matrix,v1,v2,v4,v3,HF_eigenvectors)
 end
 end

 for ja in eachindex(AmQindexset)
    AmQmatrix[ja,:]=AmQmatrixvec[ja]*1/(TDHF_Area*ϵr)
 end
 println("finishA")
 flush(stdout)


  for ja in eachindex(AmQindexset)
   v1=AmQindexset[ja][1]
   v2=AmQindexset[ja][2]
   AmQmatrix[ja,ja]+=HF_eigenvalues[v1[3],v1[2]]-HF_eigenvalues[v2[3],v2[2]]
 end




  Bmatrix=zeros(ComplexF64,length(Aindexset),length(B2indexset))
  Bmatrixvec=[zeros(ComplexF64,length(B2indexset)) for _ in eachindex(Aindexset)]

  Threads.@threads for ja in eachindex(Aindexset)
   for jb in eachindex(B2indexset)
     v1=Aindexset[ja][1]
     v2=B2indexset[jb][1]
     v3=B2indexset[jb][2]
     v4=Aindexset[ja][2]
     if v1[2]>=v3[2]
      V1matrix=fcmatrix[v1[2],v3[2]]
     else
      V1matrix=fcmatrix[v3[2],v1[2]]
     end
  
     
     if v1[2]>=v4[2]
      V2matrix=fcmatrix[v1[2],v4[2]]
     else
      V2matrix=fcmatrix[v4[2],v1[2]]
     end
     Bmatrixvec[ja][jb]=-get_Velement(V1matrix,v1,v2,v3,v4,HF_eigenvectors)+get_Velement(V2matrix,v1,v2,v4,v3,HF_eigenvectors)
   end
   end

   for ja in eachindex(Aindexset)
    Bmatrix[ja,:]=Bmatrixvec[ja]*1/(TDHF_Area*ϵr)
   end

   println("finishB")
   flush(stdout)

 return  (Amatrix+Amatrix')/2,(AmQmatrix+AmQmatrix')/2,Bmatrix
end
=#