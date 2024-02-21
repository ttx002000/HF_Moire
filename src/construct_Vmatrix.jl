


 using Combinatorics
 using LinearAlgebra
 using Arpack
 using SparseArrays


function BuildVmatrix(Nx::Int64,Ny::Int64,wave::Vector{Vector{Int64}},eigenvector_single::Matrix{ComplexF64},allowedq::Vector{Vector{Int64}},T1::Vector{Float64},T2::Vector{Float64},g1T,g3T,dimension,Lb)
  Fmatrix=Array{ComplexF64}(undef,Nx*Ny,Nx*Ny,dimension)
  invmatrix=inv([g1T g3T])
  
  for  jc in 1:Nx*Ny, jq in 1:Nx*Ny, je in 1:dimension
      k2qmeshT=[mod(allowedq[jc][1]+allowedq[jq][1],Nx),mod(allowedq[jc][2]+allowedq[jq][2],Ny)]
      qT=allowedq[jq]+wave[je]
      gk2qT=allowedq[jc]+qT-k2qmeshT
  
      qvec=(qT[1])*T1+(qT[2])*T2
      eta=1;
   
      if isinteger(round((invmatrix*gk2qT)[1]/2,digits=10))&& isinteger(round((invmatrix*gk2qT)[2]/2,digits=10))
      eta=-1
      end
  
     Fmatrix[jc,jq,je]=exp(-Lb^2/4*norm(qvec)^2)*eta*exp(im*π/(Nx*Ny)*(qT[1]*allowedq[jc][2]-qT[2]*allowedq[jc][1]))*exp(im*π/(Nx*Ny)*(gk2qT[1]*k2qmeshT[2]-gk2qT[2]*k2qmeshT[1]))
      
   
  end
  
  
  Vmatrix=zeros(ComplexF64,Nx*Ny,Nx*Ny,Nx*Ny,Nx*Ny)


 
  
  
  for jc in 1:Nx*Ny,jd in 1:(jc-1),qmesh in 1:Nx*Ny
     
        k1mesh=[mod(allowedq[jc][1]+allowedq[qmesh][1],Nx),mod(allowedq[jc][2]+allowedq[qmesh][2],Ny)]
        k2mesh=[mod(allowedq[jd][1]-allowedq[qmesh][1],Nx),mod(allowedq[jd][2]-allowedq[qmesh][2],Ny)]
        k1mesh_pos=findfirst(item->item==k1mesh,allowedq)
        k2mesh_pos=findfirst(item->item==k2mesh,allowedq)
        if k1mesh_pos≠nothing && k2mesh_pos≠nothing 
          for qg in 1:length(wave)
              mqT=-1*allowedq[qmesh]-wave[qg]
              mqmesh=[mod(mqT[1],Nx),mod(mqT[2],Ny)]
              mqgT=mqT-mqmesh
              mqmesh_pos=findfirst(item->item==mqmesh,allowedq)
              mqg_pos=findfirst(item->item==mqgT,wave)
              qvec=(allowedq[qmesh][1]+wave[qg][1])*T1+(allowedq[qmesh][2]+wave[qg][2])*T2
              if mqg_pos≠nothing && mqmesh_pos≠nothing && norm(qvec)≠0
              Vmatrix[k1mesh_pos,k2mesh_pos,jc,jd]+=1/(Nx*Ny*Lb*norm(qvec))*Fmatrix[jc,qmesh,qg]*Fmatrix[jd,mqmesh_pos,mqg_pos]
              end
          end
       end
     
   end

   return Vmatrix

end








function Construct_MBstate(Nx::Int64,Ny::Int64,Nparticle::Int64,allowedq::Vector{Vector{Int64}})



   MB_state=collect(combinations(1:Nx*Ny,Nparticle))
   MB_state_can=Array{Any}(undef,Nx*Ny)
   MB_state_integer=Array{Any}(undef,Nx*Ny)
    
    for ja=1:Nx*Ny
     MB_state_can[ja]=Vector{Int64}[]
     MB_state_integer[ja]=Int64[]
    end
    
    for ja=1:length(MB_state)
        QN=sum(allowedq[MB_state[ja]])
        QN_M=[mod(QN[1],Nx),mod(QN[2],Ny)]
        push!(MB_state_can[QN_M[1]+1+QN_M[2]*Nx],MB_state[ja])
        push!(MB_state_integer[QN_M[1]+1+QN_M[2]*Nx],sum(2 .^ (MB_state[ja].-1)))
    end
    
 
    
    for ja=1:Nx*Ny
      sortindex=sortperm(MB_state_integer[ja])
      MB_state_integer[ja][:]=MB_state_integer[ja][sortindex]
      MB_state_can[ja][:]=MB_state_can[ja][sortindex]
    end

   return MB_state_can, MB_state_integer



end


function Cdag(op_index::Int64,state_index::Int64,sign::Int64)::Tuple{Int64,Int64}
    if (sign==0) | ((2^(op_index-1)& state_index)==2^(op_index-1))
    return 0,0
    end
 
    new_sign=sign*(-1)^count_ones((2^(op_index-1)-1) & (state_index))
    new_state_index=state_index+2^(op_index-1)
   
    return new_sign,new_state_index
end





function Cann(op_index::Int64,state_index::Int64,sign::Int64)::Tuple{Int64,Int64}
    if (sign==0) | ((2^(op_index-1)& state_index)≠2^(op_index-1))
    return 0,0
    end
 
    new_sign=sign*(-1)^count_ones((2^(op_index-1)-1) & (state_index))
    new_state_index=state_index-2^(op_index-1)
   
    return new_sign,new_state_index
end




function reducedV(Vmatrix::Array{ComplexF64},Nx::Int64,Ny::Int64)
  reduced_Vcol=ComplexF64[]
  reduced_Vcoor=Vector{Int64}[]
  
  for i in 1:Nx*Ny, j in 1:i-1, k in 1:Nx*Ny, p in 1:k-1
  if abs(Vmatrix[i,j,k,p]-Vmatrix[j,i,k,p])>10^(-10)
  push!(reduced_Vcol,2*Vmatrix[i,j,k,p]-2*Vmatrix[j,i,k,p])
  push!(reduced_Vcoor,[i,j,k,p])
  end
  end
  
  return reduced_Vcol, reduced_Vcoor
end






function Construct_Manybodymatrix(reduced_Vcol::Vector{ComplexF64},reduced_Vcoor::Vector{Vector{Int64}},state_can,state_integer,eigenvalue_single)

    

  matrix_index1=Int64[]
  matrix_index2=Int64[]
  matrix_value=ComplexF64[]
   
   
     for jb=1:length(state_can), jc=1:length(reduced_Vcol)
   
        (sign,state)=Cann(reduced_Vcoor[jc][3],state_integer[jb],1)
        (sign,state)=Cann(reduced_Vcoor[jc][4],state,sign)
        (sign,state)=Cdag(reduced_Vcoor[jc][2],state,sign)
        (sign,state)=Cdag(reduced_Vcoor[jc][1],state,sign)
   
         state_index=searchsortedfirst(state_integer,state)
         
        if (state*sign)≠0 &&  (state_integer[state_index]==state) 
         push!(matrix_index1,state_index)
         push!(matrix_index2,jb)
         push!(matrix_value,sign*reduced_Vcol[jc]/2)
        end
   
     end
     
     for jb=1:length(state_can)
      
       push!(matrix_index1,jb)
       push!(matrix_index2,jb)
       push!(matrix_value,sum(eigenvalue_single[state_can[jb]]))
   
     end
   
     MB_spectrum,ζ=eigs(sparse(matrix_index1,matrix_index2,matrix_value), nev=15,which=:SR)
     #MB_spectrum=eigvals(collect(sparse(matrix_index1,matrix_index2,matrix_value)))
      
  
   
   
   
   
   return MB_spectrum
   
  end
