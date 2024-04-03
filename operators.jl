
function overlap(k::Vector{Float64},q::Vector{Float64},β::Float64)::ComplexF64
   v=q[1]^2+q[2]^2+2*im*(k[1]*q[2]-k[2]*q[1])
   #v=2*im*(k[1]*q[2]-k[2]*q[1])
   return exp(-β/4*v)
end




function Coulomb(k::Vector{Int64},T1::Vector{Float64},T2::Vector{Float64})::Float64
 

  return k==[0,0] ? 0.0 : 1/norm(k[1]*T1+k[2]*T2)
  
end

function Construct_DensityMatrix(loop_dic::Dict{Vector{Int},Any},allowedq::Vector{Vector{Int}},T1::Vector{Float64},T2::Vector{Float64},Nq::Int64,wave::Vector{Vector{Int64}},input_DensityMatrix::Vector{Matrix{ComplexF64}},single_Ham::Vector{Matrix{ComplexF64}},single_MoirePo::Vector{Matrix{ComplexF64}},constq::Float64,overlapmatrix::Array{ComplexF64,4})::Tuple{Float64,Vector{Matrix{ComplexF64}},Vector{Matrix{ComplexF64}},Vector{Vector{Float64}}}
  
 
   dimension=length(wave)
  HartreeMatrix=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
  FockMatrix=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
   output_DensityMatrix=Vector{Any}(undef,Nq^2)
  DeltaMatrix=Vector{Any}(undef,Nq^2)
  NewDensityMatrix=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
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

 bound=sort(reduce(vcat,HF_eigenvalue))[Nq^2+1]

  for ja in 1:Nq^2
    
       for jd in eachindex(HF_eigenvalue[ja])
          if HF_eigenvalue[ja][jd]<bound
             NewDensityMatrix[ja]+=HF_eigenvector[ja][:,jd]*(HF_eigenvector[ja][:,jd])'
          end
       end
       DeltaMatrix[ja]=NewDensityMatrix[ja]-input_DensityMatrix[ja]
       output_DensityMatrix[ja]=0.0*input_DensityMatrix[ja]+1.0*NewDensityMatrix[ja]
  end


  
  e1=0.0
  for ja in 1:Nq^2
    e1+=tr(DeltaMatrix[ja]'*DeltaMatrix[ja])
  end
  eout=real(e1)
  
  
 

 return  eout,output_DensityMatrix,DeltaMatrix,HF_eigenvalue
end



function metric(wavelist::Vector{Vector{Int64}},β::Float64,k::Vector{Float64},q::Vector{Float64},T1::Vector{Float64},T2::Vector{Float64})::Matrix{ComplexF64}
    Amatrix=zeros(ComplexF64,length(wavelist),length(wavelist))
    for ja in 1:length(wavelist)
    Amatrix[ja,ja]=overlap(k+wavelist[ja][1]*T1+wavelist[ja][2]*T2,q,β)
    end
    return Amatrix
end



function Construct_HFmatrix(loop_dic::Dict{Vector{Int},Any},pathpointindex::Int64,pathpoint::Vector{Int64},allowedq::Vector{Vector{Int}},T1::Vector{Float64},T2::Vector{Float64},Nq::Int64,wave::Vector{Vector{Int64}},input_DensityMatrix::Vector{Matrix{ComplexF64}},constq::Float64,chern_overlapmatrix::Array{ComplexF64,4})::Matrix{ComplexF64}
  
  
  dimension=length(wave)
  HartreeMatrix=zeros(ComplexF64,dimension,dimension) 
  FockMatrix=zeros(ComplexF64,dimension,dimension) 
 
  


    
    for jk1 in 1:Nq^2
        dmk = input_DensityMatrix[jk1]
        q=allowedq[jk1]-pathpoint
       for (dg,loop_dic_dg) in loop_dic
           CoulF1=Coulomb(q+dg,T1,T2)     
           for (gg2,loop_dic_dg_gg2) in loop_dic_dg          
               CoulF=CoulF1*chern_overlapmatrix[jk1,gg2[2],pathpointindex,gg2[1]]
           for g1g3 in loop_dic_dg_gg2
               FockMatrix[g1g3[2],gg2[1]]+=dmk[g1g3[1],gg2[2]]*CoulF*chern_overlapmatrix[pathpointindex,g1g3[2],jk1,g1g3[1]]
           end 
           end
   
       end    
    end

  

  Hartree_Density=zeros(ComplexF64,dimension,dimension)
  for jk1 in 1:Nq^2
   Hartree_Density+=input_DensityMatrix[jk1] .* transpose(chern_overlapmatrix[jk1,:,jk1,:])
  end

 

    
 
    for dg in keys(loop_dic)
        CoulH1=Coulomb(dg,T1,T2)
        for gg2 in keys(loop_dic[dg])          
            CoulH=CoulH1*chern_overlapmatrix[pathpointindex,gg2[2],pathpointindex,gg2[1]]            
        for g1g3 in loop_dic[dg][gg2]           
            HartreeMatrix[gg2[2],gg2[1]]+=Hartree_Density[g1g3[1],g1g3[2]]*CoulH               
        end 
        end

    end    
 


 

 return  constq*(HartreeMatrix-FockMatrix)
end
