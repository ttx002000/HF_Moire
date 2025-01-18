
using LinearAlgebra



function Coulomb(dis::Float64,kvec::Vector{Float64})::Float64
  
  if norm(kvec)==0.0
     return 9047.5636*(-dis)
  else
     return 9047.5636/norm(kvec)*exp(-norm(kvec)*dis)
  end
  



end




function Coulomb_matrix(l1::Int,l2::Int,z_pos::Vector{Vector{Float64}},kvec::Vector{Float64})::Matrix{Float64}
 
  
  cmatrix=[Coulomb(abs(z_pos[l1][ja]-z_pos[l2][jb]),kvec) for ja in 1:10, jb in 1:10]
  
  return cmatrix

end


function get_f(k::Vector{Float64})
  delta1=1/√3*0.246*[0,1]
  delta2=1/√3*0.246*[√3/2,-1/2]
  delta3=1/√3*0.246*[-√3/2,-1/2]

  return exp(im*dot(k,delta1))+exp(im*dot(k,delta2))+exp(im*dot(k,delta3))
end



function Hamiltonian(k::Vector{Float64},uD::Float64,valley::Int64,stacking::Int)
  NL=5
  Ham=zeros(ComplexF64,2*NL,2*NL)
  Kac=4π/(3*0.246)*[1,0]*valley
  t0=3100
  t1=380
  t2=-21
  t3=290
  t4=141
  for layer in 1:NL-1
     Ham[2*layer-1:2*layer,2*layer+1:2*layer+2]=[t4*get_f((k+Kac)*stacking) t3*conj(get_f((k+Kac)*stacking));t1 t4*get_f((k+Kac)*stacking)]
  end

  for layer in 1:NL-2
      Ham[2*layer-1:2*layer,2*layer+3:2*layer+4]=[0.0 t2/2;0.0 0.0]
  end

  Ham=Ham+Ham'

  for layer in 1:NL
      Ham[2*layer-1:2*layer,2*layer-1:2*layer]=[uD*(layer-(NL+1)/2) -t0*get_f((k+Kac)*stacking);-t0*conj(get_f((k+Kac)*stacking)) uD*(layer-(NL+1)/2)]
  end
 return Ham
end

function get_single_particle(vone::Int,vtwo::Int,radius::Float64,num_points::Int,uD::Float64,stackingtwo::Int64,CNP::Float64)



  kx_grid=collect(range(-radius/2, stop=+radius/2, length=num_points))
  ky_grid=collect(range(-radius/2, stop=+radius/2, length=num_points))
  Area=4*π^2/((kx_grid[2]-kx_grid[1])*(ky_grid[2]-ky_grid[1]))
  
  eig_set=[Float64[] for _ in 1:2]
  k_set=Vector{Float64}[]
  k_index=Vector{Int}[]
  eig_vec_set=[Vector{ComplexF64}[] for _ in 1:2]
  bandindex=5
      
  for ja in eachindex(kx_grid),jb in eachindex(ky_grid)
   
      Ham=Hamiltonian([kx_grid[ja],ky_grid[jb]],uD,vone,1)
      FFF=eigen(Ham)
      push!(eig_set[1],real(FFF.values[bandindex]))
      push!(k_set,[kx_grid[ja],ky_grid[jb]])
      push!(k_index,[ja,jb])
      push!(eig_vec_set[1],FFF.vectors[:,bandindex])
  end
  
  
  bandindex=6
      
  for ja in eachindex(kx_grid),jb in eachindex(ky_grid)
      
      Ham=Hamiltonian([kx_grid[ja],ky_grid[jb]],uD,vtwo,stackingtwo)
      FFF=eigen(Ham)
      push!(eig_set[2],real(FFF.values[bandindex])+CNP)
   
      push!(eig_vec_set[2],FFF.vectors[:,bandindex])
  end


  
  single_matrix=zeros(2,2,length(k_set))
  for ja in eachindex(k_set)
    single_matrix[:,:,ja]=[eig_set[1][ja] 0.0;0.0  eig_set[2][ja]] 
  end

  return eig_set,k_set,k_index,eig_vec_set,Area,single_matrix

end

function get_formfactors(k_set::Vector{Vector{Float64}},eig_vec_set::Vector{Vector{Vector{ComplexF64}}},ildis::Float64)
 
  z_pos=Vector{Vector{Float64}}(undef,2)
  z_pos[1]=0.335*[0,0,1,1,2,2,3,3,4,4]
  z_pos[2]=0.335*[-4,-4,-3,-3,-2,-2,-1,-1,0,0].-ildis

  formfactors=zeros(ComplexF64,2,2,length(k_set),length(k_set))
  Threads.@threads for k2 in eachindex(k_set)
  for  k1 in eachindex(k_set), l2 in 1:2, l1 in 1:2
    cmatrix=Coulomb_matrix(l2,l1,z_pos,k_set[k2]-k_set[k1])
    s2=eig_vec_set[l2][k2].*conj.(eig_vec_set[l2][k1])
    s1=eig_vec_set[l1][k2].*conj.(eig_vec_set[l1][k1])

    formfactors[l2,l1,k2,k1]=s2'*cmatrix*s1
      
    end 
  end

  Hartree_formfactors=zeros(ComplexF64,2,2,length(k_set),length(k_set))
  Threads.@threads for k1 in eachindex(k_set)
  c1=[Coulomb_matrix(l1,l2,z_pos,[0.0,0.0]) for l1 in 1:2, l2 in 1:2]
  for  k2 in eachindex(k_set), l2 in 1:2, l1 in 1:2
    #cmatrix=Coulomb_matrix(l1,l2,z_pos,[0.0,0.0])
    cmatrix=c1[l1,l2]
    s2=abs.(eig_vec_set[l2][k2]).^2
    s1=abs.(eig_vec_set[l1][k1]).^2

    Hartree_formfactors[l1,l2,k1,k2]=s1'*cmatrix*s2
      
    end 
  end

  return formfactors,Hartree_formfactors

end



function Construct_projector(k_set::Vector{Vector{Float64}},ϵr::Float64,Area::Float64,
  density_matrix::Array{ComplexF64},single_matrix::Array{Float64},
  num_particle::Int,energy_input::Float64,BG_density_matrix::Array{ComplexF64},formfactors::Array{ComplexF64},Hartree_formfactors::Array{ComplexF64})

 Fock_matrix=zeros(ComplexF64,2,2,length(k_set))
 Hartree_matrix=zeros(ComplexF64,2,2,length(k_set))
 HF_eigenvalues=zeros(Float64,2,length(k_set))
 HF_eigenvectors=zeros(ComplexF64,2,2,length(k_set))


 Threads.@threads for ja in eachindex(k_set)
 for l1 in 1:2, l2 in 1:2
 Fock_matrix[l2,l1,ja]=1/(ϵr*Area)*sum(formfactors[l2,l1,ja,:].*density_matrix[l2,l1,:])
 end
 end 
 #=
 for l1 in 1:2, l2 in 1:2
 Fock_matrix[l2,l1,:]=1/(ϵr*Area)*(formfactors[l2,l1,:,:]*density_matrix[l2,l1,:])
 end
 =#
 Threads.@threads for ja in eachindex(k_set)
 for l1 in 1:2, l2 in 1:2
 Hartree_matrix[l1,l1,ja]+=1/(ϵr*Area)*sum(Hartree_formfactors[l1,l2,ja,:].*density_matrix[l2,l2,:])
 end
 end


 Threads.@threads for ja in eachindex(k_set)
 FFF=eigen(Hartree_matrix[:,:,ja]-Fock_matrix[:,:,ja]+single_matrix[:,:,ja])
 HF_eigenvectors[:,:,ja]=FFF.vectors
 HF_eigenvalues[:,ja]=real(FFF.values)
 end


 fermi_level=(sort(vec(HF_eigenvalues))[num_particle]+sort(vec(HF_eigenvalues))[num_particle+1])/2
 density_matrix_new=zeros(ComplexF64,2,2,length(k_set))

 Threads.@threads for ja in eachindex(k_set)
 for jb in 1:2
 if HF_eigenvalues[jb,ja]<fermi_level
 density_matrix_new[:,:,ja]+=HF_eigenvectors[:,jb,ja]*(HF_eigenvectors[:,jb,ja])'
 end
 end
 end


 density_matrix_new-=BG_density_matrix 

 output_density_matrix=density_matrix_new*1.0+density_matrix*0.0

 DeltaMatrix=density_matrix_new-density_matrix
 eout=0.0
 for ja in eachindex(k_set)
   eout+=real(tr(DeltaMatrix[:,:,ja]*DeltaMatrix[:,:,ja]'))/length(k_set)
 end

 energy=0.0
 for ja in eachindex(k_set)
   energy+=real(tr(density_matrix[:,:,ja]*(Hartree_matrix[:,:,ja]/2-Fock_matrix[:,:,ja]/2+single_matrix[:,:,ja]))/length(k_set))
 end

 energy_change=real(energy-energy_input)



 return  eout,energy_change, output_density_matrix,DeltaMatrix,HF_eigenvalues,fermi_level,HF_eigenvectors,real(energy),Hartree_matrix,Fock_matrix


end






function get_initial_proj(k_set::Vector{Vector{Float64}})
 
   BG_density_matrix=zeros(ComplexF64,2,2,length(k_set))
   BG_density_matrix[1,1,:]=ones(Float64,length(k_set))

 initial_density_matrix=zeros(ComplexF64,2,2,length(k_set))

 for ja in eachindex(k_set)
   A=randn(2,2)+im*randn(2,2)
   initial_density_matrix[:,:,ja]+=(A+A')*0.5
 end

 return initial_density_matrix, BG_density_matrix
end



function iteration(formfactors::Array{ComplexF64},initial_density_matrix::Array{ComplexF64},BG_density_matrix::Array{ComplexF64},ϵr::Float64,k_set::Vector{Vector{Float64}},single_matrix::Array{Float64},Hartree_formfactors::Array{ComplexF64})
  
  eout=1.0
  itcount=0
  bad_count=0
  energy=0.0
  energy_change=0.0
  fermi_level=0.0
  num_particle=length(k_set)
  HF_eigenvalues=zeros(Float64,2,length(k_set))
  HF_eigenvectors=zeros(ComplexF64,2,2,length(k_set))

  DIIS_input_density_matrix=Vector{Array{ComplexF64}}(undef,3)
  DIIS_input_DeltaMatrix=Vector{Array{ComplexF64}}(undef,3)

  input_density_matrix=initial_density_matrix

  Fock_matrix=zeros(ComplexF64,2,2,length(k_set))
  Hartree_matrix=zeros(ComplexF64,2,2)

  itcount=0




  while (eout>1*10^(-12)) || (bad_count<4) || (energy_change>1*10^(-6))
      if eout<1*10^(-12)
       bad_count+=1
      end
      
      tic=time()

      if (itcount>30 && abs(eout)>10^(-4)) || (itcount>30 && abs(eout)<10^(-8))

        dmk=implement_DIIS(DIIS_input_density_matrix,DIIS_input_DeltaMatrix,k_set)

       

        eout,energy_change,output_density_matrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalues,fermi_level,HF_eigenvectors,energy,Hartree_matrix,Fock_matrix=Construct_projector(k_set,ϵr,Area,
                                                                                                                                                              dmk,single_matrix,
                                                                                                                                                              num_particle,energy,BG_density_matrix,formfactors,Hartree_formfactors)
        DIIS_input_density_matrix[mod(itcount,3)+1]=dmk
        input_density_matrix=output_density_matrix
        println("using DIIS")
       
      else
        eout,energy_change,output_density_matrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalues,fermi_level,HF_eigenvectors,energy,Hartree_matrix,Fock_matrix=Construct_projector(k_set,ϵr,Area,
                                                                                                                                                                                          input_density_matrix,single_matrix,
                                                                                                                                                                                           num_particle,energy,BG_density_matrix,formfactors,Hartree_formfactors)
        DIIS_input_density_matrix[mod(itcount,3)+1]=input_density_matrix
        input_density_matrix=output_density_matrix
         
     

      end

    



      itcount+=1
     
      toc=time()
      println(toc-tic,"eout=$eout","energy_change=$energy_change")
      flush(stdout)
     
    
  end
 
 


  
  return HF_eigenvalues,HF_eigenvectors,energy, DIIS_input_density_matrix,fermi_level,Hartree_matrix,Fock_matrix


end




function implement_DIIS(DIIS_input_projector::Vector{Array{ComplexF64}},DIIS_input_DeltaMatrix::Vector{Array{ComplexF64}},k_set::Vector{Vector{Float64}})


      Bmatrix=zeros(ComplexF64,4,4)
      for ja in 1:3
       Bmatrix[ja,4]=1
       Bmatrix[4,ja]=1
      end
  
      for ja in 1:3,jb in 1:3
          for jc in eachindex(k_set)
             Bmatrix[ja,jb]+=real(tr((DIIS_input_DeltaMatrix[ja][:,:,jc])'*(DIIS_input_DeltaMatrix[jb][:,:,jc])))
          end
      end
      coeff=inv(Bmatrix)*[0;0;0;1]
     
      dmk=coeff[1]*(DIIS_input_projector[1]+DIIS_input_DeltaMatrix[1])+coeff[2]*(DIIS_input_projector[2]+DIIS_input_DeltaMatrix[2])+coeff[3]*(DIIS_input_projector[3]+DIIS_input_DeltaMatrix[3])

    return dmk
end
