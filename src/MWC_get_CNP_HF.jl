# This is basically single RMG Hartree-Fock calculation. Let's  just customize it to keep the filling at charge neutrality.
using LinearAlgebra



function Coulomb(gatedis::Float64,kvec::Vector{Float64})::Float64

  if norm(kvec)==0.0
      return 9047.5636*gatedis

    else
     return 9047.5636/norm(kvec)*tanh(norm(kvec)*gatedis)
  
  end
  
end





function get_Ham(k::Vector{Float64},uD::Float64,valley::Int64,stacking::Int,NL::Int)
 
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

  if NL>2
   for layer in 1:NL-2
      Ham[2*layer-1:2*layer,2*layer+3:2*layer+4]=[0.0 t2/2;0.0 0.0]
   end
 end

  Ham=Ham+Ham'

  for layer in 1:NL
      Ham[2*layer-1:2*layer,2*layer-1:2*layer]=[uD*(layer-(NL+1)/2) -t0*get_f((k+Kac)*stacking);-t0*conj(get_f((k+Kac)*stacking)) uD*(layer-(NL+1)/2)]
  end

  
  
 return Ham
end



function get_f(k::Vector{Float64})
  delta1=1/√3*0.246*[0,1]
  delta2=1/√3*0.246*[√3/2,-1/2]
  delta3=1/√3*0.246*[-√3/2,-1/2]

  return exp(im*dot(k,delta1))+exp(im*dot(k,delta2))+exp(im*dot(k,delta3))
end






function get_single_particle(radius::Float64,num_points::Int,uD::Float64,NL::Int)

  vset=[1,-1]


  kx_grid=collect(range(-radius/2, stop=+radius/2, length=num_points))
  ky_grid=collect(range(-radius/2, stop=+radius/2, length=num_points))
  Area=4*π^2/((kx_grid[2]-kx_grid[1])*(ky_grid[2]-ky_grid[1]))

  valley_num=2
  spin_num=2
  sublattice_num=2*NL
  dimension=valley_num*spin_num*sublattice_num
  
  eig_set=[Vector{Float64}[] for _ in 1:spin_num,_ in 1:valley_num]
  k_set=Vector{Float64}[]
  k_index=Vector{Int}[]
  eig_vec_set=[Matrix{ComplexF64}[] for _ in 1:spin_num, _ in 1:valley_num]
  Ham_set=[Matrix{ComplexF64}[] for _ in 1:spin_num, _ in 1:valley_num]
      
  for ja in eachindex(kx_grid),jb in eachindex(ky_grid)
    push!(k_set,[kx_grid[ja],ky_grid[jb]])
    push!(k_index,[ja,jb])

    for vi in 1:valley_num, si in 1:spin_num
      Ham=get_Ham([kx_grid[ja],ky_grid[jb]],uD,vset[vi],1,NL)
      FFF=eigen(Ham)
      push!(eig_set[si,vi],real(FFF.values))
      push!(Ham_set[si,vi],Ham)
      push!(eig_vec_set[si,vi],FFF.vectors)

    end
  end
  
   single_matrix_complex=zeros(ComplexF64,valley_num,spin_num,sublattice_num,valley_num,spin_num,sublattice_num,length(k_set))
  

  for ja in eachindex(k_set), vi in 1:valley_num, si in 1:spin_num
    single_matrix_complex[vi,si,:,vi,si,:,ja]+=Ham_set[si,vi][ja]
  end

  single_matrix=reshape(single_matrix_complex,(dimension,dimension,length(k_set)))

  return eig_set,k_set,k_index,eig_vec_set,Area,single_matrix

end




function Construct_projector(k_set::Vector{Vector{Float64}},ϵr::Float64,
  density_matrix::Array{ComplexF64},single_matrix::Array{ComplexF64},
  energy_input::Float64,BG_density_matrix::Array{ComplexF64},fcmatrix::Array{Float64},Area::Float64,
  dimension::Int,itcount::Int,coherence::Int,NL::Int)
 
  valley_num=2
  spin_num=2
  sublattice_num=2*NL

 Fock_matrix=zeros(ComplexF64,dimension,dimension,length(k_set))
 Hartree_matrix=zeros(ComplexF64,dimension,dimension)
 HF_eigenvalues=zeros(Float64,dimension,length(k_set))
 HF_eigenvectors=zeros(ComplexF64,dimension,dimension,length(k_set))




 Threads.@threads for ja in 1:dimension
  for jb in 1:(ja - 1)
      mul!(
          @view(Fock_matrix[ja, jb, :]),
          @view(fcmatrix[:, :]),
          @view(density_matrix[ja, jb, :]),
          1,
          1,
      )
  end
 end
 

 Threads.@threads for ja in eachindex(k_set)
  axpy!(1, @view(Fock_matrix[:, :, ja])', @view(Fock_matrix[:, :, ja]))
 end


 Threads.@threads for ja in 1:dimension
  mul!(
      @view(Fock_matrix[ja, ja, :]),
      @view(fcmatrix[:, :]),
      @view(density_matrix[ja, ja, :]),
      1,
      1,
  )
 end
 Fock_matrix=Fock_matrix*1/(ϵr*Area)


 
  



 Threads.@threads for ja in eachindex(k_set)
 
    FFF=eigen(Hartree_matrix-Fock_matrix[:,:,ja]+single_matrix[:,:,ja])
    HF_eigenvectors[:,:,ja]=FFF.vectors
    HF_eigenvalues[:,ja]=real(FFF.values)
  
 
 end


  occupied_state_count = div(dimension * length(k_set), 2)

  sorted_state_indices = sortperm(vec(HF_eigenvalues))

  occupation_numbers = zeros(Float64, size(HF_eigenvalues))
  occupation_numbers[
      sorted_state_indices[1:occupied_state_count]
  ] .= 1.0

  highest_occupied_energy =
      HF_eigenvalues[sorted_state_indices[occupied_state_count]]

  lowest_unoccupied_energy =
      HF_eigenvalues[sorted_state_indices[occupied_state_count + 1]]

  fermi_level =
      (highest_occupied_energy + lowest_unoccupied_energy) / 2



  density_matrix_new =
      zeros(ComplexF64, dimension, dimension, length(k_set))

  Threads.@threads for momentum_index in eachindex(k_set)
      density_matrix_new[:, :, momentum_index] =
          HF_eigenvectors[:, :, momentum_index] *
          Diagonal(occupation_numbers[:, momentum_index]) *
          HF_eigenvectors[:, :, momentum_index]'
  end



  if coherence == 0
      coherence_mask = reshape(
          [(valley_index_1 == valley_index_2 && spin_index_1 == spin_index_2)
          for valley_index_1 in 1:valley_num, spin_index_1 in 1:spin_num,
              valley_index_2 in 1:valley_num, spin_index_2 in 1:spin_num],
          valley_num, spin_num, 1,
          valley_num, spin_num, 1,
          1
      )

      reshape(
          density_matrix_new,
          valley_num, spin_num, sublattice_num,
          valley_num, spin_num, sublattice_num,
          length(k_set)
      ) .*= coherence_mask
  end


 density_matrix_new-=BG_density_matrix 
 if itcount<15
   update_rate=0.2 
 else
    update_rate=rand() 
 end

 output_density_matrix=density_matrix_new*update_rate+density_matrix*(1-update_rate)


 DeltaMatrix=density_matrix_new-density_matrix
 eout=0.0
 for ja in eachindex(k_set)
   eout+=real(tr(DeltaMatrix[:,:,ja]*DeltaMatrix[:,:,ja]'))/length(k_set)
 end

 energy=0.0
 for ja in eachindex(k_set)
   energy+=real(tr(density_matrix[:,:,ja]*(Hartree_matrix/2-Fock_matrix[:,:,ja]/2+single_matrix[:,:,ja]))/length(k_set))
 end

 energy_change=real(energy-energy_input)



 return  eout,energy_change, output_density_matrix,DeltaMatrix,HF_eigenvalues,fermi_level,HF_eigenvectors,real(energy),Hartree_matrix,Fock_matrix


end






function get_initial_proj(k_set::Vector{Vector{Float64}},eig_vec_set::Array{Vector{Matrix{ComplexF64}}},NL::Int)
 
 valley_num=2
 spin_num=2
 sublattice_num=2*NL
 dimension=valley_num*spin_num*sublattice_num

   BG_density_matrix_complex=zeros(ComplexF64,valley_num,spin_num,sublattice_num,valley_num,spin_num,sublattice_num,length(k_set))

   for ja in eachindex(k_set), bandindex in 1:NL, vi in 1:valley_num, si in 1:spin_num
    BG_density_matrix_complex[vi,si,:,vi,si,:,ja]+=eig_vec_set[si,vi][ja][:,bandindex]*(eig_vec_set[si,vi][ja][:,bandindex])'
   end

   BG_density_matrix=reshape(BG_density_matrix_complex,(dimension,dimension,length(k_set)))

 initial_density_matrix=zeros(ComplexF64,dimension,dimension,length(k_set))

 for ja in eachindex(k_set)
   A=randn(dimension,dimension)+im*randn(dimension,dimension)

   initial_density_matrix[:,:,ja]+=(A+A')*0.1
 end



 return initial_density_matrix, BG_density_matrix
end




 



function iteration(initial_density_matrix::Array{ComplexF64},BG_density_matrix::Array{ComplexF64},
                           ϵr::Float64,k_set::Vector{Vector{Float64}},single_matrix::Array{ComplexF64},Area::Float64,NL::Int,gatedis::Float64,coherence::Int)
  valley_num=2
  spin_num=2
  sublattice_num=2*NL
  dimension=valley_num*spin_num*sublattice_num
  eout=1.0
  itcount=0
  bad_count=0
  energy=0.0
  energy_change=0.0
  fermi_level=0.0
  renormalized_density=0.0
  HF_eigenvalues=zeros(Float64,dimension,length(k_set))
  HF_eigenvectors=zeros(ComplexF64,dimension,dimension,length(k_set))


  DIIS_input_density_matrix=Vector{Array{ComplexF64}}(undef,3)
  DIIS_input_DeltaMatrix=Vector{Array{ComplexF64}}(undef,3)

  input_density_matrix=initial_density_matrix

  Fock_matrix=zeros(ComplexF64,dimension,dimension,length(k_set))
  Hartree_matrix=zeros(ComplexF64,dimension,dimension)

  itcount=0
 
  fcmatrix=zeros(Float64,length(k_set),length(k_set))



  tic=time()
   
  for ja in eachindex(k_set), jb in eachindex(k_set)
    fcmatrix[ja,jb]=Coulomb(gatedis,k_set[ja]-k_set[jb])
  end
  toc=time()
  println("formfactorstime",toc-tic)




  while (eout>1*10^(-14)) || (bad_count<4) || (abs(energy_change)>1*10^(-8))
      if eout<1*10^(-14)
       bad_count+=1
      end
      
      tic=time()

      if (itcount>60 && abs(eout)>10^(-2)) || (itcount>50 && abs(eout)<10^(-7))
      
        dmk=implement_DIIS(DIIS_input_density_matrix,DIIS_input_DeltaMatrix,k_set)
        if dmk==0
            itcount=0
            dmk=zeros(ComplexF64,dimension,dimension,length(k_set))
            for ja in eachindex(k_set)
              A=randn(dimension,dimension)+im*randn(dimension,dimension)
              dmk[:,:,ja]+=(A+A')*0.01
            end
        end
       

        eout,energy_change,output_density_matrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalues,fermi_level,HF_eigenvectors,energy,Hartree_matrix,Fock_matrix=Construct_projector(k_set,ϵr,
                                                                                                                                                              dmk,single_matrix,
                                                                                                                                                              energy,BG_density_matrix,
                                                                                                                                                              fcmatrix,Area,dimension,itcount,coherence,NL)
        DIIS_input_density_matrix[mod(itcount,3)+1]=dmk
        input_density_matrix=output_density_matrix
        println("using DIIS")
       
      else
  
        eout,energy_change,output_density_matrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalues,fermi_level,HF_eigenvectors,energy,Hartree_matrix,Fock_matrix=Construct_projector(k_set,ϵr,
                                                                                                                                                                                          input_density_matrix,single_matrix,
                                                                                                                                                                                           energy,BG_density_matrix,
                                                                                                                                                                                           fcmatrix,Area,dimension,itcount,coherence,NL)
                                                                                                                                                                                           

        DIIS_input_density_matrix[mod(itcount,3)+1]=input_density_matrix
        input_density_matrix=output_density_matrix
         
     

      end

    



      itcount+=1
     
      toc=time()
      println(toc-tic,"eout=$eout","energy_change=$energy_change","itcount=$itcount")
      flush(stdout)
     
    
  end
 

    if coherence == 0
      coherence_mask = reshape(
          [(valley_index_1 == valley_index_2 && spin_index_1 == spin_index_2)
          for valley_index_1 in 1:valley_num, spin_index_1 in 1:spin_num,
              valley_index_2 in 1:valley_num, spin_index_2 in 1:spin_num],
          valley_num, spin_num, 1,
          valley_num, spin_num, 1,
          1
      )
      for ja in eachindex(DIIS_input_density_matrix)
        reshape(
            DIIS_input_density_matrix[ja],
            valley_num, spin_num, sublattice_num,
            valley_num, spin_num, sublattice_num,
            length(k_set)
        ) .*= coherence_mask
      end

       reshape(
            Fock_matrix,
            valley_num, spin_num, sublattice_num,
            valley_num, spin_num, sublattice_num,
            length(k_set)
        ) .*= coherence_mask

      hartree_mask = reshape(
          [(valley_index_1 == valley_index_2 && spin_index_1 == spin_index_2)
          for valley_index_1 in 1:valley_num, spin_index_1 in 1:spin_num,
              valley_index_2 in 1:valley_num, spin_index_2 in 1:spin_num],
          valley_num, spin_num, 1,
          valley_num, spin_num, 1
      )

       reshape(
            Hartree_matrix,
            valley_num, spin_num, sublattice_num,
            valley_num, spin_num, sublattice_num
        ) .*=  hartree_mask
  end
 


  
  return HF_eigenvalues,HF_eigenvectors,energy, DIIS_input_density_matrix,fermi_level,Hartree_matrix,Fock_matrix,eout


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

      inB=safe_inverse(Bmatrix)
      if inB≠0
         coeff=inB*[0;0;0;1]
         dmk=coeff[1]*(DIIS_input_projector[1]+DIIS_input_DeltaMatrix[1])+coeff[2]*(DIIS_input_projector[2]+DIIS_input_DeltaMatrix[2])+coeff[3]*(DIIS_input_projector[3]+DIIS_input_DeltaMatrix[3])
         return dmk
      else
        return 0
      end
end


function safe_inverse(A)
  try
      return inv(A)  # Attempt to compute inverse
  catch e
      if isa(e, SingularException)
          println("Matrix is singular, doing randomstart again.")
          return 0  # Use pseudoinverse as an alternative
      else
          rethrow(e)  # If another error occurs, propagate it
      end
  end
end