
using LinearAlgebra



function Coulomb(dis::Float64,kvec::Vector{Float64})::Float64
  
  if norm(kvec)==0.0
     return 9047.5636*(-dis)
  else
     return 9047.5636/norm(kvec)*exp(-norm(kvec)*dis)
  end
  



end




function Coulomb_matrix(z_pos::Vector{Float64},kvec::Vector{Float64})::Matrix{Float64}
 
  
  cmatrix=kron(ones(4,4),[Coulomb(abs(z_pos[ja]-z_pos[jb]),kvec) for ja in 1:10, jb in 1:10])
  
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

function get_single_particle(radius::Float64,num_points::Int,uD::Float64)

  vone=1
  vtwo=-1

  kx_grid=collect(range(-radius/2, stop=+radius/2, length=num_points))
  ky_grid=collect(range(-radius/2, stop=+radius/2, length=num_points))
  Area=4*π^2/((kx_grid[2]-kx_grid[1])*(ky_grid[2]-ky_grid[1]))
  
  eig_set=[Vector{Float64}[] for _ in 1:4]
  k_set=Vector{Float64}[]
  k_index=Vector{Int}[]
  Ham_set=[Matrix{ComplexF64}[] for _ in 1:4]
  eig_vec_set=[Matrix{ComplexF64}[] for _ in 1:4]
  
      
  for ja in eachindex(kx_grid),jb in eachindex(ky_grid)
      
      Ham=Hamiltonian([kx_grid[ja],ky_grid[jb]],uD,vone,1)
      FFF=eigen(Ham)
      push!(eig_set[1],real(FFF.values))
      push!(eig_set[2],real(FFF.values))
      push!(k_set,[kx_grid[ja],ky_grid[jb]])
      push!(k_index,[ja,jb])
      push!(Ham_set[1],Ham)
      push!(Ham_set[2],Ham)
  
      push!(eig_vec_set[1],FFF.vectors)
      push!(eig_vec_set[2],FFF.vectors)
  end
  
  for ja in eachindex(kx_grid),jb in eachindex(ky_grid)
      
      Ham=Hamiltonian([kx_grid[ja],ky_grid[jb]],uD,vtwo,1)
      FFF=eigen(Ham)
      push!(eig_set[3],real(FFF.values))
      push!(eig_set[4],real(FFF.values))
     
      push!(Ham_set[3],Ham)
      push!(Ham_set[4],Ham)
  
      push!(eig_vec_set[3],FFF.vectors)
      push!(eig_vec_set[4],FFF.vectors)
  end
  


  
  single_matrix=zeros(ComplexF64,40,40,length(k_set))
  for findex in 1:4, ja in eachindex(k_set)
    single_matrix[10*(findex-1)+1:10*(findex-1)+10,10*(findex-1)+1:10*(findex-1)+10,ja]=Ham_set[findex][ja]
  end

  return eig_set,k_set,k_index,eig_vec_set,Area,single_matrix

end




function Construct_projector(k_set::Vector{Vector{Float64}},ϵr::Float64,Area::Float64,
  density_matrix::Array{ComplexF64},single_matrix::Array{ComplexF64},
  energy_input::Float64,BG_density_matrix::Array{ComplexF64},bg_particle_density::Float64,target_density::Float64,temp::Float64)
 dimension=40
 Fock_matrix=zeros(ComplexF64,dimension,dimension,length(k_set))
 Hartree_matrix=zeros(ComplexF64,dimension,dimension)
 HF_eigenvalues=zeros(Float64,dimension,length(k_set))
 HF_eigenvectors=zeros(ComplexF64,dimension,dimension,length(k_set))
 z_pos=reduce(vcat,[0.335*[0,0,1,1,2,2,3,3,4,4] for _ in 1:4])



 Threads.@threads for ja in eachindex(k_set)
  for jb in eachindex(k_set)
  cmatrix=Coulomb_matrix(z_pos,k_set[ja]-k_set[jb])
  Fock_matrix[:,:,ja]+=1/(ϵr*Area)*(cmatrix.*density_matrix[:,:,jb])
  end
 end



  cmatrix=Coulomb_matrix(z_pos,[0.0,0.0])
  Hartree_matrix+=1/(ϵr*Area)*diagm(cmatrix*diag(dropdims(sum(density_matrix,dims=3),dims=3)))




 Threads.@threads for ja in eachindex(k_set)
  FFF=eigen(Hartree_matrix-Fock_matrix[:,:,ja]+single_matrix[:,:,ja])
  HF_eigenvectors[:,:,ja]=FFF.vectors
  HF_eigenvalues[:,ja]=real(FFF.values)
 end

  val_s=sort(vec(HF_eigenvalues))[1]
  val_e=sort(vec(HF_eigenvalues))[end]
 fermi_level,renormalized_density=find_FL(vec(HF_eigenvalues),target_density,val_s,val_e,temp,Area,bg_particle_density)

 density_matrix_new=zeros(ComplexF64,40,40,length(k_set))

 Threads.@threads for ja in eachindex(k_set)
 for jb in 1:40

   density_matrix_new[:,:,ja]+=HF_eigenvectors[:,jb,ja]*(HF_eigenvectors[:,jb,ja])'*1/(exp((HF_eigenvalues[jb,ja]-fermi_level)/temp)+1)
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
   energy+=real(tr(density_matrix[:,:,ja]*(Hartree_matrix/2-Fock_matrix[:,:,ja]/2+single_matrix[:,:,ja])))
 end

 energy_change=real(energy-energy_input)



 return  eout,energy_change, output_density_matrix,DeltaMatrix,HF_eigenvalues,fermi_level,renormalized_density,HF_eigenvectors,real(energy),Hartree_matrix,Fock_matrix


end



function find_FL(quasi_particle_energy::Vector{Float64},target_density::Float64,val_s::Float64,val_e::Float64,temp::Float64,Area::Float64,bg_particle_density::Float64)
  try_FL=(val_s+val_e)/2

  #stan=10^(-5)*target_density
  if target_density==0.0
    stan=10^(-7)
  else
    stan=abs(10^(-5)*target_density)
  end
  fermifactor=[1/(exp((quasi_particle_energy[ja]-try_FL)/temp)+1) for ja in eachindex(quasi_particle_energy)]
 
  fl=sum(fermifactor)/Area-bg_particle_density



 if abs(fl-target_density)<stan
  println("diff",fl-target_density)
    return try_FL,fl
  elseif fl-target_density>=stan
    println("diff",fl-target_density)
    return find_FL(quasi_particle_energy,target_density,val_s, try_FL,temp,Area,bg_particle_density)
  elseif fl-target_density<=-stan
    println("diff",fl-target_density)
    return find_FL(quasi_particle_energy,target_density,try_FL,val_e,temp,Area,bg_particle_density)
   end
 

end


function get_initial_proj(k_set::Vector{Vector{Float64}},eig_vec_set::Vector{Vector{Matrix{ComplexF64}}})
 
   BG_density_matrix=zeros(ComplexF64,40,40,length(k_set))
   for ja in eachindex(k_set), findex in 1:4, bandindex in 1:5
    BG_density_matrix[10*(findex-1)+1:10*(findex-1)+10,10*(findex-1)+1:10*(findex-1)+10,ja]+=eig_vec_set[findex][ja][:,bandindex]*(eig_vec_set[findex][ja][:,bandindex])'
   end

 initial_density_matrix=zeros(ComplexF64,40,40,length(k_set))

 for ja in eachindex(k_set)
   A=randn(40,40)+im*randn(40,40)
   initial_density_matrix[:,:,ja]+=(A+A')*0.1
 end

 return initial_density_matrix, BG_density_matrix
end



function iteration(initial_density_matrix::Array{ComplexF64},BG_density_matrix::Array{ComplexF64},
                           ϵr::Float64,k_set::Vector{Vector{Float64}},single_matrix::Array{ComplexF64},target_density::Float64,temp::Float64)
  dimension=40
  eout=1.0
  itcount=0
  bad_count=0
  energy=0.0
  energy_change=0.0
  fermi_level=0.0
  bg_particle_density=length(k_set)*dimension/(2*Area)
  HF_eigenvalues=zeros(Float64,dimension,length(k_set))
  HF_eigenvectors=zeros(ComplexF64,dimension,dimension,length(k_set))
  renormalized_density=0.0

  DIIS_input_density_matrix=Vector{Array{ComplexF64}}(undef,3)
  DIIS_input_DeltaMatrix=Vector{Array{ComplexF64}}(undef,3)

  input_density_matrix=initial_density_matrix

  Fock_matrix=zeros(ComplexF64,dimension,dimension,length(k_set))
  Hartree_matrix=zeros(ComplexF64,dimension,dimension)

  itcount=0




  while (eout>1*10^(-12)) || (bad_count<4) || (energy_change>1*10^(-6))
      if eout<1*10^(-12)
       bad_count+=1
      end
      
      tic=time()

      #if (itcount>30 && abs(eout)>10^(-4)) || (itcount>30 && abs(eout)<10^(-8))
      
       # dmk=implement_DIIS(DIIS_input_density_matrix,DIIS_input_DeltaMatrix,k_set)

       

        #eout,energy_change,output_density_matrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalues,fermi_level,renormalized_density,HF_eigenvectors,energy,Hartree_matrix,Fock_matrix=Construct_projector(k_set,ϵr,Area,
                                                                                                                                                              #dmk,single_matrix,
                                                                                                                                                              #energy,BG_density_matrix,
                                                                                                                                                              #bg_particle_density,target_density,temp)
        #DIIS_input_density_matrix[mod(itcount,3)+1]=dmk
        #input_density_matrix=output_density_matrix
        #println("using DIIS")
       
      #else
        eout,energy_change,output_density_matrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalues,fermi_level,renormalized_density,HF_eigenvectors,energy,Hartree_matrix,Fock_matrix=Construct_projector(k_set,ϵr,Area,
                                                                                                                                                                                          input_density_matrix,single_matrix,
                                                                                                                                                                                           energy,BG_density_matrix,
                                                                                                                                                                                           bg_particle_density,target_density,temp)
        DIIS_input_density_matrix[mod(itcount,3)+1]=input_density_matrix
        input_density_matrix=output_density_matrix
         
     

      #end

    



      itcount+=1
     
      toc=time()
      println(toc-tic,"eout=$eout","energy_change=$energy_change")
      flush(stdout)
     
    
  end
 
 


  
  return HF_eigenvalues,HF_eigenvectors,energy, DIIS_input_density_matrix,fermi_level, renormalized_density,Hartree_matrix,Fock_matrix


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
