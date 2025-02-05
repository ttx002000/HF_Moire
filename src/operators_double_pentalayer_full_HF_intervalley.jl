
using LinearAlgebra



function Coulomb(dis::Float64,kvec::Vector{Float64})::Float64
  
  if norm(kvec)==0.0
      return 9047.5636*(-dis)
    
    else
     return 9047.5636/norm(kvec)*exp(-norm(kvec)*dis)
     
  end
  



end




function Coulomb_matrix(z_pos::Vector{Vector{Float64}},kvec::Vector{Float64})::Matrix{Float64}
 
    spin_num=2
    layer_num=2
    sublattice_num=10
    dimension=40

  s1=[Coulomb(abs(z_pos[1][ja]-z_pos[1][jb]),kvec) for ja in 1:10, jb in 1:10]
  s2=[Coulomb(abs(z_pos[1][ja]-z_pos[2][jb]),kvec) for ja in 1:10, jb in 1:10]

  cmatrix=zeros(Float64,spin_num,layer_num,sublattice_num,spin_num,layer_num,sublattice_num)
  for sindex1 in 1:2, sindex2 in 1:2
  cmatrix[sindex1,1,:,sindex2,1,:]+=s1
  cmatrix[sindex1,1,:,sindex2,2,:]+=s2
  cmatrix[sindex1,2,:,sindex2,1,:]+=s2'
  cmatrix[sindex1,2,:,sindex2,2,:]+=s1
  end
  
  return reshape(cmatrix,(dimension,dimension))

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

function get_single_particle(radius::Float64,num_points::Int,uD::Float64,CNP::Float64)

  vset=[1,-1]
  offset=[0.0,-4*uD+CNP]

  kx_grid=collect(range(-radius/2, stop=+radius/2, length=num_points))
  ky_grid=collect(range(-radius/2, stop=+radius/2, length=num_points))
  Area=4*π^2/((kx_grid[2]-kx_grid[1])*(ky_grid[2]-ky_grid[1]))

  valley_num=2
  spin_num=2
  layer_num=2
  sublattice_num=10
  dimension=40
  
  eig_set=[Vector{Float64}[] for _ in 1:spin_num, _ in 1:layer_num,_ in 1:valley_num]
  k_set=Vector{Float64}[]
  k_index=Vector{Int}[]
  eig_vec_set=[Matrix{ComplexF64}[] for _ in 1:spin_num, _ in 1:layer_num, _ in 1:valley_num]
  Ham_set=[Matrix{ComplexF64}[] for _ in 1:spin_num, _ in 1:layer_num, _ in 1:valley_num]
      
  for ja in eachindex(kx_grid),jb in eachindex(ky_grid)
    push!(k_set,[kx_grid[ja],ky_grid[jb]])
    push!(k_index,[ja,jb])

    for vi in 1:valley_num, si in 1:spin_num, li in 1:1
      Ham=Hamiltonian([kx_grid[ja],ky_grid[jb]],uD,vset[vi],1)
      FFF=eigen(Ham)
      push!(eig_set[si,li,vi],real(FFF.values).+offset[li])
      push!(Ham_set[si,li,vi],Ham+offset[li]*Matrix{Float64}(I,sublattice_num,sublattice_num))
      push!(eig_vec_set[si,li,vi],FFF.vectors)

    end
  end

  for ja in eachindex(kx_grid),jb in eachindex(ky_grid)
    

    for vi in 1:valley_num, si in 1:spin_num, li in 2:2
      Ham=Hamiltonian([kx_grid[ja],ky_grid[jb]],uD,-vset[vi],1)
      FFF=eigen(Ham)
      push!(eig_set[si,li,vi],real(FFF.values).+offset[li])
      push!(Ham_set[si,li,vi],Ham+offset[li]*Matrix{Float64}(I,sublattice_num,sublattice_num))
      push!(eig_vec_set[si,li,vi],FFF.vectors)

    end
  end
  
   single_matrix_complex=zeros(ComplexF64,spin_num,layer_num,sublattice_num,spin_num,layer_num,sublattice_num,valley_num,length(k_set))
  

  for ja in eachindex(k_set), vi in 1:valley_num, si in 1:spin_num, li in 1:layer_num
    single_matrix_complex[si,li,:,si,li,:,vi,ja]+=Ham_set[si,li,vi][ja]
  end

  single_matrix=reshape(single_matrix_complex,(dimension,dimension,valley_num,length(k_set)))

  return eig_set,k_set,k_index,eig_vec_set,Area,single_matrix

end

function find_FL(quasi_particle_energy::Vector{Float64},target_density::Float64,val_s::Float64,val_e::Float64,temp::Float64,Area::Float64,bg_particle_density::Float64)
  try_FL=(val_s+val_e)/2

  #stan=10^(-5)*target_density
  if target_density==0.0
    stan=10^(-9)
  else
    stan=abs(10^(-8)*target_density)
  end
  fermifactor=[1/(exp((quasi_particle_energy[ja]-try_FL)/temp)+1) for ja in eachindex(quasi_particle_energy)]
 
  fl=sum(fermifactor)/Area-bg_particle_density



 if abs(fl-target_density)<stan
  
    return try_FL,fl
  elseif fl-target_density>=stan
 
    return find_FL(quasi_particle_energy,target_density,val_s, try_FL,temp,Area,bg_particle_density)
  elseif fl-target_density<=-stan
 
    return find_FL(quasi_particle_energy,target_density,try_FL,val_e,temp,Area,bg_particle_density)
   end
 

end


function Construct_projector(k_set::Vector{Vector{Float64}},ϵr::Float64,
  density_matrix::Array{ComplexF64},single_matrix::Array{ComplexF64},
  energy_input::Float64,BG_density_matrix::Array{ComplexF64},fcmatrix::Array{Float64},Area::Float64,
  target_density::Float64,temp::Float64)
 
  valley_num=2
  dimension=40
 Fock_matrix=zeros(ComplexF64,dimension,dimension,valley_num,length(k_set))
 Hartree_matrix=zeros(ComplexF64,dimension,dimension)
 HF_eigenvalues=zeros(Float64,dimension,valley_num,length(k_set))
 HF_eigenvectors=zeros(ComplexF64,dimension,dimension,valley_num,length(k_set))


 #=
 
 Threads.@threads for ja in eachindex(k_set)

    for jb in 1:valley_num
        Fock_matrix[:,:,jb,ja]+=1/(ϵr*Area)*dropdims(sum(fcmatrix[:,:,ja,:].*density_matrix[:,:,jb,:],dims=3),dims=3)
    end
 end
 =#

 



 Threads.@threads for ja in 1:dimension
  for jb in 1:ja-1,vi in 1:valley_num
  Fock_matrix[ja,jb,vi,:]+=fcmatrix[ja,jb,:,:]*density_matrix[ja,jb,vi,:]
  end
 end
 Threads.@threads for ja in eachindex(k_set)
  for vi in 1:valley_num
   Fock_matrix[:,:,vi,ja]+=Fock_matrix[:,:,vi,ja]'
  end
 end

Threads.@threads for ja in 1:dimension
  for vi in 1:valley_num
    Fock_matrix[ja,ja,vi,:]+=fcmatrix[ja,ja,:,:]*density_matrix[ja,ja,vi,:]
  end
end
Fock_matrix=Fock_matrix*1/(ϵr*Area)
 



 
  Hartree_matrix+=1/(ϵr*Area)*diagm(fcmatrix[:,:,1,1]*diag(dropdims(sum(density_matrix,dims=(3,4)),dims=(3,4))))
   



 Threads.@threads for ja in eachindex(k_set)
  for jb in 1:valley_num
    FFF=eigen(Hartree_matrix-Fock_matrix[:,:,jb,ja]+single_matrix[:,:,jb,ja])
    HF_eigenvectors[:,:,jb,ja]=FFF.vectors
    HF_eigenvalues[:,jb,ja]=real(FFF.values)
  
  end
 end

  
 #fermi_level=(sort(vec(HF_eigenvalues))[Int(length(vec(HF_eigenvalues))/2)+1]+sort(vec(HF_eigenvalues))[Int(length(vec(HF_eigenvalues))/2)])/2

 val_s=sort(vec(HF_eigenvalues))[1]
 val_e=sort(vec(HF_eigenvalues))[end]
 bg_particle_density=dimension/(2*Area)*length(k_set)*valley_num

 fermi_level,renormalized_density=find_FL(vec(HF_eigenvalues),target_density,val_s,val_e,temp,Area,bg_particle_density)


 density_matrix_new=zeros(ComplexF64,dimension,dimension,valley_num,length(k_set))
 
 for jb in 1:dimension, jc in 1:valley_num
 Threads.@threads for ja in eachindex(k_set)
 
      
           density_matrix_new[:,:,jc,ja]+=HF_eigenvectors[:,jb,jc,ja]*(HF_eigenvectors[:,jb,jc,ja])'*1/(exp((HF_eigenvalues[jb,jc,ja]-fermi_level)/temp)+1)
         
    end
 end


 density_matrix_new-=BG_density_matrix 

 output_density_matrix=density_matrix_new*0.5+density_matrix*0.5


 DeltaMatrix=density_matrix_new-density_matrix
 eout=0.0
 for ja in eachindex(k_set), jb in 1:valley_num
   eout+=real(tr(DeltaMatrix[:,:,jb,ja]*DeltaMatrix[:,:,jb,ja]'))/length(k_set)
 end

 energy=0.0
 for ja in eachindex(k_set), jb in 1:valley_num
   energy+=real(tr(density_matrix[:,:,jb,ja]*(Hartree_matrix/2-Fock_matrix[:,:,jb,ja]/2+single_matrix[:,:,jb,ja]))/length(k_set))
 end

 energy_change=real(energy-energy_input)



 return  eout,energy_change, output_density_matrix,DeltaMatrix,HF_eigenvalues,fermi_level,HF_eigenvectors,real(energy),Hartree_matrix,Fock_matrix,renormalized_density


end






function get_initial_proj(k_set::Vector{Vector{Float64}},eig_vec_set::Array{Vector{Matrix{ComplexF64}}})
 dimension=40
 valley_num=2
 spin_num=2
 layer_num=2
 sublattice_num=10

   BG_density_matrix_complex=zeros(ComplexF64,spin_num,layer_num,sublattice_num,spin_num,layer_num,sublattice_num,valley_num,length(k_set))

   for ja in eachindex(k_set), bandindex in 1:5, vi in 1:2, si in 1:2, li in 1:2
    BG_density_matrix_complex[si,li,:,si,li,:,vi,ja]+=eig_vec_set[si,li,vi][ja][:,bandindex]*(eig_vec_set[si,li,vi][ja][:,bandindex])'
   end

   BG_density_matrix=reshape(BG_density_matrix_complex,(dimension,dimension,valley_num,length(k_set)))

 initial_density_matrix=zeros(ComplexF64,dimension,dimension,valley_num,length(k_set))

 for ja in eachindex(k_set), jb in 1:valley_num
   A=randn(dimension,dimension)+im*randn(dimension,dimension)

   initial_density_matrix[:,:,jb,ja]+=(A+A')*0.01
 end



 return initial_density_matrix, BG_density_matrix
end



function iteration(initial_density_matrix::Array{ComplexF64},BG_density_matrix::Array{ComplexF64},
                           ϵr::Float64,k_set::Vector{Vector{Float64}},single_matrix::Array{ComplexF64},ildis::Float64,Area::Float64)
  dimension=40
  valley_num=2
  eout=1.0
  itcount=0
  bad_count=0
  energy=0.0
  energy_change=0.0
  fermi_level=0.0
 
  HF_eigenvalues=zeros(Float64,dimension,valley_num,length(k_set))
  HF_eigenvectors=zeros(ComplexF64,dimension,dimension,valley_num,length(k_set))


  DIIS_input_density_matrix=Vector{Array{ComplexF64}}(undef,3)
  DIIS_input_DeltaMatrix=Vector{Array{ComplexF64}}(undef,3)

  input_density_matrix=initial_density_matrix

  Fock_matrix=zeros(ComplexF64,dimension,dimension,valley_num,length(k_set))
  Hartree_matrix=zeros(ComplexF64,dimension,dimension)

  itcount=0


  z_pos=[0.335*[0,0,1,1,2,2,3,3,4,4],0.335*[-4,-4,-3,-3,-2,-2,-1,-1,0,0].-ildis]
  fcmatrix=zeros(Float64,dimension,dimension,length(k_set),length(k_set))
  for ja in eachindex(k_set), jb in eachindex(k_set)
   fcmatrix[:,:,ja,jb]+=Coulomb_matrix(z_pos,k_set[ja]-k_set[jb]) #need fix
  end



  while (eout>1*10^(-12)) || (bad_count<4) || (energy_change>1*10^(-6))
      if eout<1*10^(-12)
       bad_count+=1
      end
      
      tic=time()

      if (itcount>60 && abs(eout)>10^(-3)) || (itcount>50 && abs(eout)<10^(-7))
      
        dmk=implement_DIIS(DIIS_input_density_matrix,DIIS_input_DeltaMatrix,k_set)

       

        eout,energy_change,output_density_matrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalues,fermi_level,HF_eigenvectors,energy,Hartree_matrix,Fock_matrix,renormalized_density=Construct_projector(k_set,ϵr,
                                                                                                                                                              dmk,single_matrix,
                                                                                                                                                              energy,BG_density_matrix,
                                                                                                                                                              fcmatrix,Area,target_density,temp)
        DIIS_input_density_matrix[mod(itcount,3)+1]=dmk
        input_density_matrix=output_density_matrix
        println("using DIIS")
       
      else
  
        eout,energy_change,output_density_matrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalues,fermi_level,HF_eigenvectors,energy,Hartree_matrix,Fock_matrix,renormalized_density=Construct_projector(k_set,ϵr,
                                                                                                                                                                                          input_density_matrix,single_matrix,
                                                                                                                                                                                           energy,BG_density_matrix,
                                                                                                                                                                                           fcmatrix,Area,target_density,temp)
                                                                                                                                                                                           

        DIIS_input_density_matrix[mod(itcount,3)+1]=input_density_matrix
        input_density_matrix=output_density_matrix
         
     

      end

    



      itcount+=1
     
      toc=time()
      println(toc-tic,"eout=$eout","energy_change=$energy_change","itcount=$itcount")
      flush(stdout)
     
    
  end
 
 


  
  return HF_eigenvalues,HF_eigenvectors,energy, DIIS_input_density_matrix,fermi_level,Hartree_matrix,Fock_matrix,eout,renormalized_density


end




function implement_DIIS(DIIS_input_projector::Vector{Array{ComplexF64}},DIIS_input_DeltaMatrix::Vector{Array{ComplexF64}},k_set::Vector{Vector{Float64}})

  valley_num=2

      Bmatrix=zeros(ComplexF64,4,4)
      for ja in 1:3
       Bmatrix[ja,4]=1
       Bmatrix[4,ja]=1
      end
  
      for ja in 1:3,jb in 1:3
          for jc in eachindex(k_set), jd in 1:valley_num
             Bmatrix[ja,jb]+=real(tr((DIIS_input_DeltaMatrix[ja][:,:,jd,jc])'*(DIIS_input_DeltaMatrix[jb][:,:,jd,jc])))
          end
      end
      coeff=inv(Bmatrix)*[0;0;0;1]
     
      dmk=coeff[1]*(DIIS_input_projector[1]+DIIS_input_DeltaMatrix[1])+coeff[2]*(DIIS_input_projector[2]+DIIS_input_DeltaMatrix[2])+coeff[3]*(DIIS_input_projector[3]+DIIS_input_DeltaMatrix[3])

    return dmk
end
