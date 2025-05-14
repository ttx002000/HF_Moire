
using LinearAlgebra



function Coulomb(kvec::Vector{Float64},dis::Float64)::Float64
  
  if norm(kvec)==0.0
      return 9047.5636*(-dis)

    else
     return 9047.5636/norm(kvec)*exp(-norm(kvec)*dis)

  end
  
end

function find_FL(quasi_particle_energy::Vector{Float64},target_density::Float64,val_s::Float64,val_e::Float64,temp::Float64,Area::Float64,bg_particle_density::Float64)
  try_FL=(val_s+val_e)/2

  
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




function get_f(k::Vector{Float64})
  delta1=1/√3*0.246*[0,1]
  delta2=1/√3*0.246*[√3/2,-1/2]
  delta3=1/√3*0.246*[-√3/2,-1/2]

  return exp(im*dot(k,delta1))+exp(im*dot(k,delta2))+exp(im*dot(k,delta3))
end



function Hamiltonian(k::Vector{Float64},uD::Float64,valley::Int64,stacking::Int,NL::Int)
 
  Ham=zeros(ComplexF64,2*NL,2*NL)
  Kac=4π/(3*0.246)*[1,0]*valley
  t0=3100
  t1=380
  t2=-21
  t3=290
  t4=141
  ff=(3)^(1/2)/2*0.246*(-valley*k[1]+im*k[2])
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





function get_single_particle(radius::Float64,num_points::Int,uD::Float64,NL::Int,ildis::Float64,shift::Int)

  vset=[1,-1]
  offset=[0.0,-(NL-1)*uD+CNP]

  kx_grid=collect(range(-radius/2, stop=+radius/2, length=num_points))
  ky_grid=collect(range(-radius/2, stop=+radius/2, length=num_points))
  Area=4*π^2/((kx_grid[2]-kx_grid[1])*(ky_grid[2]-ky_grid[1]))

  valley_num=2
  spin_num=2
  layer_num=2
  shiftvector=shift*[kx_grid[2]-kx_grid[1],0.0]

  sublattice_num=Int(2*NL)
  
  eig_set=[Vector{Float64}[] for  _ in 1:valley_num, _ in 1:layer_num]
  k_set=Vector{Float64}[]
  k_index=Vector{Int}[]
  eig_vec_set=[Matrix{ComplexF64}[] for  _ in 1:valley_num, _ in 1:layer_num]
  Ham_set=[Matrix{ComplexF64}[] for  _ in 1:valley_num, _ in 1:layer_num]
      
  for ja in eachindex(kx_grid),jb in eachindex(ky_grid)
    push!(k_set,[kx_grid[ja],ky_grid[jb]])
    push!(k_index,[ja,jb])

    for vi in 1:valley_num,  li in 1:layer_num
      Ham=Hamiltonian([kx_grid[ja],ky_grid[jb]]+shiftvector*vset[vi],uD,vset[vi],1,NL)+Matrix{Float64}(I,sublattice_num,sublattice_num)*offset[li]
      FFF=eigen(Ham)
      push!(eig_set[vi,li],real(FFF.values))
      push!(Ham_set[vi,li],Ham)
      push!(eig_vec_set[vi,li],FFF.vectors)

    end
  end
  

  bandindex=[NL,NL+1]
  #Fock_formfactors=[zeros(ComplexF64,length(k_set),length(k_set)) for _ in 1:spin_num, _ in 1:valley_num, _ in 1:layer_num, _ in 1:spin_num, _ in 1:valley_num, _ in 1:layer_num]
  Fock_formfactors=Matrix{Matrix{ComplexF4}}(undef,spin_num,valley_num,layer_num)
  Hartree_formfactors=[zeros(ComplexF64,length(k_set),length(k_set)) for _ in 1:spin_num, _ in 1:valley_num, _ in 1:layer_num, _ in 1:spin_num, _ in 1:valley_num, _ in 1:layer_num]
  
  z_pos=zeros(Float64,layer_num,sublattice_num)
  z_pos[1,:]=0.335*[i for i in 0:NL-1 for _ in 1:2]
  z_pos[2,:]=0.335*[i for i in -NL+1:0 for _ in 1:2].-ildis

  index=reshape(collect(1:1:valley_num*spin_num*layer_num),spin_num,valley_num,layer_num)



  
  single_matrix=zeros(ComplexF64,spin_num,valley_num,layer_num,spin_num,valley_num,layer_num,length(k_set))

  
   
    for ja in eachindex(k_set),jb in eachindex(k_set), lione in 1:layer_num, litwo in 1:layer_num
      Vmatrix=zeros(ComplexF64,sublattice_num,sublattice_num)

      for subone in 1:sublattice_num, subtwo in 1:sublattice_num
        Vmatrix[subone,subtwo]=Coulomb(k_set[ja]-k_set[jb],abs(z_pos[lione,subone]-z_pos[litwo,subtwo]))
      end
     for sione in 1:spin_num, vione in 1:valley_num, sitwo in 1:spin_num, vitwo in 1:valley_num
    
       if index[sione,vione,lione]>=index[sitwo,vitwo,litwo]
        Fock_formfactors[sione,vione,sitwo,vitwo,litwo]=zeros(ComplexF64,length(k_set),length(k_set))
        f1=transpose(conj.(eig_vec_set[vione,lione][ja][:,bandindex[lione]]).*eig_vec_set[vione,lione][jb][:,bandindex[lione]])
        f2=conj.(eig_vec_set[vitwo,litwo][jb][:,bandindex[litwo]]).*eig_vec_set[vitwo,litwo][ja][:,bandindex[litwo]]
        Fock_formfactors[sione,vione,lione,sitwo,vitwo,litwo][ja,jb]+=f1*Vmatrix*f2
       end
      
    end
    end


    for  lione in 1:layer_num, litwo in 1:layer_num
      Vmatrix=zeros(ComplexF64,sublattice_num,sublattice_num)

     for subone in 1:sublattice_num, subtwo in 1:sublattice_num
        Vmatrix[subone,subtwo]=Coulomb([0.0,0.0],abs(z_pos[lione,subone]-z_pos[litwo,subtwo]))
     end
     for sione in 1:spin_num, vione in 1:valley_num, sitwo in 1:spin_num, vitwo in 1:valley_num,ja in eachindex(k_set),jb in eachindex(k_set)
   
  
       
      h1=transpose(conj.(eig_vec_set[vione,lione][ja][:,bandindex[lione]]).*eig_vec_set[vione,lione][ja][:,bandindex[lione]])
      h2=conj.(eig_vec_set[vitwo,litwo][jb][:,bandindex[litwo]]).*eig_vec_set[vitwo,litwo][jb][:,bandindex[litwo]]


      Hartree_formfactors[sione,vione,lione,sitwo,vitwo,litwo][ja,jb]+=h1*Vmatrix*h2
    end
    end


    Hunds_formfactors=zeros(ComplexF64,valley_num,layer_num,length(k_set),length(k_set))
    for ja in eachindex(k_set),jb in eachindex(k_set)
    for vi in 1:valley_num, li in 1:layer_num
      Hunds_formfactors[vi,li,ja,jb]=eig_vec_set[vi,li][ja][:,bandindex[li]]'*eig_vec_set[vi,li][jb][:,bandindex[li]]
    end
    end
 
 


  Fock_formfactors=reshape(Fock_formfactors,spin_num*valley_num*layer_num,spin_num*valley_num*layer_num)
  Hartree_formfactors=reshape(Hartree_formfactors,spin_num*valley_num*layer_num,spin_num*valley_num*layer_num)




  for sione in 1:spin_num, vione in 1:valley_num, lione in 1:layer_num, ja in eachindex(k_set)
     single_matrix[sione,vione,lione,sione,vione,lione,ja]+=eig_set[vione,lione][ja][bandindex[lione]]
  end
  
 single_matrix=reshape(single_matrix,spin_num*valley_num*layer_num,spin_num*valley_num*layer_num,length(k_set))


  return eig_set,k_set,k_index,eig_vec_set,Area,Fock_formfactors,Hartree_formfactors,Hunds_formfactors,single_matrix

end




function Construct_projector(k_set::Vector{Vector{Float64}},ϵr::Float64,JH::Float64,
  density_matrix::Array{ComplexF64},single_matrix::Array{ComplexF64},
  energy_input::Float64,BG_density_matrix::Array{ComplexF64},Area::Float64
  ,target_density::Float64,temp::Float64,itcount::Int,Fock_formfactors::Matrix{Matrix{ComplexF64}}
  ,Hartree_formfactors::Matrix{Matrix{ComplexF64}},Hunds_formfactors::Array{ComplexF64})
 
  spin_num=2
  valley_num=2
  layer_num=2
  dimension=spin_num*valley_num*layer_num



 Fock_matrix=zeros(ComplexF64,dimension,dimension,length(k_set))
 Hartree_matrix=zeros(ComplexF64,dimension,dimension,length(k_set))


 Hunds_Fock_matrix=zeros(ComplexF64,spin_num,valley_num,layer_num,spin_num,valley_num,layer_num,length(k_set))
 Hunds_Hartree_matrix=zeros(ComplexF64,spin_num,valley_num,layer_num,spin_num,valley_num,layer_num)

 HF_eigenvalues=zeros(Float64,dimension,length(k_set))
 HF_eigenvectors=zeros(ComplexF64,dimension,dimension,length(k_set))


   
  Threads.@threads for ja in 1:dimension
  for jb in 1:(ja - 1)
      mul!(
          @view(Fock_matrix[ja, jb, :]),
          @view(Fock_formfactors[ja, jb][:, :]),
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
      @view(Fock_formfactors[ja, ja][:, :]),
      @view(density_matrix[ja, ja, :]),
      1,
      1,
  )
 end

 Fock_matrix=Fock_matrix*1/(ϵr*Area)
 for ja in 1:dimension, jb in 1:dimension
     Hartree_matrix[ja,ja,:]+=Hartree_formfactors[ja,jb]*density_matrix[jb,jb,:]*1/(ϵr*Area)
 end



  density_matrix=reshape(density_matrix,(spin_num,valley_num,layer_num,spin_num,valley_num,layer_num,length(k_set)))
 paulimatrix=[[0 1;1 0],[0 -im;im 0],[1 0;0 -1]]

 density_matrix_transformed=zeros(ComplexF64,spin_num,valley_num,layer_num,spin_num,valley_num,layer_num,length(k_set))
  for sione in 1:spin_num, sitwo in 1:spin_num, pp in 1:3, sithree in 1:spin_num, sifour in 1:spin_num
     density_matrix_transformed[sione,:,:,sitwo,:,:,:]+=paulimatrix[pp][sione,sithree]*density_matrix[sithree,:,:,sifour,:,:,:]*paulimatrix[pp][sifour,sitwo]
  end



  for  vione in 1:valley_num, vitwo in 1:valley_num, li in 1:layer_num
    if vione≠vitwo
       ff=Hunds_formfactors[vione,li,:,:].*conj.(Hunds_formfactors[vitwo,li,:,:])
        for sione in 1:spin_num, sitwo in 1:spin_num
          Hunds_Fock_matrix[sione,vione,li,sitwo,vitwo,li,:]+=(ff)*density_matrix_transformed[sione,vione,li,sitwo,vitwo,li,:]*JH/Area
      end
    end
  end


  for vione in 1:valley_num
    aa=dropdims(sum(density_matrix,dims=7),dims=7)
    if vione==1
      vitwo=2
    else
      vitwo=1
    end

    for sione in 1:spin_num, sitwo in 1:spin_num, pp in 1:3, soneprime in 1:spin_num, stwoprime in 1:spin_num, li in 1:layer_num
      Hunds_Hartree_matrix[sione,vione,li,sitwo,vione,li]+=aa[stwoprime,vitwo,li,soneprime,vitwo,li]*paulimatrix[pp][soneprime,stwoprime]*paulimatrix[pp][sione,sitwo]*JH/Area
    end
  end


   Hunds_Hartree_matrix=reshape(Hunds_Hartree_matrix,(spin_num*valley_num*layer_num,spin_num*valley_num*layer_num))
   Hunds_Fock_matrix=reshape(Hunds_Fock_matrix,(spin_num*valley_num*layer_num,spin_num*valley_num*layer_num,length(k_set)))

  for ja in eachindex(k_set)
    Hartree_matrix[:,:,ja]+=Hunds_Hartree_matrix
    Fock_matrix[:,:,ja]+=Hunds_Fock_matrix[:,:,ja]
  end

 density_matrix=reshape(density_matrix,(spin_num*valley_num*layer_num,spin_num*valley_num*layer_num,length(k_set)))



 Threads.@threads for ja in eachindex(k_set)
    totalmatrix=Hartree_matrix[:,:,ja]-Fock_matrix[:,:,ja]+single_matrix[:,:,ja]
 
    FFF=eigen(totalmatrix)
    HF_eigenvectors[:,:,ja]=FFF.vectors
    HF_eigenvalues[:,ja]=real(FFF.values)
  
 
 end

  

 val_s=sort(vec(HF_eigenvalues))[1]
 val_e=sort(vec(HF_eigenvalues))[end]

  bg_particle_density=length(k_set)*valley_num*spin_num/Area




 fermi_level,renormalized_density=find_FL(vec(HF_eigenvalues),target_density,val_s,val_e,temp,Area,bg_particle_density)



 density_matrix_new=zeros(ComplexF64,dimension,dimension,length(k_set))
 
 for jb in 1:valley_num*spin_num*layer_num
   Threads.@threads for ja in eachindex(k_set) 
      density_matrix_new[:,:,ja]+=HF_eigenvectors[:,jb,ja]*(HF_eigenvectors[:,jb,ja])'*1/(exp((HF_eigenvalues[jb,ja]-fermi_level)/temp)+1)
    end
 end


  density_matrix_new=reshape(density_matrix_new,spin_num,valley_num,layer_num,spin_num,valley_num,layer_num,length(k_set))
  density_matrix_new[1,:,:,2,:,:,:].=0.0
  density_matrix_new[2,:,:,1,:,:,:].=0.0
   density_matrix_new=reshape(density_matrix_new,dimension,dimension,length(k_set))

 density_matrix_new-=BG_density_matrix 
 if itcount<15
  update_rate=0.2
 else
    update_rate=rand() #for debug
 end

 output_density_matrix=density_matrix_new*update_rate+density_matrix*(1-update_rate)


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
 


 return  eout,energy_change, output_density_matrix,DeltaMatrix,HF_eigenvalues,fermi_level,HF_eigenvectors,real(energy),Hartree_matrix,Fock_matrix,renormalized_density


end






function get_initial_proj(k_set::Vector{Vector{Float64}})
 
 valley_num=2
 spin_num=2
 layer_num=2
 
   BG_density_matrix=zeros(ComplexF64,spin_num*valley_num,layer_num,spin_num*valley_num,layer_num,length(k_set))

    for ja in eachindex(k_set)
      BG_density_matrix[:,1,:,1,ja]=Matrix{ComplexF64}(I,valley_num*spin_num,valley_num*spin_num)
    end
  
    BG_density_matrix=reshape(BG_density_matrix,spin_num*valley_num*layer_num,spin_num*valley_num*layer_num,length(k_set))
  
  

 initial_density_matrix=zeros(ComplexF64,valley_num*spin_num*layer_num,valley_num*spin_num*layer_num,length(k_set))

 for ja in eachindex(k_set)
   A=randn(valley_num*spin_num*layer_num,valley_num*spin_num*layer_num)+im*randn(valley_num*spin_num*layer_num,valley_num*spin_num*layer_num)

   initial_density_matrix[:,:,ja]+=(A+A')*1.0
 end

 return initial_density_matrix, BG_density_matrix
end



function iteration(initial_density_matrix::Array{ComplexF64},BG_density_matrix::Array{ComplexF64},
                           ϵr::Float64,k_set::Vector{Vector{Float64}},single_matrix::Array{ComplexF64},Area::Float64,
                           target_density::Float64,temp::Float64,Fock_formfactors::Matrix{Matrix{ComplexF64}},Hartree_formfactors::Matrix{Matrix{ComplexF64}}
                           ,Hunds_formfactors::Array{ComplexF64},JH::Float64,pairing::Int64)

  valley_num=2
  spin_num=2
  layer_num=2
  dimension=valley_num*spin_num*layer_num

  eout=1.0
  itcount=0
  bad_count=0
  energy=0.0
  energy_change=0.0
  fermi_level=0.0
  renormalized_density=0.0
  HF_eigenvalues=zeros(Float64,valley_num*spin_num,length(k_set))
  HF_eigenvectors=zeros(ComplexF64,valley_num*spin_num,valley_num*spin_num,length(k_set))


  DIIS_input_density_matrix=Vector{Array{ComplexF64}}(undef,3)
  DIIS_input_DeltaMatrix=Vector{Array{ComplexF64}}(undef,3)

  input_density_matrix=initial_density_matrix

  Fock_matrix=zeros(ComplexF64,valley_num*spin_num*layer_num,valley_num*spin_num*layer_num,length(k_set))
  Hartree_matrix=zeros(ComplexF64,valley_num*spin_num*layer_num,valley_num*spin_num*layer_num,length(k_set))

  itcount=0
  








  function process_pairing(density_matrix_new::Array{ComplexF64},pairing::Int)
    density_matrix_new=reshape(density_matrix_new,spin_num,valley_num,layer_num,spin_num,valley_num,layer_num,length(k_set))
    
   if pairing==-1 #no interlayer coherence
      density_matrix_new[:,:,1,:,:,2,:].=0.0
      density_matrix_new[:,:,2,:,:,1,:].=0.0
   end
   density_matrix_new=reshape(density_matrix_new,dimension,dimension,length(k_set))
   
   return density_matrix_new

  end

  while (eout>1*10^(-14)) || (bad_count<4) || (energy_change>1*10^(-6))
      if eout<1*10^(-14)
       bad_count+=1
      end
      
    
   
      tic=time()

      if (itcount>60 && abs(eout)>10^(-2)) || (itcount>70 && abs(eout)<10^(-7))
      
        dmk=implement_DIIS(DIIS_input_density_matrix,DIIS_input_DeltaMatrix,k_set)
        if dmk==0
            itcount=0
            dmk=zeros(ComplexF64,valley_num*spin_num,valley_num*spin_num,length(k_set))
            for ja in eachindex(k_set)
              A=randn(valley_num*spin_num,valley_num*spin_num)+im*randn(valley_num*spin_num,valley_num*spin_num)
              dmk[:,:,ja]+=(A+A')*0.01
            end
        end
       
       dmk=process_pairing(dmk,pairing)
        eout,energy_change,output_density_matrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalues,fermi_level,HF_eigenvectors,energy,Hartree_matrix,Fock_matrix,renormalized_density=Construct_projector(k_set,ϵr,JH,
                                                                                                                                                                                                         dmk,single_matrix,
                                                                                                                                                                                                         energy,BG_density_matrix,Area,
                                                                                                                                                                                                        target_density,temp,itcount,Fock_formfactors,Hartree_formfactors,Hunds_formfactors)
        DIIS_input_density_matrix[mod(itcount,3)+1]=dmk
        input_density_matrix=output_density_matrix
        println("using DIIS")
       
      else
        input_density_matrix=process_pairing(input_density_matrix,pairing)
  
        eout,energy_change,output_density_matrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalues,fermi_level,HF_eigenvectors,energy,Hartree_matrix,Fock_matrix,renormalized_density=Construct_projector(k_set,ϵr,JH,
                                                                                                                                                                                                                input_density_matrix,single_matrix,
                                                                                                                                                                                                               energy,BG_density_matrix,
                                                                                                                                                                                                                Area,target_density,temp,itcount,Fock_formfactors,Hartree_formfactors,Hunds_formfactors)
                                                                                                                                                                                           

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