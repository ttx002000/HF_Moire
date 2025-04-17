
using LinearAlgebra



function Coulomb(kvec::Vector{Float64})::Float64
  Dgate=40.0
  if norm(kvec)==0.0
      return 9047.5636*Dgate
   
    else
      q=norm(kvec)
     return 9047.5636/q*tanh(Dgate*q)
    
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


#=
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
=#


function Hamiltonian(k::Vector{Float64},uD::Float64,valley::Int64,stacking::Int,NL::Int)
 
  Ham=zeros(ComplexF64,2*NL,2*NL)
  #Kac=4π/(3*0.246)*[1,0]*valley
  t0=3100
  t1=380
  t2=-15 # There is something weird about this parameters here
  t3=290
  t4=141
  ff=(3)^(1/2)/2*0.246*(-valley*k[1]+im*k[2])*stacking
  Δ2=2.0 #I take this from the four layer paper
  δ=10.5
  for layer in 1:NL-1
     Ham[2*layer-1:2*layer,2*layer+1:2*layer+2]+=[t4*ff t3*conj(ff);t1 t4*ff]
  end

  if NL>2
   for layer in 1:NL-2
      Ham[2*layer-1:2*layer,2*layer+3:2*layer+4]+=[0.0 t2/2;0.0 0.0]
   end
 end

  Ham=Ham+Ham'

  for layer in 1:NL
      Ham[2*layer-1:2*layer,2*layer-1:2*layer]+=[uD*(layer-(NL+1)/2) -t0*ff;-t0*conj(ff) uD*(layer-(NL+1)/2)]
  end
  
  for layer in 2:NL-1
    Ham[2*layer-1:2*layer,2*layer-1:2*layer]+=[-2*Δ2 0;0 -2*Δ2]
  end
  Ham[1,1]+=-δ
  Ham[2*NL,2*NL]+=-δ

 return Ham
end

function get_single_particle(radius::Float64,num_points::Int,uD::Float64,NL::Int,whichside::Int,SOCcoef::Float64)

  vset=[1,-1]


  kx_grid=collect(range(-radius/2, stop=+radius/2, length=num_points))
  ky_grid=collect(range(-radius/2, stop=+radius/2, length=num_points))
  Area=4*π^2/((kx_grid[2]-kx_grid[1])*(ky_grid[2]-ky_grid[1]))

  valley_num=2
  spin_num=2

 
  
  eig_set=[Vector{Float64}[] for _ in 1:valley_num]
  k_set=Vector{Float64}[]
  k_index=Vector{Int}[]
  eig_vec_set=[Matrix{ComplexF64}[] for  _ in 1:valley_num]
  Ham_set=[Matrix{ComplexF64}[] for  _ in 1:valley_num]
      
  for ja in eachindex(kx_grid),jb in eachindex(ky_grid)
    push!(k_set,[kx_grid[ja],ky_grid[jb]])
    push!(k_index,[ja,jb])

    for vi in 1:valley_num
      Ham=Hamiltonian([kx_grid[ja],ky_grid[jb]],uD,vset[vi],1,NL)
      FFF=eigen(Ham)
      push!(eig_set[vi],real(FFF.values))
      push!(Ham_set[vi],Ham)
      push!(eig_vec_set[vi],FFF.vectors)

    end
  end
  


  if whichside==1
     bandindex=NL+1
  elseif whichside==2
     bandindex=NL
  end
  
  
  formfactors=zeros(ComplexF64,valley_num,length(k_set),length(k_set))
  Coulombmatrix=zeros(Float64,length(k_set),length(k_set))
  single_matrix=zeros(ComplexF64,spin_num,valley_num,spin_num,valley_num,length(k_set))
  
  for ja in eachindex(k_set),jb in eachindex(k_set)
    for vi in 1:valley_num
      formfactors[vi,ja,jb]=eig_vec_set[vi][ja][:,bandindex]'*eig_vec_set[vi][jb][:,bandindex]
    end

    Coulombmatrix[ja,jb]=Coulomb(k_set[ja]-k_set[jb])
  end
  pz=[1.0,-1.0]
  sublattice_operators=diagm(repeat([1, -1], NL))


  for sione in 1:spin_num, vione in 1:valley_num, ja in eachindex(k_set)
     single_matrix[sione,vione,sione,vione,ja]=eig_set[vione][ja][bandindex]+pz[sione]*pz[vione]*(eig_vec_set[vione][ja][:,bandindex]'*sublattice_operators*eig_vec_set[vione][ja][:,bandindex])*SOCcoef/2
  end
  
 single_matrix=reshape(single_matrix,valley_num*spin_num,valley_num*spin_num,length(k_set))


  return eig_set,k_set,k_index,eig_vec_set,Area,formfactors,Coulombmatrix,single_matrix

end




function Construct_projector(k_set::Vector{Vector{Float64}},ϵr::Float64,JH::Float64,
  density_matrix::Array{ComplexF64},single_matrix::Array{ComplexF64},
  energy_input::Float64,BG_density_matrix::Array{ComplexF64},Area::Float64
  ,target_density::Float64,temp::Float64,itcount::Int,Coulombmatrix::Matrix{Float64},formfactors::Array{ComplexF64},whichside::Int)
 
  spin_num=2
  valley_num=2

  density_matrix=reshape(density_matrix,(spin_num,valley_num,spin_num,valley_num,length(k_set)))

 Fock_matrix=zeros(ComplexF64,spin_num,valley_num,spin_num,valley_num,length(k_set))
 Hartree_matrix=zeros(ComplexF64,spin_num,valley_num,spin_num,valley_num)
 HF_eigenvalues=zeros(Float64,spin_num*valley_num,length(k_set))
 HF_eigenvectors=zeros(ComplexF64,spin_num*valley_num,spin_num*valley_num,length(k_set))


  for  vione in 1:valley_num, vitwo in 1:valley_num
    ff=formfactors[vione,:,:].*conj.(formfactors[vitwo,:,:]) .*Coulombmatrix
    for sione in 1: spin_num, sitwo in 1: spin_num
    Fock_matrix[sione,vione,sitwo,vitwo,:]+=(ff)*density_matrix[sione,vione,sitwo,vitwo,:]*1/(ϵr*Area)
    end    
  end

  aaone=reshape(dropdims(sum(density_matrix[:,:,:,:,:],dims=5),dims=5),valley_num*spin_num,valley_num*spin_num)
 for sione in 1:valley_num, vione in 1:valley_num
     Hartree_matrix[sione,vione,sione,vione]+=tr(aaone)*1/(ϵr*Area)*Coulomb([0.0,0.0])
 end




 paulimatrix=[[0 1;1 0],[0 -im;im 0],[1 0;0 -1]]

 density_matrix_transformed=zeros(ComplexF64,spin_num,valley_num,spin_num,valley_num,length(k_set))
  for sione in 1:spin_num, sitwo in 1:spin_num, pp in 1:3, sithree in 1:spin_num, sifour in 1:spin_num
     density_matrix_transformed[sione,:,sitwo,:,:]+=paulimatrix[pp][sione,sithree]*density_matrix[sithree,:,sifour,:,:]*paulimatrix[pp][sifour,sitwo]
  end



  for  vione in 1:valley_num, vitwo in 1:valley_num
    if vione≠vitwo
       ff=formfactors[vione,:,:].*conj.(formfactors[vitwo,:,:])
        for sione in 1:spin_num, sitwo in 1:spin_num
          Fock_matrix[sione,vione,sitwo,vitwo,:]+=(ff)*density_matrix_transformed[sione,vione,sitwo,vitwo,:]*JH/Area
      end
    end
  end


  for vione in 1:valley_num
    aa=dropdims(sum(density_matrix[:,:,:,:,:],dims=5),dims=5)
    if vione==1
      vitwo=2
    else
      vitwo=1
    end

    for sione in 1:spin_num, sitwo in 1:spin_num, pp in 1:3, soneprime in 1:spin_num, stwoprime in 1:spin_num
      Hartree_matrix[sione,vione,sitwo,vione]+=aa[stwoprime,vitwo,soneprime,vitwo]*paulimatrix[pp][soneprime,stwoprime]*paulimatrix[pp][sione,sitwo]*JH/Area
    end
  end


  Hartree_matrix=reshape(Hartree_matrix,(spin_num*valley_num,spin_num*valley_num))
  Fock_matrix=reshape(Fock_matrix,(spin_num*valley_num,spin_num*valley_num,length(k_set)))



 density_matrix=reshape(density_matrix,(spin_num*valley_num,spin_num*valley_num,length(k_set)))



 Threads.@threads for ja in eachindex(k_set)
    totalmatrix=Hartree_matrix-Fock_matrix[:,:,ja]+single_matrix[:,:,ja]
 
    FFF=eigen(totalmatrix)
    HF_eigenvectors[:,:,ja]=FFF.vectors
    HF_eigenvalues[:,ja]=real(FFF.values)
  
 
 end

  

 val_s=sort(vec(HF_eigenvalues))[1]
 val_e=sort(vec(HF_eigenvalues))[end]
 if whichside==1
  bg_particle_density=0.0
 elseif whichside==2
  bg_particle_density=length(k_set)*4/Area
 end



 fermi_level,renormalized_density=find_FL(vec(HF_eigenvalues),target_density,val_s,val_e,temp,Area,bg_particle_density)



 density_matrix_new=zeros(ComplexF64,valley_num*spin_num,valley_num*spin_num,length(k_set))
 
 for jb in 1:valley_num*spin_num
   Threads.@threads for ja in eachindex(k_set) 
           density_matrix_new[:,:,ja]+=HF_eigenvectors[:,jb,ja]*(HF_eigenvectors[:,jb,ja])'*1/(exp((HF_eigenvalues[jb,ja]-fermi_level)/temp)+1)
    end
 end


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
   energy+=real(tr(density_matrix[:,:,ja]*(Hartree_matrix/2-Fock_matrix[:,:,ja]/2+single_matrix[:,:,ja]))/length(k_set))
 end

 energy_change=real(energy-energy_input)
 


 return  eout,energy_change, output_density_matrix,DeltaMatrix,HF_eigenvalues,fermi_level,HF_eigenvectors,real(energy),Hartree_matrix,Fock_matrix,renormalized_density


end






function get_initial_proj(k_set::Vector{Vector{Float64}},whichside::Int)
 
 valley_num=2
 spin_num=2
 
   BG_density_matrix=zeros(ComplexF64,valley_num*spin_num,valley_num*spin_num,length(k_set))
   if whichside==2
    for ja in eachindex(k_set)
      BG_density_matrix[:,:,ja]=Matrix{ComplexF64}(I,valley_num*spin_num,valley_num*spin_num)
    end
   end
  
  

 initial_density_matrix=zeros(ComplexF64,valley_num*spin_num,valley_num*spin_num,length(k_set))

 for ja in eachindex(k_set)
   A=randn(valley_num*spin_num,valley_num*spin_num)+im*randn(valley_num*spin_num,valley_num*spin_num)

   initial_density_matrix[:,:,ja]+=(A+A')*1.0
 end

 if rand()>0.5
  for ja in eachindex(k_set)

 
    initial_density_matrix[:,:,ja]=initial_density_matrix[:,:,ja].*Matrix{Float64}(I,valley_num*spin_num,valley_num*spin_num)
  end
 
 end


 return initial_density_matrix, BG_density_matrix
end



function iteration(initial_density_matrix::Array{ComplexF64},BG_density_matrix::Array{ComplexF64},
                           ϵr::Float64,k_set::Vector{Vector{Float64}},single_matrix::Array{ComplexF64},Area::Float64,
                           target_density::Float64,temp::Float64,Coulombmatrix::Matrix{Float64},formfactors::Array{ComplexF64},JH::Float64,whichside::Int)

  valley_num=2
  spin_num=2

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

  Fock_matrix=zeros(ComplexF64,valley_num*spin_num,valley_num*spin_num,length(k_set))
  Hartree_matrix=zeros(ComplexF64,valley_num*spin_num,valley_num*spin_num)

  itcount=0


  while (eout>1*10^(-12)) || (bad_count<4) || (energy_change>1*10^(-6))
      if eout<1*10^(-12)
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
       

        eout,energy_change,output_density_matrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalues,fermi_level,HF_eigenvectors,energy,Hartree_matrix,Fock_matrix,renormalized_density=Construct_projector(k_set,ϵr,JH,
                                                                                                                                                                                                         dmk,single_matrix,
                                                                                                                                                                                                         energy,BG_density_matrix,Area,
                                                                                                                                                                                                        target_density,temp,itcount,Coulombmatrix,formfactors,whichside)
        DIIS_input_density_matrix[mod(itcount,3)+1]=dmk
        input_density_matrix=output_density_matrix
        println("using DIIS")
       
      else
  
        eout,energy_change,output_density_matrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalues,fermi_level,HF_eigenvectors,energy,Hartree_matrix,Fock_matrix,renormalized_density=Construct_projector(k_set,ϵr,JH,
                                                                                                                                                                                                                input_density_matrix,single_matrix,
                                                                                                                                                                                                               energy,BG_density_matrix,
                                                                                                                                                                                                                Area,target_density,temp,itcount,Coulombmatrix,formfactors,whichside)
                                                                                                                                                                                           

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