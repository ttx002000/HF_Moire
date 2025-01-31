using LinearAlgebra



function get_Hamiltonian(k::Vector{Float64},tper::Float64,tpa::Float64,tNNN::Float64)
  H=zeros(ComplexF64,2,2)
  kx=k[1]
  ky=k[2]
  H[1,1]-=tpa*exp(-im*kx)+tper*exp(-im*ky)
  H[2,2]-=tpa*exp(-im*ky)+tper*exp(-im*kx)
  H[2,1]-=tNNN*exp(-im*(kx+ky))+tNNN*exp(im*(kx+ky))-tNNN*exp(-im*(kx-ky))-tNNN*exp(im*(kx-ky))
  H+=H'
  return H
end


function construct_Ham(px_xbond::Vector{Vector{Int}},px_ybond::Vector{Vector{Int}},py_xbond::Vector{Vector{Int}},py_ybond::Vector{Vector{Int}},phonon_coor::Vector{Float64},Nx::Int,Ny::Int,α::Float64,β::Float64)::Matrix{ComplexF64}
  Hphonon=zeros(ComplexF64,2*Nx*Ny,2*Nx*Ny)


  
  for ja in eachindex(px_xbond)
  Hphonon[px_xbond[ja][1],px_xbond[ja][2]]-=α*(phonon_coor[px_xbond[ja][3]]-phonon_coor[px_xbond[ja][4]])
  Hphonon[py_xbond[ja][1],py_xbond[ja][2]]-=β*(phonon_coor[py_xbond[ja][3]]-phonon_coor[py_xbond[ja][4]])
  end

  for ja in eachindex(px_ybond)
    Hphonon[px_ybond[ja][1],px_ybond[ja][2]]-=β*(phonon_coor[px_ybond[ja][3]]-phonon_coor[px_ybond[ja][4]])
    Hphonon[py_ybond[ja][1],py_ybond[ja][2]]-=α*(phonon_coor[py_ybond[ja][3]]-phonon_coor[py_ybond[ja][4]])
  end

  Hphonon+=Hphonon';
return Hphonon
end


function get_Lambdamatrix(k1::Vector{Float64},k2::Vector{Float64})

  Λmatrix=zeros(ComplexF64,4,2,2)# the index is d,o1,o2,the d direction is px_xbond, py_xbond,py_ybond, px_ybond
  
      xvec=[1,0]
      yvec=[0,1]
     Λmatrix[1,1,1]+=2*im*(sin(dot(k2,xvec))-sin(dot(k1,xvec)))
     Λmatrix[2,2,2]+=2*im*(sin(dot(k2,xvec))-sin(dot(k1,xvec)))
     Λmatrix[3,2,2]+=2*im*(sin(dot(k2,yvec))-sin(dot(k1,yvec)))
     Λmatrix[4,1,1]+=2*im*(sin(dot(k2,yvec))-sin(dot(k1,yvec)))

 return  Λmatrix
end


function get_Lambdamatrix_bandbasis(k1::Vector{Float64},k2::Vector{Float64},tper::Float64,tpa::Float64,tNNN::Float64)
    Λmatrix_bandbasis=zeros(ComplexF64,4,2,2) #the index order is d, n1,n2,q=k_1-k_2


  FFF=eigen(get_Hamiltonian(k1,tper,tpa,tNNN))
  vectork1=deepcopy(FFF.vectors)
  evalk1=deepcopy(FFF.values)

  FFF=eigen(get_Hamiltonian(k2,tper,tpa,tNNN))
  vectork2=deepcopy(FFF.vectors)
  evalk2=deepcopy(FFF.values)

  Λmatrix=get_Lambdamatrix(k1,k2)
  for di in 1:4
    Λmatrix_bandbasis[di,:,:]+=vectork1'*Λmatrix[di,:,:]*vectork2
  end
  return Λmatrix_bandbasis, evalk1, evalk2

end


function findFL(Nelec::Float64,spectrum::Vector{Float64},temp::Float64,val_s::Float64,val_e::Float64)
  fl=0
  stan=10^(-6)

  try_FL=(val_s+val_e)/2
  for ja in eachindex(spectrum)
    fd=1/(1+exp((spectrum[ja]- try_FL)/temp))
    fl+=real(fd)
  end


   if abs(fl-Nelec)<stan
    return try_FL
  elseif fl-Nelec>=stan
    return findFL(Nelec,spectrum,temp,val_s, try_FL)
  elseif fl-Nelec<=-stan
    return findFL(Nelec,spectrum,temp, try_FL,val_e)
   end

end



function  get_pertur_factor(ea::Float64,eb::Float64,FL::Float64,temp::Float64)::Float64
  x=((ea+eb)/2-FL)/temp
  y=(ea-eb)/temp
  f1=1/(exp((ea-FL)/temp)+1)

  f2=1/(exp((eb-FL)/temp)+1)
  
  return abs(y)>10^(-5) ? (f1-f2)/(ea-eb) : -1/(2*temp)*1/(cosh(x)+cosh(y/2))


end



function get_Keff(qvec::Vector{Float64},k1set::Vector{Vector{Float64}},Nsites::Int,FL::Float64,temp::Float64,tper::Float64,tpa::Float64,tNNN::Float64)
  Keff_momentum_redundant=zeros(ComplexF64,4,4)
  for k1vec in k1set
      k2vec=k1vec+qvec
      Λ1, evalk1, evalk2=get_Lambdamatrix_bandbasis(k1vec,k2vec,tper,tpa,tNNN)
      #Λ2, evalk3, evalk4=get_Lambdamatrix_bandbasis(k2vec,k1vec,tper,tpa,tNNN)
       
      Λ2=zeros(ComplexF64,4,2,2)
      for ja in 1:4
        Λ2[ja,:,:]=Λ1[ja,:,:]'
      end
      
      for  n1 in 1:2, n2 in 1:2
         
          e2=real(evalk2[n2])
          e1=real(evalk1[n1])
          factor=get_pertur_factor(e2,e1,FL,temp)
          
          Keff_momentum_redundant+=(Λ1[:,n1,n2]*transpose(Λ2[:,n2,n1])*factor)/Nsites
     
      end
  end

      return Keff_momentum_redundant

end



function runrunrun(α::Float64,β::Float64,Nx::Int64,Ny::Int64,
                   Nqx::Int,Nqy::Int64,temp::Float64,tper::Float64,tpa::Float64,tNNN::Float64,
                   K::Float64,KNNN::Float64,shearstrength::Float64,filling::Float64)
  gmatrix=zeros(Float64,4,2)
  gmatrix[1,1]=α
  gmatrix[2,1]=β
  
  gmatrix[3,2]=α
  gmatrix[4,2]=β

  Nsites=Nx*Ny
  
  
  kxspace=collect(0:1:Nx-1)*2*π/Nx
  kyspace=collect(0:1:Ny-1)*2*π/Ny
  k1set=Vector{Float64}[]
  for ja in eachindex(kxspace), jb in eachindex(kyspace)
    push!(k1set,[kxspace[ja],kyspace[jb]])
  end

  electron_spectrum=Float64[]
  for k1vec in k1set
  
      FFF=eigen(get_Hamiltonian(k1vec,tper,tpa,tNNN))
      push!(electron_spectrum,real.(FFF.values)[1])
      push!(electron_spectrum,real.(FFF.values)[2])
  end
  
  FL=findFL(Nsites*filling,sort(electron_spectrum),temp,sort(vec(electron_spectrum))[1],sort(vec(electron_spectrum))[end])
  




  qxspace=collect(0:1:Nqx-1)*2*π/Nqx
  qyspace=collect(0:1:Nqy-1)*2*π/Nqy
  
  Kbare_momentum_set=Matrix{Matrix{ComplexF64}}(undef,length(qxspace),length(qyspace))
  Threads.@threads for ja in eachindex(qxspace)
    for jb in eachindex(qyspace)
      Kbare_momentum_set[ja,jb]=get_Kbare(K,KNNN,shearstrength,[qxspace[ja],qyspace[jb]])
    end
  end

  Keff_set=Matrix{Matrix{ComplexF64}}(undef,length(qxspace),length(qyspace))
  

  Threads.@threads for ja in eachindex(qxspace)
    for jb in eachindex(qyspace)
    Keff_set[ja,jb]=get_Keff([qxspace[ja],qyspace[jb]],k1set,Nsites,FL,temp,tper,tpa,tNNN)
    end
  end


  return gmatrix,Keff_set, Kbare_momentum_set
end


function get_Kbare(K::Float64,KNNN::Float64,shearstrength::Float64,qvec::Vector{Float64})
  Kbare_momentum=zeros(ComplexF64,2,2)
  Kbare_momentum[:,:]+=2*KNNN*[1 -1;-1 1]*(1-cos(dot(qvec,[1,-1])))
  Kbare_momentum[:,:]+=2*KNNN*ones(Float64,2,2)*(1-cos(dot(qvec,[1,1])))
  Kbare_momentum[:,:]+=2*K*diagm([1-cos(dot(qvec,[1,0])),1-cos(dot(qvec,[0,1]))])
  Kbare_momentum[:,:]+=2*shearstrength*diagm([1-cos(dot(qvec,[0,1])),1-cos(dot(qvec,[1,0]))]) # note that the convention difference
  return Kbare_momentum
end