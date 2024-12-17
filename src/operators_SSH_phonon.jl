using LinearAlgebra



function get_Hamiltonian(k::Vector{Float64},tper::Float64,tpa::Float64,tNNN::Float64)::Matrix{ComplexF64}
  H=zeros(ComplexF64,2,2)
  kx=k[1]
  ky=k[2]
  
  H[1,1]-=tpa*exp(-im*kx)+tper*exp(-im*ky)
  H[2,2]-=tpa*exp(-im*ky)+tper*exp(-im*kx)
  H[2,1]-=tNNN*exp(-im*(kx+ky))+tNNN*exp(im*(kx+ky))-tNNN*exp(-im*(kx-ky))-tNNN*exp(im*(kx-ky))
  H+=H'
  return H
end



function initialize(Nx::Int64,Ny::Int64)
  
  #The lattice constant is around 4A
  kxspace=collect(0:1:Nx-1)*2*π/Nx
  kyspace=collect(0:1:Ny-1)*2*π/Ny
  
  px_xbond=Vector{Vector{Int}}[]
  px_ybond=Vector{Vector{Int}}[]
  py_xbond=Vector{Vector{Int}}[]
  py_ybond=Vector{Vector{Int}}[]

  shear_xbond=Vector{Vector{Int}}[]
  shear_ybond=Vector{Vector{Int}}[]
     
  NNN_sp_d1=Vector{Vector{Int}}[]
  NNN_sp_d2=Vector{Vector{Int}}[]
  
  
  for ja in 1:Nx,jb in 1:Ny
     push!(px_xbond,[[mod(ja,Nx)+1,jb,1],[ja,jb,1],[mod(ja,Nx)+1,jb,1],[ja,jb,1]]) #The first two are orbitals, the last two are phonons
     push!(py_xbond,[[mod(ja,Nx)+1,jb,2],[ja,jb,2],[mod(ja,Nx)+1,jb,1],[ja,jb,1]])
  end
  
  for ja in 1:Nx,jb in 1:Ny
    push!(px_ybond,[[ja,mod(jb,Ny)+1,1],[ja,jb,1],[ja,mod(jb,Ny)+1,2],[ja,jb,2]])
    push!(py_ybond,[[ja,mod(jb,Ny)+1,2],[ja,jb,2],[ja,mod(jb,Ny)+1,2],[ja,jb,2]])
  end
  

  
  for ja in 1:Nx, jb in 1:Ny
    push!(NNN_sp_d1,[[mod(ja,Nx)+1,mod(jb,Ny)+1,1],[mod(ja,Nx)+1,mod(jb,Ny)+1,2],[ja,jb,1],[ja,jb,2]])
    push!(NNN_sp_d2,[[mod(ja,Nx)+1,mod(jb-2,Ny)+1,1],[ja,jb,2],[mod(ja,Nx)+1,mod(jb-2,Ny)+1,2],[ja,jb,1]])
  end
  for ja in 1:Nx,jb in 1:Ny
    push!(shear_ybond,[[ja,mod(jb,Ny)+1,1],[ja,jb,1],[ja,mod(jb,Ny)+1,1],[ja,jb,1]])
  end
  for ja in 1:Nx,jb in 1:Ny
    push!(shear_xbond,[[mod(ja,Nx)+1,jb,1],[ja,jb,1],[mod(ja,Nx)+1,jb,2],[ja,jb,2]]) #The first two are orbitals, the last two are phonons
 end

  
  
  
   return  px_xbond, px_ybond, py_xbond, py_ybond, NNN_sp_d1, NNN_sp_d2,kxspace,kyspace,shear_xbond,shear_ybond
end





function get_Lambdaset(Nx::Int,Ny::Int,α::Float64,β::Float64)::Vector{Any}
  px_xbond, px_ybond, py_xbond, py_ybond, NNN_sp_d1, NNN_sp_d2,kxspace,kyspace,shear_xbond,shear_ybond=initialize(Nx,Ny)
  Λset=[]
  for ja in eachindex(px_xbond), ui in 3:4
   if px_xbond[ja][ui][1:2]==[1,1] #Λset records the real space Lambda
      d=px_xbond[ja][ui][3]
      oi=px_xbond[ja][1][3]
      oj=px_xbond[ja][2][3]
      xi=px_xbond[ja][1][1]
      yi=px_xbond[ja][1][2]
      xj=px_xbond[ja][2][1]
      yj=px_xbond[ja][2][2]
      indexset=[d,oi,oj,xi,yi,xj,yj]
       push!(Λset,[(-1)^(ui)*α,indexset])
      indexset=[d,oj,oi,xj,yj,xi,yi]
      push!(Λset,[(-1)^(ui)*α,indexset])
   end
  
  end

  for ja in eachindex(px_ybond), ui in 3:4
    if px_ybond[ja][ui][1:2]==[1,1]
       d=px_ybond[ja][ui][3]
       oi=px_ybond[ja][1][3]
       oj=px_ybond[ja][2][3]
       xi=px_ybond[ja][1][1]
       yi=px_ybond[ja][1][2]
       xj=px_ybond[ja][2][1]
       yj=px_ybond[ja][2][2]
       indexset=[d,oi,oj,xi,yi,xj,yj]
        push!(Λset,[(-1)^(ui)*β,indexset])
       indexset=[d,oj,oi,xj,yj,xi,yi]
       push!(Λset,[(-1)^(ui)*β,indexset])
    end

  end

  for ja in eachindex(py_ybond), ui in 3:4
    if py_ybond[ja][ui][1:2]==[1,1]
       d=py_ybond[ja][ui][3]
       oi=py_ybond[ja][1][3]
       oj=py_ybond[ja][2][3]
       xi=py_ybond[ja][1][1]
       yi=py_ybond[ja][1][2]
       xj=py_ybond[ja][2][1]
       yj=py_ybond[ja][2][2]
       indexset=[d,oi,oj,xi,yi,xj,yj]
        push!(Λset,[(-1)^(ui)*α,indexset])
       indexset=[d,oj,oi,xj,yj,xi,yi]
       push!(Λset,[(-1)^(ui)*α,indexset])
    end
   
  end

  for ja in eachindex(py_xbond), ui in 3:4
    if py_xbond[ja][ui][1:2]==[1,1] #We are considering terms coupling to phonons living on (1,1 site)
       d=py_xbond[ja][ui][3]
       oi=py_xbond[ja][1][3]
       oj=py_xbond[ja][2][3]
       xi=py_xbond[ja][1][1]
       yi=py_xbond[ja][1][2]
       xj=py_xbond[ja][2][1]
       yj=py_xbond[ja][2][2]
       indexset=[d,oi,oj,xi,yi,xj,yj]
        push!(Λset,[(-1)^(ui)*β,indexset])
       indexset=[d,oj,oi,xj,yj,xi,yi]
       push!(Λset,[(-1)^(ui)*β,indexset])
    end   
 end

 return  Λset

end





function get_Lambdaset_reduandant(Nx::Int,Ny::Int,α::Float64,β::Float64)
  px_xbond, px_ybond, py_xbond, py_ybond, NNN_sp_d1, NNN_sp_d2,kxspace,kyspace,shear_xbond,shear_ybond=initialize(Nx,Ny)
  Λset=[]
  for ja in eachindex(px_xbond), ui in 3:4
   if px_xbond[ja][ui][1:2]==[1,1] #Λset records the real space Lambda
      d=1
      oi=px_xbond[ja][1][3]
      oj=px_xbond[ja][2][3]
      xi=px_xbond[ja][1][1]
      yi=px_xbond[ja][1][2]
      xj=px_xbond[ja][2][1]
      yj=px_xbond[ja][2][2]
      indexset=[d,oi,oj,xi,yi,xj,yj]
       push!(Λset,[(-1)^(ui),indexset])
      indexset=[d,oj,oi,xj,yj,xi,yi]
      push!(Λset,[(-1)^(ui),indexset])
   end
  
  end

  for ja in eachindex(px_ybond), ui in 3:4
    if px_ybond[ja][ui][1:2]==[1,1]
       d=4
       oi=px_ybond[ja][1][3]
       oj=px_ybond[ja][2][3]
       xi=px_ybond[ja][1][1]
       yi=px_ybond[ja][1][2]
       xj=px_ybond[ja][2][1]
       yj=px_ybond[ja][2][2]
       indexset=[d,oi,oj,xi,yi,xj,yj]
        push!(Λset,[(-1)^(ui),indexset])
       indexset=[d,oj,oi,xj,yj,xi,yi]
       push!(Λset,[(-1)^(ui),indexset])
    end

  end

  for ja in eachindex(py_ybond), ui in 3:4
    if py_ybond[ja][ui][1:2]==[1,1]
       d=3
       oi=py_ybond[ja][1][3]
       oj=py_ybond[ja][2][3]
       xi=py_ybond[ja][1][1]
       yi=py_ybond[ja][1][2]
       xj=py_ybond[ja][2][1]
       yj=py_ybond[ja][2][2]
       indexset=[d,oi,oj,xi,yi,xj,yj]
        push!(Λset,[(-1)^(ui),indexset])
       indexset=[d,oj,oi,xj,yj,xi,yi]
       push!(Λset,[(-1)^(ui),indexset])
    end
   
  end

  for ja in eachindex(py_xbond), ui in 3:4
    if py_xbond[ja][ui][1:2]==[1,1] #We are considering terms coupling to phonons living on (1,1 site)
       d=2
       oi=py_xbond[ja][1][3]
       oj=py_xbond[ja][2][3]
       xi=py_xbond[ja][1][1]
       yi=py_xbond[ja][1][2]
       xj=py_xbond[ja][2][1]
       yj=py_xbond[ja][2][2]
       indexset=[d,oi,oj,xi,yi,xj,yj]
        push!(Λset,[(-1)^(ui),indexset])
       indexset=[d,oj,oi,xj,yj,xi,yi]
       push!(Λset,[(-1)^(ui),indexset])
    end   
 end


 gmatrix=zeros(Float64,4,2)
 gmatrix[1,1]=α
 gmatrix[2,1]=β

 gmatrix[3,2]=α
 gmatrix[4,2]=β

 return  Λset,gmatrix

end















function  get_pertur_factor(ea::Float64,eb::Float64,FL::Float64,temp::Float64)::ComplexF64
   
  f1=1/(exp((ea-FL)/temp)+1)

  f2=1/(exp((eb-FL)/temp)+1)
 

 return   (f1-f2)/(ea-eb+im*10^(-8))
end



function Lambdaset_to_matrix(Nx::Int,Ny::Int,Λset::Vector{Any},tper::Float64,tpa::Float64,tNNN::Float64,temp::Float64)
  px_xbond, px_ybond, py_xbond, py_ybond, NNN_sp_d1, NNN_sp_d2,kxspace,kyspace= initialize(Nx,Ny)

  Λmatrix=zeros(ComplexF64,2,2,2,Nx,Ny,Nx,Ny) # the index is d,o1,o2,k1,k2
  for k1x in 1:Nx, k1y in 1:Ny, k2x in 1:Nx, k2y in 1:Ny
     k1=[kxspace[k1x],kyspace[k1y]]
     k2=[kxspace[k2x],kyspace[k2y]]
     for ja in eachindex(Λset)
       R1=[Λset[ja][2][4],Λset[ja][2][5]]-[1,1]
       R2=[Λset[ja][2][6],Λset[ja][2][7]]-[1,1]
       d=Λset[ja][2][1]
       o1=Λset[ja][2][2]
       o2=Λset[ja][2][3]
       Λmatrix[d,o1,o2,k1x,k1y,k2x,k2y]+=1/sqrt(Nx*Ny)*exp(-im*dot(k1,R1))*exp(im*dot(k2,R2))*Λset[ja][1]
     end
  end

  Λmatrix_bandbasis=zeros(ComplexF64,2,2,2,Nx,Ny,Nx,Ny) #the index order is d, n1,n2,k1x,k1y,k2x,k2y,q=k_1-k_2
 for k1x in 1:Nx, k1y in 1:Ny, k2x in 1:Nx, k2y in 1:Ny
  k1=[kxspace[k1x],kyspace[k1y]]
  k2=[kxspace[k2x],kyspace[k2y]]
  FFF=eigen(get_Hamiltonian(k1,tper,tpa,tNNN))
  vectork1=FFF.vectors
  FFF=eigen(get_Hamiltonian(k2,tper,tpa,tNNN))
  vectork2=FFF.vectors
  for di in 1:2
    Λmatrix_bandbasis[di,:,:,k1x,k1y,k2x,k2y]+=vectork1'*Λmatrix[di,:,:,k1x,k1y,k2x,k2y]*vectork2
  end
 end


 electron_spectrum=zeros(Float64,2,Nx,Ny)
 for k1x in 1:Nx, k1y in 1:Ny
    k1=[kxspace[k1x],kyspace[k1y]]
    FFF=eigen(get_Hamiltonian(k1,tper,tpa,tNNN))
    electron_spectrum[:,k1x,k1y]=real.(FFF.values)
 end

 FL=sort(vec(electron_spectrum))[Int(round(Nx*Ny*filling))]
  return Λmatrix_bandbasis,electron_spectrum,FL

end



function Lambdaset_to_matrix_redundant(Nx::Int,Ny::Int,Λset::Vector{Any},tper::Float64,tpa::Float64,tNNN::Float64,temp::Float64)
  px_xbond, px_ybond, py_xbond, py_ybond, NNN_sp_d1, NNN_sp_d2,kxspace,kyspace= initialize(Nx,Ny)

  Λmatrix=zeros(ComplexF64,4,2,2,Nx,Ny,Nx,Ny) # the index is d,o1,o2,k1,k2
  for k1x in 1:Nx, k1y in 1:Ny, k2x in 1:Nx, k2y in 1:Ny
     k1=[kxspace[k1x],kyspace[k1y]]
     k2=[kxspace[k2x],kyspace[k2y]]
     for ja in eachindex(Λset)
       R1=[Λset[ja][2][4],Λset[ja][2][5]]-[1,1]
       R2=[Λset[ja][2][6],Λset[ja][2][7]]-[1,1]
       d=Λset[ja][2][1]
       o1=Λset[ja][2][2]
       o2=Λset[ja][2][3]
       Λmatrix[d,o1,o2,k1x,k1y,k2x,k2y]+=1/sqrt(Nx*Ny)*exp(-im*dot(k1,R1))*exp(im*dot(k2,R2))*Λset[ja][1]
     end
  end

  Λmatrix_bandbasis=zeros(ComplexF64,4,2,2,Nx,Ny,Nx,Ny) #the index order is d, n1,n2,k1x,k1y,k2x,k2y,q=k_1-k_2
 for k1x in 1:Nx, k1y in 1:Ny, k2x in 1:Nx, k2y in 1:Ny
  k1=[kxspace[k1x],kyspace[k1y]]
  k2=[kxspace[k2x],kyspace[k2y]]
  FFF=eigen(get_Hamiltonian(k1,tper,tpa,tNNN))
  vectork1=FFF.vectors
  FFF=eigen(get_Hamiltonian(k2,tper,tpa,tNNN))
  vectork2=FFF.vectors
  for di in 1:4
    Λmatrix_bandbasis[di,:,:,k1x,k1y,k2x,k2y]+=vectork1'*Λmatrix[di,:,:,k1x,k1y,k2x,k2y]*vectork2
  end
 end


 electron_spectrum=zeros(Float64,2,Nx,Ny)
 for k1x in 1:Nx, k1y in 1:Ny
    k1=[kxspace[k1x],kyspace[k1y]]
    FFF=eigen(get_Hamiltonian(k1,tper,tpa,tNNN))
    electron_spectrum[:,k1x,k1y]=real.(FFF.values)
 end

 FL=sort(vec(electron_spectrum))[Int(round(Nx*Ny*filling))]
  return Λmatrix_bandbasis,electron_spectrum,FL

end






function Lambda_to_Keff(Λmatrix_bandbasis::Array{ComplexF64},electron_spectrum::Array{Float64},FL::Float64,temp::Float64)::Array{ComplexF64}
  Keff_momentum=zeros(ComplexF64,2,2,Nx,Ny)
 for  qx in 1:Nx, qy in 1:Ny
    for k1x in 1:Nx, k1y in 1:Ny, n1 in 1:2, n2 in 1:2
        k2x=mod((k1x-1)+(qx-1),Nx)+1 #k2 is kalpha, k1 is kgamma
        k2y=mod((k1y-1)+(qy-1),Ny)+1
        e2=electron_spectrum[n2,k2x,k2y]
        e1=electron_spectrum[n1,k1x,k1y]
        factor=get_pertur_factor(e2,e1,FL,temp)
        
        Keff_momentum[:,:,qx,qy]+=Λmatrix_bandbasis[:,n1,n2,k1x,k1y,k2x,k2y]*transpose(Λmatrix_bandbasis[:,n2,n1,k2x,k2y,k1x,k1y])*factor
   
    end
 end
   return Keff_momentum
end



function Lambda_to_Keff_redundant(Λmatrix_bandbasis::Array{ComplexF64},electron_spectrum::Array{Float64},FL::Float64,temp::Float64)::Array{ComplexF64}
  Keff_momentum=zeros(ComplexF64,4,4,Nx,Ny)
 for  qx in 1:Nx, qy in 1:Ny
    for k1x in 1:Nx, k1y in 1:Ny, n1 in 1:2, n2 in 1:2
        k2x=mod((k1x-1)+(qx-1),Nx)+1 #k2 is kalpha, k1 is kgamma
        k2y=mod((k1y-1)+(qy-1),Ny)+1
        e2=electron_spectrum[n2,k2x,k2y]
        e1=electron_spectrum[n1,k1x,k1y]
        factor=get_pertur_factor(e2,e1,FL,temp)
        
        Keff_momentum[:,:,qx,qy]+=Λmatrix_bandbasis[:,n1,n2,k1x,k1y,k2x,k2y]*transpose(Λmatrix_bandbasis[:,n2,n1,k2x,k2y,k1x,k1y])*factor
   
    end
 end
   return Keff_momentum
end




function barephonon(Nx::Int64,Ny::Int64,K::Float64,KNNN::Float64,shearstrength::Float64)::Array{ComplexF64}

  px_xbond, px_ybond, py_xbond, py_ybond, NNN_sp_d1, NNN_sp_d2,kxspace,kyspace,shear_xbond,shear_ybond=initialize(Nx,Ny)

  Kbare_real=zeros(Float64,Nx,Ny,2,Nx,Ny,2)
  for ja in eachindex(px_xbond)
     bond=px_xbond[ja]
     for jc in 3:4, jd in 3:4
     Kbare_real[bond[jc][1],bond[jc][2],bond[jc][3],bond[jd][1],bond[jd][2],bond[jd][3]]+=K*(-1)^(jc+jd)
     end
  end
  
  
  for ja in eachindex(px_ybond)
     bond=px_ybond[ja]
     for jc in 3:4, jd in 3:4
     Kbare_real[bond[jc][1],bond[jc][2],bond[jc][3],bond[jd][1],bond[jd][2],bond[jd][3]]+=K*(-1)^(jc+jd)
     end
  end
  
  for ja in eachindex(shear_xbond)
    bond=shear_xbond[ja]
    for jc in 3:4, jd in 3:4
    Kbare_real[bond[jc][1],bond[jc][2],bond[jc][3],bond[jd][1],bond[jd][2],bond[jd][3]]+=shearstrength*(-1)^(jc+jd)
    end
 end
 
 for ja in eachindex(shear_ybond)
    bond=shear_ybond[ja]
    for jc in 3:4, jd in 3:4
    Kbare_real[bond[jc][1],bond[jc][2],bond[jc][3],bond[jd][1],bond[jd][2],bond[jd][3]]+=shearstrength*(-1)^(jc+jd)
    end
 end

  
  for ja in eachindex(NNN_sp_d1)
     bond=NNN_sp_d1[ja]
     sig=[1,1,-1,-1]
     for jc in 1:4, jd in 1:4
     Kbare_real[bond[jc][1],bond[jc][2],bond[jc][3],bond[jd][1],bond[jd][2],bond[jd][3]]+=KNNN*sig[jc]*sig[jd]
     end
  end
  
  for ja in eachindex(NNN_sp_d2)
     bond=NNN_sp_d2[ja]
     sig=[1,1,-1,-1]
     for jc in 1:4, jd in 1:4
     Kbare_real[bond[jc][1],bond[jc][2],bond[jc][3],bond[jd][1],bond[jd][2],bond[jd][3]]+=KNNN*sig[jc]*sig[jd]
     end
  end
  
  Kbare_momentum=zeros(ComplexF64,2,2,Nx,Ny)
  
  for ja in 1:Nx, jb in 1:Ny, jc in 1:Nx, jd in 1:Ny
      qvec=[kxspace[ja],kyspace[jb]]
      Rvec=[jc,jd]
      RR1vecpos=[1+mod(jc,Nx),1+mod(jd,Ny)]
      Kbare_momentum[:,:,ja,jb]+=exp(im*dot(qvec,Rvec))*Kbare_real[1,1,:,RR1vecpos[1],RR1vecpos[2],:]
  end
  Kbare_real=nothing #realease the RAM here

  return Kbare_momentum

end

function get_spectrum(Kbare_momentum::Array{ComplexF64},Keff_momentum::Array{ComplexF64},Nx::Int,Ny::Int)::Tuple{Array{ComplexF64},Array{ComplexF64},Array{ComplexF64}}
  spectrum=zeros(ComplexF64,2,Nx,Ny)
 for ja in 1:Nx, jb in 1:Ny
   FFF=eigen(Kbare_momentum[:,:,ja,jb]+Keff_momentum[:,:,ja,jb])
   spectrum[:,ja,jb]=FFF.values
 end



 unperturbed_spectrum=zeros(ComplexF64,2,Nx,Ny)
 for ja in 1:Nx, jb in 1:Ny
   FFF=eigen(Kbare_momentum[:,:,ja,jb])
   unperturbed_spectrum[:,ja,jb]=FFF.values
 end

 χspectrum=zeros(ComplexF64,2,Nx,Ny)

 for ja in 1:Nx, jb in 1:Ny
  FFF=eigen(-Keff_momentum[:,:,ja,jb])
  χspectrum[:,ja,jb]=FFF.values
 end

 



 return spectrum,unperturbed_spectrum,χspectrum
end


function get_spectrum_redundant(Kbare_momentum::Array{ComplexF64},Keff_momentum::Array{ComplexF64},Nx::Int,Ny::Int,gmatrix::Matrix{Float64})::Tuple{Array{ComplexF64},Array{ComplexF64},Array{ComplexF64}}
  spectrum=zeros(ComplexF64,2,Nx,Ny)
 for ja in 1:Nx, jb in 1:Ny
   FFF=eigen(Kbare_momentum[:,:,ja,jb]+gmatrix'*Keff_momentum[:,:,ja,jb]*gmatrix)
   spectrum[:,ja,jb]=FFF.values
 end



 unperturbed_spectrum=zeros(ComplexF64,2,Nx,Ny)
 for ja in 1:Nx, jb in 1:Ny
   FFF=eigen(Kbare_momentum[:,:,ja,jb])
   unperturbed_spectrum[:,ja,jb]=FFF.values
 end

 χspectrum=zeros(ComplexF64,2,Nx,Ny)

 for ja in 1:Nx, jb in 1:Ny
  FFF=eigen(-gmatrix'*Keff_momentum[:,:,ja,jb]*gmatrix)
  χspectrum[:,ja,jb]=FFF.values
 end

 



 return spectrum,unperturbed_spectrum,χspectrum
end