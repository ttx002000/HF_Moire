using LinearAlgebra


function construct_Ham(px_xbond::Vector{Vector{Int}},px_ybond::Vector{Vector{Int}},phonon_coor::Vector{Float64},Nx::Int,Ny::Int,α::Float64,β::Float64)::Matrix{ComplexF64}
    Hphonon=zeros(ComplexF64,Nx*Ny,Nx*Ny)


    
    for ja in eachindex(px_xbond) #alpha=\beta for this case
    Hphonon[px_xbond[ja][1],px_xbond[ja][2]]-=α*(phonon_coor[px_xbond[ja][3]]-phonon_coor[px_xbond[ja][4]])
   
    end

    for ja in eachindex(px_ybond)
      Hphonon[px_ybond[ja][1],px_ybond[ja][2]]-=β*(phonon_coor[px_ybond[ja][3]]-phonon_coor[px_ybond[ja][4]])
     
    end

    Hphonon+=Hphonon';
  return Hphonon
end


function findFL(Nelec::Int,spectrum::Vector{Float64},temp::Float64,val_s::Float64,val_e::Float64)
  fl=0

  try_FL=(val_s+val_e)/2
  for ja in eachindex(spectrum)
    fd=1/(1+exp((spectrum[ja]- try_FL)/temp))
    fl+=real(fd)
  end


   if abs(fl-Nelec)<10^(-8)
    return try_FL
  elseif fl-Nelec>=10^(-8)
    return findFL(Nelec,spectrum,temp,val_s, try_FL)
  elseif fl-Nelec<=-10^(-8)
    return findFL(Nelec,spectrum,temp, try_FL,val_e)
   end

end


function calculate_gradient(K::Float64,KNNN::Float64,NNN_sp_d1::Vector{Vector{Int}},NNN_sp_d2::Vector{Vector{Int}},px_xbond::Vector{Vector{Int}},px_ybond::Vector{Vector{Int}},Htotal::Matrix{ComplexF64},Nx::Int,Ny::Int,orbital_id::Array{Int},phonon_id::Array{Int},phonon_coor::Vector{Float64},Nelec::Int,α::Float64,β::Float64)::Tuple{Vector{Float64},Float64,Float64,Float64}

    

    FFF=eigen(Htotal)
    spectrum=real.(FFF.values)
    #FL=spectrum[Nelec]

    temp=10^(-7)
    FL=findFL(Nelec,spectrum,temp,spectrum[Nelec]-0.1,spectrum[Nelec]+0.1)

    


    Egap=real(FFF.values[Nelec+1]-FFF.values[Nelec])
    println("gap=",Egap)

    E0=0
   
    for ja in eachindex(spectrum)
     fd=1/(1+exp((spectrum[ja]-FL)/temp))
     E0+=real(spectrum[ja]*fd)
     
    end

    E_elec=copy(E0)
    for ja in eachindex(px_xbond)
      E0+=K/2*(phonon_coor[px_xbond[ja][3]]-phonon_coor[px_xbond[ja][4]])^2
    end
    for ja in eachindex(px_ybond)
      E0+=K/2*(phonon_coor[px_ybond[ja][3]]-phonon_coor[px_ybond[ja][4]])^2
    end
  
    for ja in eachindex(NNN_sp_d1)
      E0+=KNNN/2*(phonon_coor[NNN_sp_d1[ja][1]]+phonon_coor[NNN_sp_d1[ja][2]]-phonon_coor[NNN_sp_d1[ja][3]]-phonon_coor[NNN_sp_d1[ja][4]])^2
      E0+=KNNN/2*(phonon_coor[NNN_sp_d2[ja][1]]+phonon_coor[NNN_sp_d2[ja][2]]-phonon_coor[NNN_sp_d2[ja][3]]-phonon_coor[NNN_sp_d2[ja][4]])^2
    end
    
    
    #This is the phonon part
    gradient=zeros(Float64,2*Nx*Ny)
    for ja in eachindex(px_xbond)
      gradient[px_xbond[ja][3]]+=K*(phonon_coor[px_xbond[ja][3]]-phonon_coor[px_xbond[ja][4]])
      gradient[px_xbond[ja][4]]+=K*(phonon_coor[px_xbond[ja][4]]-phonon_coor[px_xbond[ja][3]])
    end
  
    for ja in eachindex(px_ybond)
      gradient[px_ybond[ja][3]]+=K*(phonon_coor[px_ybond[ja][3]]-phonon_coor[px_ybond[ja][4]])
      gradient[px_ybond[ja][4]]+=K*(phonon_coor[px_ybond[ja][4]]-phonon_coor[px_ybond[ja][3]])
    end
   
    for ja in eachindex(NNN_sp_d1)
      ss=KNNN*(phonon_coor[NNN_sp_d1[ja][1]]+phonon_coor[NNN_sp_d1[ja][2]]-phonon_coor[NNN_sp_d1[ja][3]]-phonon_coor[NNN_sp_d1[ja][4]])
      gradient[NNN_sp_d1[ja][1]]+=ss
      gradient[NNN_sp_d1[ja][2]]+=ss
      gradient[NNN_sp_d1[ja][3]]-=ss
      gradient[NNN_sp_d1[ja][4]]-=ss
    end
  
    for ja in eachindex(NNN_sp_d2)
      ss=KNNN*(phonon_coor[NNN_sp_d2[ja][1]]+phonon_coor[NNN_sp_d2[ja][2]]-phonon_coor[NNN_sp_d2[ja][3]]-phonon_coor[NNN_sp_d2[ja][4]])
      gradient[NNN_sp_d2[ja][1]]+=ss
      gradient[NNN_sp_d2[ja][2]]+=ss
      gradient[NNN_sp_d2[ja][3]]-=ss
      gradient[NNN_sp_d2[ja][4]]-=ss
    end
   
  
    #This is the part for the ux 
   for jelec in eachindex(spectrum)
    fd=1/(1+exp((spectrum[jelec]-FL)/temp))
    for ja in eachindex(px_xbond)
     gradient[px_xbond[ja][3]]-=2*α*real(conj(FFF.vectors[px_xbond[ja][1],jelec])*FFF.vectors[px_xbond[ja][2],jelec])*fd
     gradient[px_xbond[ja][4]]+=2*α*real(conj(FFF.vectors[px_xbond[ja][1],jelec])*FFF.vectors[px_xbond[ja][2],jelec])*fd
    end
   end
  
   for jelec in eachindex(spectrum)
    fd=1/(1+exp((spectrum[jelec]-FL)/temp))
    for ja in eachindex(px_ybond)
     gradient[px_ybond[ja][3]]-=2*β*real(conj(FFF.vectors[px_ybond[ja][1],jelec])*FFF.vectors[px_ybond[ja][2],jelec])*fd
     gradient[px_ybond[ja][4]]+=2*β*real(conj(FFF.vectors[px_ybond[ja][1],jelec])*FFF.vectors[px_ybond[ja][2],jelec])*fd
    end
  end
  
  
  
  
    return gradient,E0/(Nx*Ny),E_elec/(Nx*Ny),Egap
  
end


function initialize(Nx::Int64,Ny::Int64,tpa::Float64)
   
   
    
    orbital_id=zeros(Int,Nx,Ny)
    for ja in 1:Nx, jb in 1:Ny
        orbital_id[ja,jb]=ja+(jb-1)*Nx
    end
    
    
    phonon_id=zeros(Int,Nx,Ny,2)
    
    for ja in 1:Nx, jb in 1:Ny, jo in 1:2
        phonon_id[ja,jb,jo]=ja+(jb-1)*Nx+(jo-1)*Nx*Ny
    end
    
    
    
  
    #The lattice constant is around 4A

    
    px_xbond=Vector{Int}[]
    px_ybond=Vector{Int}[]
    
       
    NNN_sp_d1=Vector{Int}[]
    NNN_sp_d2=Vector{Int}[]
    
    
    for ja in 1:Nx,jb in 1:Ny
       push!(px_xbond,[orbital_id[mod(ja,Nx)+1,jb],orbital_id[ja,jb],phonon_id[mod(ja,Nx)+1,jb,1],phonon_id[ja,jb,1]])
    end
    
    for ja in 1:Nx,jb in 1:Ny
      push!(px_ybond,[orbital_id[ja,mod(jb,Ny)+1],orbital_id[ja,jb],phonon_id[ja,mod(jb,Ny)+1,2],phonon_id[ja,jb,2]])
    end
    
 
    
    for ja in 1:Nx, jb in 1:Ny
      push!(NNN_sp_d1,[phonon_id[mod(ja,Nx)+1,mod(jb,Ny)+1,1],phonon_id[mod(ja,Nx)+1,mod(jb,Ny)+1,2],phonon_id[ja,jb,1],phonon_id[ja,jb,2]])
      push!(NNN_sp_d2,[phonon_id[mod(ja,Nx)+1,mod(jb-2,Ny)+1,1],phonon_id[ja,jb,2],phonon_id[mod(ja,Nx)+1,mod(jb-2,Ny)+1,2],phonon_id[ja,jb,1]])
    end
    
    
    
    H0=zeros(ComplexF64,Nx*Ny,Nx*Ny)
    
    for ja in eachindex(px_xbond)
       H0[px_xbond[ja][1],px_xbond[ja][2]]-=tpa;
       H0[px_ybond[ja][1],px_ybond[ja][2]]-=tpa;
    
    end
    
    
  
    H0+=H0';
    
    
     return H0, orbital_id, phonon_id, px_xbond, px_ybond,NNN_sp_d1, NNN_sp_d2
end


function resh_phonon(phonon_coor,phonon_id,Nx,Ny)


  atom_x=zeros(Float64,Nx,Ny)
  atom_y=zeros(Float64,Nx,Ny)
  dis_x=zeros(Float64,Nx,Ny)
  dis_y=zeros(Float64,Nx,Ny)
  
  for ja in 1:Nx, jb in 1:Ny
      atom_x[ja,jb]=ja
      atom_y[ja,jb]=jb
      dis_x[ja,jb]=phonon_coor[phonon_id[ja,jb,1]]
      dis_y[ja,jb]=phonon_coor[phonon_id[ja,jb,2]]
      
  end
  COM_x=sum(dis_x)/(Nx*Ny)
  COM_y=sum(dis_y)/(Nx*Ny)
  dis_x=dis_x .- COM_x
  dis_y=dis_y .- COM_y
  max_record=sort(vec(sqrt.(dis_x.^2+dis_y.^2)))[Nx*Ny]
  average_record=sum(vec(sqrt.(dis_x.^2+dis_y.^2)))/(Nx*Ny)
  return dis_x,dis_y,max_record,average_record
end



function iteration(Nx::Int,Ny::Int,Nelec::Int,px_xbond::Vector{Vector{Int}},px_ybond::Vector{Vector{Int}},NNN_sp_d1::Vector{Vector{Int}},NNN_sp_d2::Vector{Vector{Int}},orbital_id::Array{Int},phonon_id::Array{Int},α::Float64,β::Float64,K::Float64,KNNN::Float64,H0::Matrix{ComplexF64})
    #phonon_coor=randn(2*Nx*Ny)*10^(-1)
    phonon_coor=zeros(2*Nx*Ny)
    E_old=10^8
    E_new=0.0
    grad_old=zeros(Float64,2*Nx*Ny).+10
    grad_new=zeros(Float64,2*Nx*Ny).+10

    update_rate=1.0
    Eelec_new=0.0
    
    Hph=0.0
    Egap=0.0
    
    
    
    
    
    itcount=0
     while norm(grad_new)>1*10^(-6) && update_rate>10^(-5)
      itcount+=1
       println("iterations",itcount)  
        Hph=construct_Ham(px_xbond,px_ybond,phonon_coor,Nx,Ny,α,β) #I modified the order between py_xbond and py_ybond
          grad_new,E_new,Eelec_new,Egap=calculate_gradient(K,KNNN,NNN_sp_d1,NNN_sp_d2,px_xbond,px_ybond,H0+Hph,Nx,Ny,orbital_id,phonon_id,phonon_coor,Nelec,α,β)
          
          println("norm=",norm(grad_new))
           
          println("Etotal=",E_new,"Eelec=",Eelec_new)
      
          if (E_new>E_old)
            update_rate=0.8*update_rate
            println("update rate adjusted to be","$(update_rate)")
          elseif ((E_new<E_old)&&(norm(grad_old)>norm(grad_new)))||(abs(E_new-E_old)<10^(-10)) 
            update_rate=update_rate*1.01
            println("update rate adjusted to be","$(update_rate)")
      
          end
      
        
      
          E_old=E_new
          grad_old=grad_new
      
         
          phonon_coor=phonon_coor-grad_new*update_rate
          COM_x=sum(phonon_coor[vec(phonon_id[:,:,1])])/(Nx*Ny)
          COM_y=sum(phonon_coor[vec(phonon_id[:,:,2])])/(Nx*Ny)
          phonon_coor[vec(phonon_id[:,:,1])]=phonon_coor[vec(phonon_id[:,:,1])] .- COM_x
          phonon_coor[vec(phonon_id[:,:,2])]=phonon_coor[vec(phonon_id[:,:,2])] .- COM_y
      end


      temp=10^(-7)
      FFF=eigen(H0+Hph)
      spectrum=FFF.values
      #FL=spectrum[Nelec]
      FL=findFL(Nelec,spectrum,temp,spectrum[Nelec]-0.1,spectrum[Nelec]+0.1)
   
   
  
       ave_npa=0
     
      for ja in eachindex(spectrum)
       fd=1/(1+exp((spectrum[ja]-FL)/temp))
       ave_npa+=fd

      end













       
      return phonon_coor, Hph, grad_old, E_old, Eelec_new, Egap,ave_npa
end


