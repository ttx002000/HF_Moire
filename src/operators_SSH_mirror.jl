using LinearAlgebra


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

function get_mirror_oneminusone(phonon_id,phonon_coor,Nx,Ny,center_x,center_y)
    if Nx≠Ny
        println("error")
       return 10^8
    end
    # Let's first do 1,1 get_mirror
    phonon_id_mirror=zeros(Int,Nx,Ny,2)
    for ja in 1:Nx, jb in 1:Ny
      mirror_x=mod(center_x+center_y-jb-1,Nx)+1
      mirror_y=mod(center_x+center_y-ja-1,Ny)+1

      phonon_id_mirror[ja,jb,1]=phonon_id[mirror_x,mirror_y,2]
      phonon_id_mirror[ja,jb,2]=phonon_id[mirror_x,mirror_y,1]


    end
    phonon_coor_mirror=zeros(Float64,Nx*Ny*2)
    for ja in 1:Nx, jb in 1:Ny, jc in 1:2
       phonon_coor_mirror[phonon_id[ja,jb,jc]]=-phonon_coor[phonon_id_mirror[ja,jb,jc]]
    end

  return phonon_coor_mirror,sum(abs.(phonon_coor_mirror-phonon_coor))

end



function get_mirror_oneone(phonon_id,phonon_coor,Nx,Ny,center_x,center_y)
    if Nx≠Ny
        println("error")
       return 10^8
    end
    # Let's first do 1,1 get_mirror
    phonon_id_mirror=zeros(Int,Nx,Ny,2)
    for ja in 1:Nx, jb in 1:Ny
      mirror_x=mod(center_x-center_y+jb-1,Nx)+1
      mirror_y=mod(center_y-center_x+ja-1,Ny)+1

      phonon_id_mirror[ja,jb,1]=phonon_id[mirror_x,mirror_y,2]
      phonon_id_mirror[ja,jb,2]=phonon_id[mirror_x,mirror_y,1]
    end
    phonon_coor_mirror=zeros(Float64,Nx*Ny*2)
    for ja in 1:Nx, jb in 1:Ny, jc in 1:2
       phonon_coor_mirror[phonon_id[ja,jb,jc]]=phonon_coor[phonon_id_mirror[ja,jb,jc]]
    end

  return phonon_coor_mirror,sum(abs.(phonon_coor_mirror-phonon_coor))

end






function calculate_freeenergy(K::Float64,KNNN::Float64,NNN_sp_d1::Vector{Vector{Int}},
  NNN_sp_d2::Vector{Vector{Int}},px_xbond::Vector{Vector{Int}},
  px_ybond::Vector{Vector{Int}},py_xbond::Vector{Vector{Int}},
  py_ybond::Vector{Vector{Int}},Htotal::Matrix{ComplexF64},
  Nx::Int,Ny::Int,orbital_id::Array{Int},phonon_id::Array{Int},
  phonon_coor::Vector{Float64},Nelec::Float64,α::Float64,
  β::Float64,gshear::Float64,temp::Float64)



 FFF=eigen(Htotal)
 spectrum=real.(FFF.values)

 Egap=real(FFF.values[Int(round(Nelec))+1]-FFF.values[Int(round(Nelec))])
 FL=findFL(Nelec,spectrum,temp,spectrum[1],spectrum[length(spectrum)])



 FL_list=zeros(Float64,2*Nx*Ny)

 for ja in eachindex(spectrum)

 FL_list[ja]=1/(1+exp((spectrum[ja]-FL)/temp))

 end
 E0=real(sum(spectrum .* FL_list))
 Eelec=copy(E0)


 for ja in eachindex(px_xbond)
 E0+=K/2*(phonon_coor[px_xbond[ja][3]]-phonon_coor[px_xbond[ja][4]])^2+gshear*(phonon_coor[px_xbond[ja][5]]-phonon_coor[px_xbond[ja][6]])^2
 end
 for ja in eachindex(px_ybond)
 E0+=K/2*(phonon_coor[px_ybond[ja][3]]-phonon_coor[px_ybond[ja][4]])^2+gshear*(phonon_coor[px_ybond[ja][5]]-phonon_coor[px_ybond[ja][6]])^2
 end

 for ja in eachindex(NNN_sp_d1)
 E0+=KNNN/2*(phonon_coor[NNN_sp_d1[ja][1]]+phonon_coor[NNN_sp_d1[ja][2]]-phonon_coor[NNN_sp_d1[ja][3]]-phonon_coor[NNN_sp_d1[ja][4]])^2
 E0+=KNNN/2*(phonon_coor[NNN_sp_d2[ja][1]]+phonon_coor[NNN_sp_d2[ja][2]]-phonon_coor[NNN_sp_d2[ja][3]]-phonon_coor[NNN_sp_d2[ja][4]])^2
 end

 grand_po=0.0
 for ja in eachindex(spectrum)
 if spectrum[ja]<FL
 grand_po+=(spectrum[ja]-FL)-temp*log(1+exp((spectrum[ja]-FL)/temp))
 else
 grand_po+=-temp*log(1+exp(-(spectrum[ja]-FL)/temp))
 end   
 end

 free_energy=(E0-Eelec)+grand_po+FL*Nelec













  return free_energy/(Nx*Ny)

end


function initialize(Nx::Int64,Ny::Int64,tper::Float64,tpa::Float64,tNNN::Float64)
   
   
    
    orbital_id=zeros(Int,Nx,Ny,2)
    for ja in 1:Nx, jb in 1:Ny, jo in 1:2
        orbital_id[ja,jb,jo]=ja+(jb-1)*Nx+(jo-1)*Nx*Ny
    end
    
    
    phonon_id=zeros(Int,Nx,Ny,2)
    
    for ja in 1:Nx, jb in 1:Ny, jo in 1:2
        phonon_id[ja,jb,jo]=ja+(jb-1)*Nx+(jo-1)*Nx*Ny
    end
    
    
    
  
    #The lattice constant is around 4A

    
    px_xbond=Vector{Int}[]
    px_ybond=Vector{Int}[]
    py_xbond=Vector{Int}[]
    py_ybond=Vector{Int}[]
       
    NNN_sp_d1=Vector{Int}[]
    NNN_sp_d2=Vector{Int}[]
    
    
    for ja in 1:Nx,jb in 1:Ny
       push!(px_xbond,[orbital_id[mod(ja,Nx)+1,jb,1],orbital_id[ja,jb,1],phonon_id[mod(ja,Nx)+1,jb,1],phonon_id[ja,jb,1],phonon_id[mod(ja,Nx)+1,jb,2],phonon_id[ja,jb,2]])
       push!(py_xbond,[orbital_id[mod(ja,Nx)+1,jb,2],orbital_id[ja,jb,2],phonon_id[mod(ja,Nx)+1,jb,1],phonon_id[ja,jb,1],phonon_id[mod(ja,Nx)+1,jb,2],phonon_id[ja,jb,2]])
    end
    
    for ja in 1:Nx,jb in 1:Ny
      push!(px_ybond,[orbital_id[ja,mod(jb,Ny)+1,1],orbital_id[ja,jb,1],phonon_id[ja,mod(jb,Ny)+1,2],phonon_id[ja,jb,2],phonon_id[ja,mod(jb,Ny)+1,1],phonon_id[ja,jb,1]])
      push!(py_ybond,[orbital_id[ja,mod(jb,Ny)+1,2],orbital_id[ja,jb,2],phonon_id[ja,mod(jb,Ny)+1,2],phonon_id[ja,jb,2],phonon_id[ja,mod(jb,Ny)+1,1],phonon_id[ja,jb,1]])
    end
    
 
    
    for ja in 1:Nx, jb in 1:Ny
      push!(NNN_sp_d1,[phonon_id[mod(ja,Nx)+1,mod(jb,Ny)+1,1],phonon_id[mod(ja,Nx)+1,mod(jb,Ny)+1,2],phonon_id[ja,jb,1],phonon_id[ja,jb,2]])
      push!(NNN_sp_d2,[phonon_id[mod(ja,Nx)+1,mod(jb-2,Ny)+1,1],phonon_id[ja,jb,2],phonon_id[mod(ja,Nx)+1,mod(jb-2,Ny)+1,2],phonon_id[ja,jb,1]])
    end
    
    
    
    H0=zeros(ComplexF64,2*Nx*Ny,2*Nx*Ny)
    
    for ja in eachindex(px_xbond)
       H0[px_xbond[ja][1],px_xbond[ja][2]]-=tpa;
      
       H0[py_xbond[ja][1],py_xbond[ja][2]]-=tper;
     
    
    end

    for ja in eachindex(px_ybond)
  
      H0[py_ybond[ja][1],py_ybond[ja][2]]-=tpa;
    
      H0[px_ybond[ja][1],px_ybond[ja][2]]-=tper;
   
   end
    
    
    
    for ja in 1:Nx,jb in 1:Ny
    
      H0[orbital_id[mod(ja,Nx)+1,mod(jb,Ny)+1,2],orbital_id[ja,jb,1]]-=tNNN
      H0[orbital_id[mod(ja,Nx)+1,mod(jb-2,Ny)+1,2],orbital_id[ja,jb,1]]-=(-tNNN) #for mirrow symmetry
      H0[orbital_id[mod(ja-2,Nx)+1,mod(jb,Ny)+1,2],orbital_id[ja,jb,1]]-=(-tNNN)
      H0[orbital_id[mod(ja-2,Nx)+1,mod(jb-2,Ny)+1,2],orbital_id[ja,jb,1]]-=tNNN
    
    end
    
    
    
    H0+=H0';

    return H0, orbital_id, phonon_id, px_xbond, px_ybond, py_xbond, py_ybond, NNN_sp_d1, NNN_sp_d2
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



function construct_sym(unsym_phonon_coor,Nx,Ny,tper,tpa,tNNN,α,β,K,KNNN,gshear,temp,Nelec)

  diff_try=zeros(Float64,Nx,Ny)
  free_energy_symmetric=zeros(Float64,Nx,Ny)
  symmetrized_phonon_coor=Array{Vector{Float64}}(undef,Nx,Ny)



  H0, orbital_id, phonon_id, px_xbond, px_ybond, py_xbond, py_ybond, NNN_sp_d1, NNN_sp_d2=initialize(Nx,Ny,tper,tpa,tNNN)

    
    for jb in 1:Nx, jc in 1:Ny
      println("jb=",jb," jc=",jc)
      flush(stdout)
        center_x=jb
        center_y=jc
      
     tryconfig1,_=get_mirror_oneone(phonon_id,unsym_phonon_coor,Nx,Ny,center_x,center_y)
     phonon_coord_try1=(unsym_phonon_coor+tryconfig1)/2

      tryconfig2,_=get_mirror_oneminusone(phonon_id,phonon_coord_try1,Nx,Ny,center_x,center_y)
     phonon_coord_try2=( phonon_coord_try1+tryconfig2)/2
  
      
      diff_try[jb,jc]=sum(abs.(phonon_coord_try2-unsym_phonon_coor).^2)
     symmetrized_phonon_coor[jb,jc]=phonon_coord_try2
      Hph=construct_Ham(px_xbond,px_ybond,py_xbond,py_ybond, phonon_coord_try2,Nx,Ny,α,β) #I modified the order between py_xbond and py_ybond
      Htotal=Hph+H0
  



    

      free_energy_symmetric[jb,jc]=calculate_freeenergy(K,KNNN,NNN_sp_d1,
      NNN_sp_d2,px_xbond,
      px_ybond,py_xbond,
      py_ybond,Htotal,
      Nx,Ny,orbital_id,phonon_id,
       phonon_coord_try2,Nelec,α,
      β,gshear,temp)
    end

    return free_energy_symmetric,diff_try, symmetrized_phonon_coor


end





function findFL(Nelec::Float64,spectrum::Vector{Float64},temp::Float64,val_s::Float64,val_e::Float64)
  fl=0
  stan=10^(-7)

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