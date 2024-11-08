using LinearAlgebra



function get_Moire_Ham(wave::Vector{Vector{Int}},wAA::Float64,wAB::Float64,vF::Float64,qset::Vector{Vector{Float64}},Kset::Vector{Vector{Float64}},dt::Vector{Float64},db::Vector{Float64},T1::Vector{Float64},T2::Vector{Float64},vi::Int,kvec::Vector{Float64},g1mT::Vector{Int},g2mT::Vector{Int})
 
  num_layer=3
  num_sub=2
  

  σx=[0.0 1.0;1.0 0.0]
  σy=[0.0 -im;im 0.0]

 Moire=zeros(ComplexF64,num_layer,num_sub,length(wave),num_layer,num_sub,length(wave))
 if vi==1
   for jb in eachindex(wave)

     Moire[1,:,jb,2,:,jb]+=[wAA wAB;wAB wAA]*exp(-im*dot(qset[1],dt))
     Moire[2,:,jb,3,:,jb]+=[wAA wAB;wAB wAA]*exp(-im*dot(qset[1],db))


    pos=findfirst(item->item==wave[jb]+g1mT,wave)
    if pos≠nothing
       Moire[1,:,pos,2,:,jb]+=[wAA wAB*exp(-im*2*π/3);wAB*exp(im*2*π/3) wAA]*exp(-im*dot(qset[2],dt))
       Moire[2,:,pos,3,:,jb]+=[wAA wAB*exp(-im*2*π/3);wAB*exp(im*2*π/3) wAA]*exp(-im*dot(qset[2],db))
    end

    pos=findfirst(item->item==wave[jb]+g2mT,wave)
    if pos≠nothing
      Moire[1,:,pos,2,:,jb]+=[wAA wAB*exp(im*2*π/3);wAB*exp(-im*2*π/3) wAA]*exp(-im*dot(qset[3],dt))
      Moire[2,:,pos,3,:,jb]+=[wAA wAB*exp(im*2*π/3);wAB*exp(-im*2*π/3) wAA]*exp(-im*dot(qset[3],db))

     end
   end
 end


 if vi==-1
   for jb in eachindex(wave)

      Moire[2,:,jb,1,:,jb]+=[wAA wAB;wAB wAA]*exp(-im*dot(qset[1],dt))
      Moire[3,:,jb,2,:,jb]+=[wAA wAB;wAB wAA]*exp(-im*dot(qset[1],db))


     pos=findfirst(item->item==wave[jb]+g1mT,wave)
     if pos≠nothing
        Moire[2,:,pos,1,:,jb]+=[wAA wAB*exp(-im*2*π/3);wAB*exp(im*2*π/3) wAA]*exp(-im*dot(qset[2],dt))
        Moire[3,:,pos,2,:,jb]+=[wAA wAB*exp(-im*2*π/3);wAB*exp(im*2*π/3) wAA]*exp(-im*dot(qset[2],db))
     end

     pos=findfirst(item->item==wave[jb]+g2mT,wave)
     if pos≠nothing
       Moire[2,:,pos,1,:,jb]+=[wAA wAB*exp(im*2*π/3);wAB*exp(-im*2*π/3) wAA]*exp(-im*dot(qset[3],dt))
       Moire[3,:,pos,2,:,jb]+=[wAA wAB*exp(im*2*π/3);wAB*exp(-im*2*π/3) wAA]*exp(-im*dot(qset[3],db))

      end
    end
 end

 Hamiltonian=zeros(ComplexF64,num_layer,num_sub,length(wave),num_layer,num_sub,length(wave))
 
 for jb in eachindex(wave)
   kvec1=kvec+wave[jb][1]*T1+wave[jb][2]*T2-vi*Kset[1]
   Hamiltonian[1,:,jb,1,:,jb]+=vF*(σx*kvec1[1]+σy*kvec1[2])
   kvec2=kvec+wave[jb][1]*T1+wave[jb][2]*T2-vi*Kset[2]
   Hamiltonian[2,:,jb,2,:,jb]+=vF*(σx*kvec2[1]+σy*kvec2[2])
   kvec3=kvec+wave[jb][1]*T1+wave[jb][2]*T2-vi*Kset[3]
   Hamiltonian[3,:,jb,3,:,jb]+=vF*(σx*kvec3[1]+σy*kvec3[2])
 end


 return Moire, Hamiltonian


end


function Coulomb(k::Vector{Int},T1::Vector{Float64},T2::Vector{Float64})::Float64
    D=25
    return k==[0,0] ? D : tanh(norm(k[1]*T1+k[2]*T2)*D)/norm(k[1]*T1+k[2]*T2)
end


function sendtomesh(Minv::Matrix{Int64},q1::Vector{Int})::Vector{Int}
    Qvec=q1'*inv(Minv)
    return q1 .-vec((Int.(floor.(round.(Qvec,digits=5)))*Minv)')
end


function shuffle_vector(eg_vec::Matrix{ComplexF64},shuff_vec::Vector{Int},wave::Vector{Vector{Int}},wave_dic::Dict{Vector{Int64},Int64},Nband::Int)::Array{ComplexF64}
    num_layer=3 # shuff_vec is g, I am out putting eigvec*exp(igr)
    num_sub=2
    res_eigvec=reshape(eg_vec,num_layer,num_sub,length(wave),Nband)
    shuff_eigvec=zeros(ComplexF64,num_layer,num_sub,length(wave),Nband)
    for ja in eachindex(wave)
        if haskey(wave_dic,wave[ja]-shuff_vec)
            pos=wave_dic[wave[ja]-shuff_vec]
          shuff_eigvec[:,:,ja,:]=res_eigvec[:,:,pos,:]
        end
    end
   return reshape(shuff_eigvec,num_layer*num_sub*length(wave),Nband)
end


function Geometry(geonum::Int)
  if geonum==1
    l1=[18,0]
    l2=[0,18]
    Nx=18;
    Ny=18;
    
  end

  if geonum==2
    l1=[6,0]
    l2=[0,6]
    Nx=6;
    Ny=6;
    
  end

  return l1,l2,Nx,Ny

end


function single_particle(geonum::Int64,θ::Float64,wAA::Float64,wAB::Float64,vF::Float64,ϵr::Float64,Nband::Int64)


  #a0=0.246
  a0=0.142*√3

  
  KGr=4π/(3*a0)
  kθ=2*KGr*sin(θ/2)
  
  
  Kset=[kθ*[0.0,1.0],[0.0,0.0],kθ*[0.0,-1.0]]
  qset=[kθ*[0.0,-1.0],kθ*[√3/2,1/2],kθ*[-√3/2,1/2]]
  
  
  g1m=kθ*√3*[1/2,√3/2]
  g2m=kθ*√3*[-1/2,√3/2]
  
  
  a1m=4π/(3*kθ)*[√3/2,1/2]
  a2m=4π/(3*kθ)*[-√3/2,1/2]
  db=[0.0,0.0]
  dt=1/3*(a2m-a1m)
  num_sub=2
  num_layer=3
  num_spin=2
  num_valley=2
  
  
  
  println("moirelength",norm(a1m))
  

  l1,l2,Nx,Ny=Geometry(geonum)
  L1=l1[1]*a1m+l1[2]*a2m;
  L2=l2[1]*a1m+l2[2]*a2m;
  area=abs(L1[1]*L2[2]-L2[1]*L1[2]);
  Rotminus90=[0 1;-1 0]
  T1=2*π/area*Rotminus90*L2
  T2=-2*π/area*Rotminus90*L1


  
  
  g1mT=Int.(round.(inv([T1 T2])*g1m))
  g2mT=Int.(round.(inv([T1 T2])*g2m))
  Minv=[g1mT';g2mT']
  
  allowedq=Vector{Int}[]
  for jb in 0:Ny-1, ja in 0:Nx-1
    push!(allowedq,sendtomesh(Minv,[ja,jb]))
  end
  allowedq_dic=Dict{Vector{Int},Int}()
  
  for ja in eachindex(allowedq)
    allowedq_dic[allowedq[ja]]=ja
  end

  #constq=1/(√3/2*norm(a1m)^2*ϵr*length(allowedq))*9047.5636
  constq=1/(√3/2*norm(a1m)^2*ϵr*length(allowedq))*9035.4642
  wave=Vector{Int64}[]
  cutoff=18
  cutoffstandard=3.01*norm(g1m)
  for ja in -cutoff:cutoff, jb in -cutoff:cutoff
      gtest=ja*g1m+jb*g2m;
      if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
          push!(wave,ja*g1mT+jb*g2mT)
      end
  end
  wave_dic=Dict{Vector{Int},Int}()
  for ja in eachindex(wave)
    wave_dic[wave[ja]]=ja
  end
  
  allowedq_dic=Dict{Vector{Int},Int}()
  for ja in eachindex(allowedq)
    allowedq_dic[allowedq[ja]]=ja
  end
  
  
  dimension=num_layer*length(wave)*num_sub
  wave_diff=Vector{Int64}[]
  cutoff=18
  cutoffstandard=6.01*norm(g1m)
  for ja in -cutoff:cutoff, jb in -cutoff:cutoff
      gtest=ja*g1m+jb*g2m;
      if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
          push!(wave_diff,ja*g1mT+jb*g2mT)
      end
  end

 
  eigenvalue=[zeros(Float64,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in 1:Nx*Ny] 
  eigenvector=[zeros(ComplexF64,dimension,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in 1:Nx*Ny]
  
  
  
  for ja in eachindex(allowedq), spin_i in 1:num_spin, valley in  1:num_spin
      vset=[1,-1]
      kvec=allowedq[ja][1]*T1+allowedq[ja][2]*T2
      Moire, Ham=get_Moire_Ham(wave,wAA,wAB,vF,qset,Kset,dt,db,T1,T2,vset[valley],kvec,g1mT,g2mT)
      totalHam=reshape(Moire,dimension,dimension)+reshape(Moire,dimension,dimension)'+reshape(Ham,dimension,dimension)
      FFF=eigen(totalHam)
  
      eigenvalue[spin_i,valley,ja]=real.(FFF.values[Int(dimension/2):Int(dimension/2)+Nband-1])
      eigenvector[spin_i,valley,ja]=FFF.vectors[:,Int(dimension/2):Int(dimension/2)+Nband-1]
  end
  
  

 return  eigenvector,eigenvalue,wave,wave_diff,wave_dic,allowedq,allowedq_dic,T1,T2,constq,Minv,g1mT,g2mT
  
  

end


function get_initial_proj(allowedq::Vector{Vector{Int}},eigenvalue::Array{Vector{Float64}},Nband::Int)
  num_spin=2
  num_valley=2

  initial_projector=[zeros(ComplexF64,Nband,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)]
  bg_projector=[zeros(ComplexF64,Nband,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)]


 single_Ham=[zeros(ComplexF64,Nband,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)]
 single_eigenvector=[zeros(ComplexF64,Nband,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)]


 for spin_i in 1:num_spin, valley in 1:num_valley,ja in eachindex(allowedq)
  single_Ham[spin_i,valley,ja]=diagm(eigenvalue[spin_i,valley,ja])
  FFF=eigen(single_Ham[spin_i,valley,ja])
  single_eigenvector=FFF.vectors
 end



 for spin_i in 1:num_spin, valley in 1:num_valley, jc in eachindex(allowedq)
   A=randn(Nband,Nband)+im*randn(Nband,Nband)
   initial_projector[spin_i,valley,jc]+=(A+A')*1.0
   bg_projector[spin_i,valley,jc]+=1/2*Matrix{ComplexF64}(I,Nband,Nband)
 end 

 return initial_projector, bg_projector, single_Ham
end




function iteration(formfactors::Array{Matrix{ComplexF64}},initial_projector::Array{Matrix{ComplexF64}},bg_projector::Array{Matrix{ComplexF64}},constq::Float64,Nband::Int,wave_diff::Vector{Vector{Int}},allowedq::Vector{Vector{Int}},allowedq_dic::Dict{Vector{Int},Int},T1::Vector{Float64},T2::Vector{Float64},single_Ham::Array{Matrix{ComplexF64}},Npa::Int)
  num_spin=2
  num_valley=2

  eout=1.0
  itcount=0
  bad_count=0
  energy=0.0
  energy_change=0.0
  bound=0.0
  HF_eigenvalue=[zeros(ComplexF64,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in 1:length(allowedq)]
  HF_eigenvector=[zeros(ComplexF64,Nband,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in 1:length(allowedq)]
  DIIS_input_projector=Vector{Array{Matrix{ComplexF64}}}(undef,3)
  DIIS_input_DeltaMatrix=Vector{Array{Matrix{ComplexF64}}}(undef,3)
  input_projector=initial_projector
  
  
  
  while (eout>1*10^-15) || (bad_count<4) || (energy_change>1*10^-8)
      if eout<1*10^-15
       bad_count+=1
      end
      tic=time()
      eout,energy_change,output_projector,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue,bound,HF_eigenvector,energy=Construct_projector(formfactors,input_projector,bg_projector,constq,Nband,wave_diff,allowedq,allowedq_dic,T1,T2,single_Ham,energy,Npa)
      DIIS_input_projector[mod(itcount,3)+1]=input_projector
      input_projector=output_projector
      itcount+=1
     
      toc=time()
      println(toc-tic,"eout=$eout","energy_change=$energy_change")
      flush(stdout)
     
    
    end
  
  return HF_eigenvalue,HF_eigenvector,energy, DIIS_input_projector,bound


end






function get_formfactors(allowedq::Vector{Vector{Int}},wave::Vector{Vector{Int}},wave_diff::Vector{Vector{Int}},wave_dic::Dict{Vector{Int64},Int64},Minv::Matrix{Int},Nband::Int,eigenvector::Array{Matrix{ComplexF64}})::Array{Matrix{ComplexF64}}
 num_spin=2
 num_valley=2

  formfactors=Array{Matrix{ComplexF64}}(undef,num_spin,num_valley,length(allowedq),length(allowedq),length(wave_diff)) #The last two are q
  formfactors_threaded=[Array{Matrix{ComplexF64}}(undef,num_spin,num_valley,length(allowedq),length(wave_diff)) for _ in 1:length(allowedq)] 
#=
  for ja in eachindex(allowedq), valley in  1:num_spin
      for qvec in eachindex(allowedq), gqindex in eachindex(wave_diff)
          kplusq_pos=allowedq_dic[sendtomesh(Minv,allowedq[qvec]+allowedq[ja])]
          gkplusq=allowedq[ja]+allowedq[qvec]+wave_diff[gqindex]-allowedq[kplusq_pos]
  
          prod=shuffle_vector(eigenvector[1,valley,ja],gkplusq,wave,wave_dic,Nband)'*eigenvector[1,valley,kplusq_pos]
          formfactors[:, valley, ja, qvec, gqindex]=[prod for _ in 1:num_spin]
    
      end
  end
=#
 Threads.@threads for ja in eachindex(allowedq)
 for  valley in  1:num_valley, spin_i in 1:num_spin
   for qvec in eachindex(allowedq), gqindex in eachindex(wave_diff)
      kplusq_pos=allowedq_dic[sendtomesh(Minv,allowedq[qvec]+allowedq[ja])]
      gkplusq=allowedq[ja]+allowedq[qvec]+wave_diff[gqindex]-allowedq[kplusq_pos]
      s1=shuffle_vector(eigenvector[spin_i,valley,ja],gkplusq,wave,wave_dic,Nband)
      prod=s1'*eigenvector[spin_i,valley,kplusq_pos]
      formfactors_threaded[ja][spin_i, valley, qvec, gqindex]=prod

    end
   end
 end

 for ja in eachindex(allowedq)
  formfactors[:, :, ja, :, :]=formfactors_threaded[ja]
 end

  return formfactors

end


function Construct_projector(form_factor::Array{Matrix{ComplexF64}},projector::Array{Matrix{ComplexF64}},bg_projector::Array{Matrix{ComplexF64}},constq::Float64,Nband::Int64,wave_diff::Vector{Vector{Int64}},allowedq::Vector{Vector{Int}},allowedq_dic::Dict{Vector{Int},Int},T1::Vector{Float64},T2::Vector{Float64},single_Ham::Array{Matrix{ComplexF64}},energy_input::Float64,Npa::Int64)
  num_spin=2
  num_valley=2
  
  
 # Hartree=[zeros(ComplexF64,Nband,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)]
 Hartree_threaded=[[zeros(ComplexF64,Nband,Nband) for _ in 1:num_spin,_ in 1:num_valley] for _ in eachindex(allowedq)]
 # Fock=[zeros(ComplexF64,Nband,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)]
  Fock_threaded=[[zeros(ComplexF64,Nband,Nband) for _ in 1:num_spin,_ in 1:num_valley] for _ in eachindex(allowedq)]
  New_projector=[zeros(ComplexF64,Nband,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)]
  
Threads.@threads for k2 in eachindex(allowedq)
  for jqmesh in eachindex(allowedq)
   k1pos=allowedq_dic[sendtomesh(Minv,allowedq[k2]+allowedq[jqmesh])]
   for jqg in eachindex(wave_diff), spin_i in 1:num_spin, valley in 1:num_valley
     Fock_threaded[k2][spin_i,valley]+=Coulomb(allowedq[jqmesh]+wave_diff[jqg],T1,T2)*form_factor[spin_i,valley,k2,jqmesh,jqg]*projector[spin_i,valley,k1pos]*(form_factor[spin_i,valley,k2,jqmesh,jqg])'
   end
  end
 end 
 
 HartreeDensity=zeros(ComplexF64,length(wave_diff))
 zeropos=allowedq_dic[[0,0]]
 for jqg in eachindex(wave_diff), jb in eachindex(allowedq),spin_i in 1:num_spin, valley in 1:num_valley
    HartreeDensity[jqg]+=tr(projector[spin_i,valley,jb]*(form_factor[spin_i,valley,jb,zeropos,jqg])')
 end

 Threads.@threads for ja in eachindex(allowedq)
 for spin_i in 1:num_spin, valley in 1:num_valley, jqg in eachindex(wave_diff)
      Hartree_threaded[ja][spin_i,valley]+=Coulomb(wave_diff[jqg],T1,T2)*form_factor[spin_i,valley,ja,zeropos,jqg]*HartreeDensity[jqg]
 end
 end

 HF_eigenvalue=[zeros(ComplexF64,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)]
 HF_eigenvector=[zeros(ComplexF64,Nband,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)]

 for spin_i in 1:num_spin, valley in 1:num_valley,ja in eachindex(allowedq)
    FFF=eigen(constq*Hartree_threaded[ja][spin_i,valley]+single_Ham[spin_i,valley,ja]-constq*Fock_threaded[ja][spin_i,valley])
    HF_eigenvalue[spin_i,valley,ja]=(FFF.values)
    HF_eigenvector[spin_i,valley,ja]=FFF.vectors
 end



 sorted=sort(reduce(vcat,reduce(vcat,real.(HF_eigenvalue)))) 
 bound=(sorted[Npa+1]+sorted[Npa])/2

  for spin_i in 1:num_spin, valley in 1:num_valley,ja in eachindex(allowedq)
       for jd in eachindex(HF_eigenvalue[spin_i,valley,ja])
          if real(HF_eigenvalue[spin_i,valley,ja][jd])<bound
             New_projector[spin_i,valley,ja]+=HF_eigenvector[spin_i,valley,ja][:,jd]*(HF_eigenvector[spin_i,valley,ja][:,jd])'
          end
       end
  end

 energy=0.0
  for spin_i in 1:num_spin, valley in 1:num_valley,ja in eachindex(allowedq)
    energy+=tr((constq/2*Hartree_threaded[ja][spin_i,valley]+single_Ham[spin_i,valley,ja]-constq/2*Fock_threaded[ja][spin_i,valley])*New_projector[spin_i,valley,ja])
  end
  for spin_i in 1:num_spin, valley in 1:num_valley,ja in eachindex(allowedq)
    energy-=tr((constq/2*Hartree_threaded[ja][spin_i,valley]-constq/2*Fock_threaded[ja][spin_i,valley])*bg_projector[spin_i,valley,ja])
  end

  New_projector-=bg_projector

 

  DeltaMatrix=New_projector-projector
  output_projector=0.0*projector+1.0*New_projector

  e1=0.0
  for spin_i in 1:num_spin, valley in 1:num_valley,ja in eachindex(allowedq)
    e1+=tr(DeltaMatrix[spin_i,valley,ja]'*DeltaMatrix[spin_i,valley,ja])
  end
  eout=real(e1)/length(allowedq)
   
 #energy=0.0
 #for spin_i in 1:num_spin, valley in 1:num_valley,ja in eachindex(allowedq)
    #energy+=tr((constq/2*Hartree_threaded[ja][spin_i,valley]+single_Ham[spin_i,valley,ja]-constq/2*Fock_threaded[ja][spin_i,valley])*output_projector[spin_i,valley,ja])
 #end


 energy_change=real(energy-energy_input)

  return  eout,energy_change,output_projector,DeltaMatrix,HF_eigenvalue,bound, HF_eigenvector,real(energy)
 

end