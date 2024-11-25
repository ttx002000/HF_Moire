using LinearAlgebra



function get_Moire_Ham(wave::Vector{Vector{Int}},wAA::Float64,wAB::Float64,vF::Float64,qset::Vector{Vector{Float64}},Kset::Vector{Vector{Float64}},dt::Vector{Float64},db::Vector{Float64},T1::Vector{Float64},T2::Vector{Float64},vi::Int,kvec::Vector{Float64},g1m_ps_int::Vector{Int},g2m_ps_int::Vector{Int},lambda_MDT::Float64,Dfield::Float64)
 
  num_layer=3
  num_sub=2
  
  qper_set=[[1.0,0.0],[-1/2,√3/2],[-1/2,-√3/2]]

  σx=[0.0 1.0;1.0 0.0]
  σy=[0.0 -im;im 0.0]

 Moire=zeros(ComplexF64,num_layer,num_sub,length(wave),num_layer,num_sub,length(wave))
 if vi==1
  Moire=zeros(ComplexF64,num_layer,num_sub,length(wave),num_layer,num_sub,length(wave))
   
    for jb in eachindex(wave)
     absolute_k=kvec+wave[jb][1]*T1+wave[jb][2]*T2
     Moire[2,:,jb,1,:,jb]+=[wAA wAB;wAB wAA]*exp(-im*dot(qset[1],dt))*(1+lambda_MDT*dot(qper_set[1],absolute_k-Kset[1]))
 

     pos=findfirst(item->item==wave[jb]+g1m_ps_int,wave)
     if pos≠nothing
       Moire[2,:,pos,1,:,jb]+=[wAA wAB*exp(-im*2*π/3);wAB*exp(im*2*π/3) wAA]*exp(-im*dot(qset[2],dt))*(1+lambda_MDT*dot(qper_set[2],absolute_k-Kset[1]))
     end

     pos=findfirst(item->item==wave[jb]+g2m_ps_int,wave)
     if pos≠nothing
       Moire[2,:,pos,1,:,jb]+=[wAA wAB*exp(im*2*π/3);wAB*exp(-im*2*π/3) wAA]*exp(-im*dot(qset[3],dt))*(1+lambda_MDT*dot(qper_set[3],absolute_k-Kset[1]))
     end
   end


   for jb in eachindex(wave)

     absolute_k=kvec+wave[jb][1]*T1+wave[jb][2]*T2
      Moire[2,:,jb,3,:,jb]+=[wAA wAB;wAB wAA]*exp(-im*2*dot(qset[1],db))*(1+lambda_MDT*dot(qper_set[1],absolute_k-Kset[3]))


      pos=findfirst(item->item==wave[jb]+2*g1m_ps_int,wave)
     if pos≠nothing
         Moire[2,:,pos,3,:,jb]+=[wAA wAB*exp(-im*2*π/3);wAB*exp(im*2*π/3) wAA]*exp(-im*2*dot(qset[2],db))*(1+lambda_MDT*dot(qper_set[2],absolute_k-Kset[3]))
     end

      pos=findfirst(item->item==wave[jb]+2*g2m_ps_int,wave)
      if pos≠nothing
         Moire[2,:,pos,3,:,jb]+=[wAA wAB*exp(im*2*π/3);wAB*exp(-im*2*π/3) wAA]*exp(-im*2*dot(qset[3],db))*(1+lambda_MDT*dot(qper_set[3],absolute_k-Kset[3]))
     end
   end
 end
    

 if vi==-1
  Moire=zeros(ComplexF64,num_layer,num_sub,length(wave),num_layer,num_sub,length(wave))
   
    for jb in eachindex(wave)
     absolute_k=kvec+wave[jb][1]*T1+wave[jb][2]*T2
     Moire[1,:,jb,2,:,jb]+=[wAA wAB;wAB wAA]*exp(-im*dot(qset[1],dt))*(1+lambda_MDT*dot(qper_set[1],absolute_k+Kset[2]))
 

     pos=findfirst(item->item==wave[jb]+g1m_ps_int,wave)
     if pos≠nothing
       Moire[1,:,pos,2,:,jb]+=[wAA wAB*exp(-im*2*π/3);wAB*exp(im*2*π/3) wAA]*exp(-im*dot(qset[2],dt))*(1+lambda_MDT*dot(qper_set[2],absolute_k+Kset[2]))
     end

     pos=findfirst(item->item==wave[jb]+g2m_ps_int,wave)
     if pos≠nothing
       Moire[1,:,pos,2,:,jb]+=[wAA wAB*exp(im*2*π/3);wAB*exp(-im*2*π/3) wAA]*exp(-im*dot(qset[3],dt))*(1+lambda_MDT*dot(qper_set[3],absolute_k+Kset[2]))
     end
   end


   for jb in eachindex(wave)

     absolute_k=kvec+wave[jb][1]*T1+wave[jb][2]*T2
      Moire[3,:,jb,2,:,jb]+=[wAA wAB;wAB wAA]*exp(-im*2*dot(qset[1],db))*(1+lambda_MDT*dot(qper_set[1],absolute_k+Kset[2]))


      pos=findfirst(item->item==wave[jb]+2*g1m_ps_int,wave)
     if pos≠nothing
         Moire[3,:,pos,2,:,jb]+=[wAA wAB*exp(-im*2*π/3);wAB*exp(im*2*π/3) wAA]*exp(-im*2*dot(qset[2],db))*(1+lambda_MDT*dot(qper_set[2],absolute_k+Kset[2]))
     end

      pos=findfirst(item->item==wave[jb]+2*g2m_ps_int,wave)
      if pos≠nothing
         Moire[3,:,pos,2,:,jb]+=[wAA wAB*exp(im*2*π/3);wAB*exp(-im*2*π/3) wAA]*exp(-im*2*dot(qset[3],db))*(1+lambda_MDT*dot(qper_set[3],absolute_k+Kset[2]))
     end
   end
 end



 Hamiltonian=zeros(ComplexF64,num_layer,num_sub,length(wave),num_layer,num_sub,length(wave))
 
 for jb in eachindex(wave)
   kvec1=kvec+wave[jb][1]*T1+wave[jb][2]*T2-vi*Kset[1]
   Hamiltonian[1,:,jb,1,:,jb]+=vF*(σx*kvec1[1]+σy*kvec1[2])-Dfield*Matrix{Float64}(I,2,2)
   kvec2=kvec+wave[jb][1]*T1+wave[jb][2]*T2-vi*Kset[2]
   Hamiltonian[2,:,jb,2,:,jb]+=vF*(σx*kvec2[1]+σy*kvec2[2])
   kvec3=kvec+wave[jb][1]*T1+wave[jb][2]*T2-vi*Kset[3]
   Hamiltonian[3,:,jb,3,:,jb]+=vF*(σx*kvec3[1]+σy*kvec3[2])+Dfield*Matrix{Float64}(I,2,2)
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


function shuffle_vector(eg_vec::Matrix{ComplexF64},shuff_vec::Vector{Int},wave::Vector{Vector{Int}},wave_dic::Dict{Vector{Int64},Int64},Nband::Int)::Matrix{ComplexF64}
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



function shuffle_vector_singlevector(eg_vec::Vector{ComplexF64},shuff_vec::Vector{Int},wave::Vector{Vector{Int}},wave_dic::Dict{Vector{Int64},Int64})::Vector{ComplexF64}
  num_layer=3 # shuff_vec is g, I am out putting eigvec*exp(igr)
  num_sub=2
  res_eigvec=reshape(eg_vec,num_layer,num_sub,length(wave))
  shuff_eigvec=zeros(ComplexF64,num_layer,num_sub,length(wave))
  for ja in eachindex(wave)
      if haskey(wave_dic,wave[ja]-shuff_vec)
          pos=wave_dic[wave[ja]-shuff_vec]
        shuff_eigvec[:,:,ja]=res_eigvec[:,:,pos]
      end
  end
 return reshape(shuff_eigvec,num_layer*num_sub*length(wave))
end


function Geometry(geonum::Int)
  if geonum==1
    l1=[6,0]
    l2=[0,6]
    Nx=6;
    Ny=6;
    
  end

  if geonum==2
    l1=[9,0]
    l2=[0,9]
    Nx=9;
    Ny=9;
    
  end



  return l1,l2,Nx,Ny

end


function single_particle(geonum::Int64,θ::Float64,wAA::Float64,wAB::Float64,vF::Float64,ϵr::Float64,Nband::Int64,lambda_MDT::Float64,Nb_down::Int,Nb_up::Int,Dfield::Float64,shift::Int64)


  a0=0.246
  

  
  KGr=4π/(3*a0)
  kθ=2*KGr*sin(θ/2)
  
  
  Kset=[kθ*[0.0,0.0],kθ*[0.0,1.0],kθ*[0.0,-1.0]]
  qset=[kθ*[0.0,-1.0],kθ*[√3/2,1/2],kθ*[-√3/2,1/2]]
  
  
  g1m_ps=kθ*√3*[1/2,√3/2] #stands for pristine
  g2m_ps=kθ*√3*[-1/2,√3/2]
  
  
  a1m_ps=4π/(3*kθ)*[√3/2,1/2]
  a2m_ps=4π/(3*kθ)*[-√3/2,1/2]
  



  if shift==1
    db=[0.0,0.0]
    dt=1/3*(a2m_ps-a1m_ps)
    gridshift=[0.0,0.0]
  elseif  shift==2
    db=[0.0,0.0]
    dt=[0.0,0.0]
    gridshift=[0.0,0.0]
  elseif  shift==3
    db=[0.0,0.0]
    dt=1/2*(a2m_ps)
    gridshift=[0.0,0.0]
  end

  a1m=2*a1m_ps+a2m_ps
  a2m=2*a2m_ps+a1m_ps
  #a1m=a1m_ps
  #a2m=a2m_ps
  CC=2π*inv([a1m a2m])
  g1m=CC[1,:]
  g2m=CC[2,:]
  num_sub=2
  num_layer=3
  num_spin=2
  num_valley=2
  
  
  

  

  l1,l2,Nx,Ny=Geometry(geonum)
  L1=l1[1]*a1m+l1[2]*a2m;
  L2=l2[1]*a1m+l2[2]*a2m;
  area=abs(L1[1]*L2[2]-L2[1]*L1[2]);
  Rotminus90=[0 1;-1 0]
  T1=2*π/area*Rotminus90*L2
  T2=-2*π/area*Rotminus90*L1


  
  
  g1mT=Int.(round.(inv([T1 T2])*g1m))
  g2mT=Int.(round.(inv([T1 T2])*g2m))

  g1m_ps_int=Int.(round.(inv([T1 T2])*g1m_ps))
  g2m_ps_int=Int.(round.(inv([T1 T2])*g2m_ps))

  Minv=[g1mT';g2mT']
  
  allowedq=Vector{Int}[]
  for jb in 0:Ny-1, ja in 0:Nx-1
    push!(allowedq,sendtomesh(Minv,[ja,jb]))
  end
  allowedq_dic=Dict{Vector{Int},Int}()
  
  for ja in eachindex(allowedq)
    allowedq_dic[allowedq[ja]]=ja
  end

  constq=1/(√3/2*norm(a1m)^2*ϵr*length(allowedq))*9047.5636
 #constq=0.0
  wave=Vector{Int64}[]
  cutoff=20*5
  g_cutoff=6.01
  cutoffstandard=g_cutoff*norm(g1m_ps)
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
  cutoff=20*5
  q_cutoff=9.01
  cutoffstandard=q_cutoff*norm(g1m_ps)
  for ja in -cutoff:cutoff, jb in -cutoff:cutoff
      gtest=ja*g1m+jb*g2m;
      if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
          push!(wave_diff,ja*g1mT+jb*g2mT)
      end
  end

 
  eigenvalue=[zeros(Float64,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)] 
  eigenvector=[zeros(ComplexF64,dimension,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)]
  
  
  
  Threads.@threads for ja in eachindex(allowedq)
    for spin_i in 1:num_spin, valley in  1:num_spin
      vset=[1,-1]
      kvec=allowedq[ja][1]*T1+allowedq[ja][2]*T2+gridshift
      Moire, Ham=get_Moire_Ham(wave,wAA,wAB,vF,qset,Kset,dt,db,T1,T2,vset[valley],kvec,g1m_ps_int, g2m_ps_int,lambda_MDT,Dfield)
      totalHam=reshape(Moire,dimension,dimension)+reshape(Moire,dimension,dimension)'+reshape(Ham,dimension,dimension)
      FFF=eigen(totalHam)
  
      eigenvalue[spin_i,valley,ja]=real.(FFF.values[Int(dimension/2)-Nb_down+1:Int(dimension/2)+Nb_up])
      eigenvector[spin_i,valley,ja]=FFF.vectors[:,Int(dimension/2)-Nb_down+1:Int(dimension/2)+Nb_up]
    end
  end
  
  

 return  eigenvector,eigenvalue,wave,wave_diff,wave_dic,allowedq,allowedq_dic,T1,T2,constq,Minv,g1mT,g2mT,a1m,a2m,g_cutoff,q_cutoff
  
  

end


function  get_bias(g1mT,g2mT,T1,T2,eigenvector,allowedq,a1m,a2m,Nband)
  num_layer=3
  num_sub=2
  num_valley=2
  num_spin=2
  sub_pick=rand([1,2]) 
  valley_pick=rand([1,2])
  spin_pick=rand([1,2])
  phi=rand([2/3*π,4/3*π,0.0])
  g1m=g1mT[1]*T1+g1mT[2]*T2
  g2m=g2mT[1]*T1+g2mT[2]*T2

    perturb=zeros(ComplexF64,num_layer,num_sub,length(wave),num_layer,num_sub,length(wave))
    for jb in eachindex(wave)
 
      pos=findfirst(item->item==wave[jb]-g1mT-g2mT,wave)
      if pos≠nothing
        perturb[:,sub_pick,pos,:,sub_pick,jb]+=5*Matrix{Float64}(I,num_layer,num_layer)*exp(im*phi)*exp(-im*dot(-g1m-g2m,1/3*a1m))
      end
  
 
      pos=findfirst(item->item==wave[jb]+g1mT,wave)
      if pos≠nothing
        perturb[:,sub_pick,pos,:,sub_pick,jb]+=5*Matrix{Float64}(I,num_layer,num_layer)*exp(im*phi)*exp(-im*dot(g1m,1/3*a1m))
      end
 
      pos=findfirst(item->item==wave[jb]+g2mT,wave)
      if pos≠nothing
        perturb[:,sub_pick,pos,:,sub_pick,jb]+=5*Matrix{Float64}(I,num_layer,num_layer)*exp(im*phi)*exp(-im*dot(g2m,1/3*a1m))
      end
    end
     reshaped_perturb=reshape(perturb,num_layer*num_sub*length(wave),num_layer*num_sub*length(wave))+reshape(perturb,num_layer*num_sub*length(wave),num_layer*num_sub*length(wave))'

    perturb_Ham=[zeros(ComplexF64,Nband,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)]
    for ja in eachindex(allowedq)
       perturb_Ham[spin_pick,valley_pick,ja]=eigenvector[spin_pick,valley_pick,ja]'*reshaped_perturb*eigenvector[spin_pick,valley_pick,ja]
    end

    if only(rand(1))>0.4
      return  [zeros(ComplexF64,Nband,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)]
    else
      return  perturb_Ham
    end


end


function get_initial_proj(allowedq::Vector{Vector{Int}},eigenvalue::Array{Vector{Float64}},Nband::Int,perturb_Ham::Array{Matrix{ComplexF64}},Npa::Int)
  num_spin=2
  num_valley=2

  initial_projector=[zeros(ComplexF64,Nband,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)]
  bg_projector=[1/2*Matrix{ComplexF64}(I,Nband,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)]


 single_Ham=[zeros(ComplexF64,Nband,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)]
 single_eigenvector=[zeros(ComplexF64,Nband,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)]

 single_eigenvalue=[zeros(ComplexF64,Nband) for _ in 1:num_spin,_ in 1:num_valley, _ in eachindex(allowedq)]

 for spin_i in 1:num_spin, valley in 1:num_valley,ja in eachindex(allowedq)
  single_Ham[spin_i,valley,ja]=diagm(eigenvalue[spin_i,valley,ja])
  FFF=eigen(single_Ham[spin_i,valley,ja]+perturb_Ham[spin_i,valley,ja])
  single_eigenvector[spin_i,valley,ja]=FFF.vectors
  single_eigenvalue[spin_i,valley,ja]=FFF.values
 end



    
 sorted=sort(reduce(vcat,reduce(vcat,real.(single_eigenvalue)))) 
 bound=(sorted[Npa+1]+sorted[Npa])/2

 Threads.@threads for ja in eachindex(allowedq)
  for spin_i in 1:num_spin, valley in 1:num_valley
       for jd in eachindex(single_eigenvalue[spin_i,valley,ja])
          if real(single_eigenvalue[spin_i,valley,ja][jd])<bound
            initial_projector[spin_i,valley,ja]+=single_eigenvector[spin_i,valley,ja][:,jd]*(single_eigenvector[spin_i,valley,ja][:,jd])'
          end
       end
  end
 end




 for spin_i in 1:num_spin, valley in 1:num_valley, jc in eachindex(allowedq)
   A=randn(Nband,Nband)+im*randn(Nband,Nband)
   initial_projector[spin_i,valley,jc]+=(A+A')*0.3
 end 

 return initial_projector, bg_projector, single_Ham
end




function iteration(formfactors::Array{Matrix{ComplexF64}},initial_projector::Array{Matrix{ComplexF64}},bg_projector::Array{Matrix{ComplexF64}},constq::Float64,Nband::Int,wave_diff::Vector{Vector{Int}},allowedq::Vector{Vector{Int}},allowedq_dic::Dict{Vector{Int},Int},T1::Vector{Float64},T2::Vector{Float64},single_Ham::Array{Matrix{ComplexF64}},perturb_Ham::Array{Matrix{ComplexF64}},Npa::Int)
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
  
  
  while itcount<20
  
    
    tic=time()

      eout,energy_change,output_projector,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue,bound,HF_eigenvector,energy=Construct_projector(formfactors,input_projector,bg_projector,constq,Nband,wave_diff,allowedq,allowedq_dic,T1,T2,single_Ham+perturb_Ham,energy,Npa)
      DIIS_input_projector[mod(itcount,3)+1]=input_projector
      input_projector=output_projector 

    itcount+=1
   
    toc=time()
    println(toc-tic,"eout=$eout","energy_change=$energy_change")
    flush(stdout)
   
  
 end


  itcount=0







  
  while (eout>1*10^(-12)) || (bad_count<4) || (energy_change>1*10^(-6))
      if eout<1*10^(-12)
       bad_count+=1
      end
      
      tic=time()
      if (itcount>30 && abs(energy_change)>0.1) || (itcount>30 && abs(eout)<10^(-8))

        dmk=implement_DIIS(DIIS_input_projector,DIIS_input_DeltaMatrix,allowedq)
        eout,energy_change,output_projector,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue,bound,HF_eigenvector,energy=Construct_projector(formfactors,dmk,bg_projector,constq,Nband,wave_diff,allowedq,allowedq_dic,T1,T2,single_Ham,energy,Npa)
        DIIS_input_projector[mod(itcount,3)+1]=dmk
        input_projector=output_projector
        println("using DIIS")
       
      else
        eout,energy_change,output_projector,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue,bound,HF_eigenvector,energy=Construct_projector(formfactors,input_projector,bg_projector,constq,Nband,wave_diff,allowedq,allowedq_dic,T1,T2,single_Ham,energy,Npa)
        DIIS_input_projector[mod(itcount,3)+1]=input_projector
        input_projector=output_projector 
         
     

      end

    



      itcount+=1
     
      toc=time()
      println(toc-tic,"eout=$eout","energy_change=$energy_change")
      flush(stdout)
     
    
  end
 
 


  
  return HF_eigenvalue,HF_eigenvector,energy, DIIS_input_projector,bound


end

function implement_DIIS(DIIS_input_projector,DIIS_input_DeltaMatrix,allowedq)

  num_spin=2
  num_valley=2

      Bmatrix=zeros(ComplexF64,4,4)
      for ja in 1:3
       Bmatrix[ja,4]=1
       Bmatrix[4,ja]=1
      end
  
      for ja in 1:3,jb in 1:3
          for spin_i in 1:num_spin, valley in 1:num_valley,jc in eachindex(allowedq)
             Bmatrix[ja,jb]+=tr((DIIS_input_DeltaMatrix[ja][spin_i,valley,jc])'*(DIIS_input_DeltaMatrix[jb][spin_i,valley,jc]))
          end
      end
      coeff=inv(Bmatrix)*[0;0;0;1]
     
      dmk=coeff[1]*(DIIS_input_projector[1]+DIIS_input_DeltaMatrix[1])+coeff[2]*(DIIS_input_projector[2]+DIIS_input_DeltaMatrix[2])+coeff[3]*(DIIS_input_projector[3]+DIIS_input_DeltaMatrix[3])

    return dmk
end




function get_polarization(HF_eigenvector,eigenvector,allowedq,wave,Nband)
  num_spin=2
  num_valley=2
  num_layer=3
  num_sub=2
  PW_basis_eig=[zeros(ComplexF64,num_layer,num_sub,length(wave),Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)]
 
  for spin_i in 1:num_spin, valley in 1:num_valley, ja in eachindex(allowedq)
   PW_basis_eig[spin_i,valley,ja]=reshape(eigenvector[spin_i,valley,ja]*HF_eigenvector[spin_i,valley,ja],num_layer,num_sub,length(wave),Nband)
  end

  layer_pol=[zeros(ComplexF64,num_layer,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)]
  sub_pol=[zeros(ComplexF64,num_sub,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)]
 
  for spin_i in 1:num_spin, valley in 1:num_valley, ja in eachindex(allowedq)
    for layer_index in 1:num_layer, bandi in 1:Nband
    layer_pol[spin_i,valley,ja][layer_index,bandi]=sum(abs.(vec(PW_basis_eig[spin_i,valley,ja][layer_index,:,:,bandi])).^2)
    end
  end

  for spin_i in 1:num_spin, valley in 1:num_valley, ja in eachindex(allowedq)
    for sub_index in 1:num_sub, bandi in 1:Nband
      sub_pol[spin_i,valley,ja][sub_index,bandi]=sum(abs.(vec(PW_basis_eig[spin_i,valley,ja][:,sub_index ,:,bandi])).^2)
    end
  end

  return layer_pol,sub_pol


end



function get_chernsub(HF_eigenvector,eigenvector,allowedq,wave,Nband)
  num_spin=2
  num_valley=2
  num_layer=3
  num_sub=2
  PW_basis_eig=[zeros(ComplexF64,num_layer*num_sub*length(wave),Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)]
 
  for spin_i in 1:num_spin, valley in 1:num_valley, ja in eachindex(allowedq)
   PW_basis_eig[spin_i,valley,ja]=reshape(eigenvector[spin_i,valley,ja]*HF_eigenvector[spin_i,valley,ja],num_layer*num_sub*length(wave),Nband)
  end

  sub_op=zeros(ComplexF64,num_layer,num_sub,length(wave),num_layer,num_sub,length(wave))
  for ja in eachindex(wave), jb in 1:num_layer
     sub_op[jb,:,ja,jb,:,ja]=[1.0 0.0;0.0 -1.0]
  end
  
  reshape_subop=reshape(sub_op,num_layer*num_sub*length(wave),num_layer*num_sub*length(wave))

  sub_exp_HF=[zeros(ComplexF64,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)]
  sub_eig_HF=[zeros(ComplexF64,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)]
 
  for spin_i in 1:num_spin, valley in 1:num_valley, ja in eachindex(allowedq)
   sd=PW_basis_eig[spin_i,valley,ja]'*reshape_subop*PW_basis_eig[spin_i,valley,ja]
   sub_exp_HF[spin_i,valley,ja]=diag(sd,0) #This can be checked with the sub_pol
   FFF=eigen(sd)
   sub_eig_HF[spin_i,valley,ja]=FFF.values #It would be a good check that this agrees in the single-particle case.
  end
  
 chern_sub_operator=[zeros(ComplexF64,num_layer*num_sub*length(wave),num_layer*num_sub*length(wave)) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)]
 
 for spin_i in 1:num_spin, valley in 1:num_valley, ja in eachindex(allowedq)
  projected_subop=eigenvector[spin_i,valley,ja]'*reshape_subop*eigenvector[spin_i,valley,ja] #sigma z in the 6 by 6 basis
  FFF=eigen(projected_subop)
  chern_sub_eigenvectors=eigenvector[spin_i,valley,ja]*(FFF.vectors)
  for bandi in eachindex(FFF.values)
    chern_sub_operator[spin_i,valley,ja]+=sign(FFF.values[bandi])*chern_sub_eigenvectors[:,bandi]*chern_sub_eigenvectors[:,bandi]'
  end
 end

 
 chern_exp_HF=[zeros(ComplexF64,Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)]
 
 for spin_i in 1:num_spin, valley in 1:num_valley, ja in eachindex(allowedq)
  FFF=PW_basis_eig[spin_i,valley,ja]'*chern_sub_operator[spin_i,valley,ja]*PW_basis_eig[spin_i,valley,ja]
  chern_exp_HF[spin_i,valley,ja]=diag(FFF,0)
 end
 

 return sub_exp_HF,chern_exp_HF,sub_eig_HF

end



function get_chernnumber(HF_eigenvector::Array{Matrix{ComplexF64}},eigenvector::Array{Matrix{ComplexF64}},allowedq::Vector{Vector{Int}},allowedq_dic,wave::Vector{Vector{Int}},wave_dic,Nband::Int,geonum::Int,Minv::Matrix{Int})
  num_spin=2
  num_valley=2
 num_layer=3
 num_sub=2
 l1,l2,Nx,Ny=Geometry(geonum)

  PW_basis_eig=[zeros(ComplexF64,num_layer*num_sub*length(wave),Nband) for _ in 1:num_spin,_ in 1:num_valley,_ in eachindex(allowedq)]
 
  for spin_i in 1:num_spin, valley in 1:num_valley, ja in eachindex(allowedq)
   PW_basis_eig[spin_i,valley,ja]=eigenvector[spin_i,valley,ja]*HF_eigenvector[spin_i,valley,ja]
  end

  chern_num=[zeros(ComplexF64,Nband) for _ in 1:num_spin,_ in 1:num_valley]
  chern_num_nonabelian=[zeros(ComplexF64,Nband) for _ in 1:num_spin,_ in 1:num_valley]

  for spin_i in 1:num_spin, valley in 1:num_valley, bandin in 1:Nband
  
  
    chern_eigenvector=zeros(ComplexF64,num_layer*num_sub*length(wave),Nx+1,Ny+1)
    for ja in 1:Nx+1, jb in 1:Ny+1
        
        wavevec=sendtomesh(Minv,[ja-1,jb-1])
        if wavevec==[ja-1,jb-1]
          chern_eigenvector[:,ja,jb]=PW_basis_eig[spin_i,valley,allowedq_dic[wavevec]][:,bandin]
        
        else
            shuffle_vec=[ja-1,jb-1]-wavevec
            chern_eigenvector[:,ja,jb]=shuffle_vector_singlevector(PW_basis_eig[spin_i,valley,allowedq_dic[wavevec]][:,bandin],-shuffle_vec,wave,wave_dic)
        end
    end


    Uonelink=zeros(ComplexF64,Nx,Ny+1)
    Utwolink=zeros(ComplexF64,Nx+1,Ny)
 
    for ja in 1:Nx, jb in 1:Ny+1
        Uonelink[ja,jb]=dot(chern_eigenvector[:,ja,jb],chern_eigenvector[:,ja+1,jb])/abs(dot(chern_eigenvector[:,ja,jb],chern_eigenvector[:,ja+1,jb]))
    end
    
        
        
    for ja in 1:Nx+1, jb in 1:Ny
        Utwolink[ja,jb]=dot(chern_eigenvector[:,ja,jb],chern_eigenvector[:,ja,jb+1])/abs(dot(chern_eigenvector[:,ja,jb],chern_eigenvector[:,ja,jb+1]))
    end
        
    Flink=zeros(ComplexF64,Nx,Ny)
    for ja in 1:Nx, jb in 1:Ny
     Flink[ja,jb]=log(Uonelink[ja,jb]*Utwolink[ja+1,jb]/(Uonelink[ja,jb+1]*Utwolink[ja,jb]))
    end
    chern_num[spin_i,valley][bandin]=sum(Flink)/(2*π*im)
  
  
  
  end



  for spin_i in 1:num_spin, valley in 1:num_valley, bandin in 1:Nband
  
  
    chern_eigenvector=zeros(ComplexF64,num_layer*num_sub*length(wave),bandin,Nx+1,Ny+1)
    for ja in 1:Nx+1, jb in 1:Ny+1
        
        wavevec=sendtomesh(Minv,[ja-1,jb-1])
        if wavevec==[ja-1,jb-1]
          chern_eigenvector[:,:,ja,jb]=PW_basis_eig[spin_i,valley,allowedq_dic[wavevec]][:,1:bandin]
        
        else
            shuffle_vec=[ja-1,jb-1]-wavevec
            for ss in 1:bandin
            chern_eigenvector[:,ss,ja,jb]=shuffle_vector_singlevector(PW_basis_eig[spin_i,valley,allowedq_dic[wavevec]][:,ss],-shuffle_vec,wave,wave_dic)
            end
        end
    end


    Uonelink=zeros(ComplexF64,Nx,Ny+1)
    Utwolink=zeros(ComplexF64,Nx+1,Ny)
 
    for ja in 1:Nx, jb in 1:Ny+1
        Uonelink[ja,jb]=det(chern_eigenvector[:,:,ja,jb]'*chern_eigenvector[:,:,ja+1,jb])/abs(det(chern_eigenvector[:,:,ja,jb]'*chern_eigenvector[:,:,ja+1,jb]))
    end
    
        
        
    for ja in 1:Nx+1, jb in 1:Ny
        Utwolink[ja,jb]=det(chern_eigenvector[:,:,ja,jb]'*chern_eigenvector[:,:,ja,jb+1])/abs(det(chern_eigenvector[:,:,ja,jb]'*chern_eigenvector[:,:,ja,jb+1]))
    end
        
    Flink=zeros(ComplexF64,Nx,Ny)
    for ja in 1:Nx, jb in 1:Ny
     Flink[ja,jb]=log(Uonelink[ja,jb]*Utwolink[ja+1,jb]/(Uonelink[ja,jb+1]*Utwolink[ja,jb]))
    end
    chern_num_nonabelian[spin_i,valley][bandin]=sum(Flink)/(2*π*im)
  
  
  
  end




  return chern_num,chern_num_nonabelian


end






function get_formfactors(allowedq::Vector{Vector{Int}},wave::Vector{Vector{Int}},wave_diff::Vector{Vector{Int}},wave_dic::Dict{Vector{Int64},Int64},Minv::Matrix{Int},Nband::Int,eigenvector::Array{Matrix{ComplexF64}})::Array{Matrix{ComplexF64}}
 num_spin=2
 num_valley=2
 println("start FF")
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
  tic=time()
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
 println("finish FF")
 toc=time()
 println("form factors takes time",toc-tic)
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
   for jqg in eachindex(wave_diff)
    Cq=Coulomb(allowedq[jqmesh]+wave_diff[jqg],T1,T2)
    for spin_i in 1:num_spin, valley in 1:num_valley
     Fock_threaded[k2][spin_i,valley]+=Cq*form_factor[spin_i,valley,k2,jqmesh,jqg]*projector[spin_i,valley,k1pos]*(form_factor[spin_i,valley,k2,jqmesh,jqg])'
    end
   end
  end
 end 
 
 HartreeDensity=zeros(ComplexF64,length(wave_diff))
 zeropos=allowedq_dic[[0,0]]

 Threads.@threads for jqg in eachindex(wave_diff)
 for  jb in eachindex(allowedq),spin_i in 1:num_spin, valley in 1:num_valley
    HartreeDensity[jqg]+=tr(projector[spin_i,valley,jb]*(form_factor[spin_i,valley,jb,zeropos,jqg])')
 end
 end

 Threads.@threads for ja in eachindex(allowedq)
 for spin_i in 1:num_spin, valley in 1:num_valley, jqg in eachindex(wave_diff)
      Hartree_threaded[ja][spin_i,valley]+=Coulomb(wave_diff[jqg],T1,T2)*form_factor[spin_i,valley,ja,zeropos,jqg]*HartreeDensity[jqg]
 end
 end

 HF_eigenvalue=[zeros(ComplexF64,Nband) for _ in 1:num_spin,_ in 1:num_valley, _ in eachindex(allowedq)]
 HF_eigenvector=[zeros(ComplexF64,Nband,Nband) for _ in 1:num_spin,_ in 1:num_valley, _ in eachindex(allowedq)]

 Threads.@threads for ja in eachindex(allowedq)
 for spin_i in 1:num_spin, valley in 1:num_valley
    FFF=eigen(constq*Hartree_threaded[ja][spin_i,valley]+single_Ham[spin_i,valley,ja]-constq*Fock_threaded[ja][spin_i,valley])
    HF_eigenvalue[spin_i,valley,ja]=(FFF.values)
    HF_eigenvector[spin_i,valley,ja]=FFF.vectors
 end
end



 sorted=sort(reduce(vcat,reduce(vcat,real.(HF_eigenvalue)))) 
 bound=(sorted[Npa+1]+sorted[Npa])/2

 Threads.@threads for ja in eachindex(allowedq)
  for spin_i in 1:num_spin, valley in 1:num_valley
       for jd in eachindex(HF_eigenvalue[spin_i,valley,ja])
          if real(HF_eigenvalue[spin_i,valley,ja][jd])<bound
             New_projector[spin_i,valley,ja]+=HF_eigenvector[spin_i,valley,ja][:,jd]*(HF_eigenvector[spin_i,valley,ja][:,jd])'
          end
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




function plot_Chargedensity(eigenvector,total_projector,a1m,a2m,wave,T1,T2,allowedq)
  N3=50
  num_spin=2
  num_valley=2
  num_sub=2
  num_layer=3
  dimension=length(wave)*num_sub*num_layer

  zgrid=zeros(Float64,N3,N3,num_spin,num_valley,num_layer,num_sub)

  xgrid=zeros(Float64,N3,N3)
  ygrid=zeros(Float64,N3,N3)
  
  H_Density=[zeros(ComplexF64,dimension,dimension) for _ in 1:num_spin,_ in 1:num_valley]
  for k_index in eachindex(allowedq), spin_i in 1:num_spin, valley in 1:num_valley
     H_Density[spin_i,valley]+=eigenvector[spin_i,valley,k_index]*total_projector[spin_i,valley,k_index]*eigenvector[spin_i,valley,k_index]'
  end

  reshaped_Hartree_Density=[zeros(ComplexF64,num_layer,num_sub,length(wave),num_layer,num_sub,length(wave)) for _ in 1:num_spin,_ in 1:num_valley]
  
  for spin_i in 1:num_spin, valley in 1:num_valley
    reshaped_Hartree_Density[spin_i,valley]=reshape(H_Density[spin_i,valley],num_layer,num_sub,length(wave),num_layer,num_sub,length(wave))
  end
  #=
  zgrid_threaded=[zeros(Float64,N3,num_spin,num_valley,num_layer,num_sub) for _ in 1:N3]

 Threads.@threads for ja in 1:50
  for jb in 1:50, spin_i in 1:num_spin, valley in 1:num_valley, sub_index in 1:num_sub, layer_index in 1:num_layer
    rvec=ja/50*a1m+jb/50*a2m
    for jc in eachindex(wave), jd in eachindex(wave)
      gvec=[T1 T2]*(wave[jc]-wave[jd])
     zgrid[ja,jb,spin_i,valley,layer_index,sub_index]+=real(reshaped_Hartree_Density[spin_i,valley][layer_index,sub_index,jc,layer_index,sub_index,jd]*exp(im*(gvec[1]*rvec[1]+gvec[2]*rvec[2])))
    end
  
  end
 end
 =#

 Threads.@threads for ja in 1:50
  for jb in 1:50
    rvec=ja/50*a1m+jb/50*a2m
    v1=[exp(-im*dot([T1 T2]*wave[jv],rvec)) for jv in eachindex(wave)]
    for  spin_i in 1:num_spin, valley in 1:num_valley, sub_index in 1:num_sub, layer_index in 1:num_layer
    
     zgrid[ja,jb,spin_i,valley,layer_index,sub_index]+=real(v1'*reshaped_Hartree_Density[spin_i,valley][layer_index,sub_index,:,layer_index,sub_index,:]*v1)

    end
  
  end
 end

#=
 for ja in 1:50
  zgrid[ja,:,:,:,:,:]=zgrid_threaded[ja]
 end
 =#

  for ja in 1:50, jb in 1:50
    rvec=ja/50*a1m+jb/50*a2m
    xgrid[ja,jb]=rvec[1]
    ygrid[ja,jb]=rvec[2]
  end




  return xgrid,ygrid,zgrid
end