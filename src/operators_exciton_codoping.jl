using LinearAlgebra
using Arpack
using Combinatorics

using Random

#The unit cell length is different here


 
 
 
function Coulomb(qab::Float64,k::Vector{Int64},Ld::Float64)::Float64
  
   #return k==[0,0] ? 0.0 : (1/qab)
   #return qab<0.5 ? 0.0 : (exp(-qab*Ld)/qab)
   #return exp(-qab*Ld)/(qab+1.0)
  return k==[0,0] ? 30.0 : (tanh(qab*30)/qab)
 
end


function CoulombMatrix(k::Vector{Int64},T1::Vector{Float64},T2::Vector{Float64})::Matrix{Float64}
  
  #vp=[2.73,0.73,0.0]
  #vp=[0.0 2.0 2.0;2.0 0.0 0.73;2.0 0.73 0.0]
  
  Cmatrix=Matrix{Float64}(undef,3,3)
  qab=norm(k[1]*T1+k[2]*T2)
  for L1 in 1:3, L2 in 1:3
    Cmatrix[L1,L2]=Coulomb(qab,k,0.0)
  end
  return Cmatrix
  
end


function triangle_initial_Densitymatrix_control(parameters::Vector{Float64},geonum::Int64)
 

  mt=parameters[1]
  mm=parameters[2]
  mb=parameters[3]
  Vt=parameters[4]
  ϕt=parameters[5]/180*π
  Vm=parameters[6]
  ϕm=parameters[7]/180*π
  Vb=parameters[8]
  ϕb=parameters[9]/180*π
  ϵr=parameters[10]
  Eg=parameters[11]
  θ=parameters[12]/180*π
  w=parameters[13]


  Nx,Ny,l1,l2=Geometry(geonum)




  alattice=0.352# In units of NM
  blattice=4*π/(√3*alattice)

  
  
  bm=blattice*θ
  am=4*π/(√3*bm);
  constt=38.09981949*1/(mt)
  constm=38.09981949*1/(mm)
  constb=38.09981949*1/(mb) # the 38 is \hbar^2/(2me*nm^2) in units of meV
  constq=1/(Nx*Ny*3^(1/2)/2*am^2)*1/ϵr*9047.5636

  
  b1=bm*[1,0]
  b2=bm*[-1/2,√3/2]
  

  a1m=am*[√3/2,1/2]
  a2m=am*[0,1]

  L1=l1[1]*a1m+l1[2]*a2m;
  L2=l2[1]*a1m+l2[2]*a2m;
  area=abs(L1[1]*L2[2]-L2[1]*L1[2]);
  Rotminus90=[0 1;-1 0]
  T1=2*π/area*Rotminus90*L2
  T2=-2*π/area*Rotminus90*L1

  b1T=Int.(round.(inv([T1 T2])*b1))
  b2T=Int.(round.(inv([T1 T2])*b2))
 
  
  κp=bm*[-1/2,1/(2√3)]
  κm=bm*[-1/2,-1/(2√3)]
  
  
 
  Minv=[b1T';b2T']
  
  allowedq=Vector{Int64}[]
  for jb in 0:Ny-1,ja in 0:Nx-1
      push!(allowedq,sendtomesh(Minv,[ja,jb]))
  end

  
  
  
  wave=Vector{Int64}[]
  cutoff=18
  cutoffstandard=3.51*bm
  for ja in -cutoff:cutoff, jb in -cutoff:cutoff
      gtest=ja*b1+jb*b2;
      if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
          push!(wave,ja*b1T+jb*b2T)
      end
  end
  dimension=3*length(wave)

  #shift=1/2*T1+1/2*T2
 shift=[0,0]



  single_Ham=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nx*Ny]
  single_MoirePo=zeros(ComplexF64,dimension,dimension)
  tunnel=zeros(ComplexF64,dimension,dimension)
  M_tunnel=zeros(ComplexF64,dimension,dimension)
 
  
  single_eigenvalue=[[zeros(Float64,dimension) for _ in 1:2] for _ in 1:Nx*Ny]
  single_eigenvector=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nx*Ny]
  uncoupled_eigenvector=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nx*Ny]
 
 
  
  for jc in eachindex(wave)
    pos=findfirst(item->item==wave[jc]-b1T,wave)
    if pos≠nothing
      single_MoirePo[3*(jc-1)+1:3*jc,3*(pos-1)+1:3*pos]=diagm([Vt*exp(im*ϕt),Vm*exp(im*ϕm),Vb*exp(im*ϕb)])
    end
    
  
    pos=findfirst(item->item==wave[jc]+b2T+b1T,wave)
    if pos≠nothing
      single_MoirePo[3*(jc-1)+1:3*jc,3*(pos-1)+1:3*pos]=diagm([Vt*exp(im*ϕt),Vm*exp(im*ϕm),Vb*exp(im*ϕb)])
    end
  
    pos=findfirst(item->item==wave[jc]-b2T,wave)
    if pos≠nothing
      single_MoirePo[3*(jc-1)+1:3*jc,3*(pos-1)+1:3*pos]=diagm([Vt*exp(im*ϕt),Vm*exp(im*ϕm),Vb*exp(im*ϕb)])
    end
  end
  
  
  single_MoirePo=single_MoirePo+single_MoirePo'
  
  for jc in eachindex(wave)
   
    tunnel[3*(jc-1)+3,3*(jc-1)+2]=w
  
    pos=findfirst(item->item==wave[jc]+b2T+b1T,wave)
    if pos≠nothing
   
      tunnel[3*(jc-1)+3,3*(pos-1)+2]=w
    end
  
    pos=findfirst(item->item==wave[jc]+b2T,wave)
    if pos≠nothing
     
      tunnel[3*(jc-1)+3,3*(pos-1)+2]=w
    end
  end
  
  tunnel=tunnel+tunnel'
  
  for jc in eachindex(wave)

    M_tunnel[3*(jc-1)+3,3*(jc-1)+2]=w
  
    pos=findfirst(item->item==wave[jc]-b2T-b1T,wave)
    if pos≠nothing
     
      M_tunnel[3*(jc-1)+3,3*(pos-1)+2]=w
    end
  
    pos=findfirst(item->item==wave[jc]-b2T,wave)
    if pos≠nothing
 
      M_tunnel[3*(jc-1)+3,3*(pos-1)+2]=w
    end
  end
  
  M_tunnel=M_tunnel+M_tunnel'
  
  
  Threads.@threads for ja in 1:Nx*Ny
    
    
    k=allowedq[ja][1]*T1+allowedq[ja][2]*T2
  
    for jb in eachindex(wave)
     single_Ham[ja][1][3*(jb-1)+1,3*(jb-1)+1]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2-κp+shift)^2*constt+Eg #This is the configuration
     single_Ham[ja][1][3*(jb-1)+2,3*(jb-1)+2]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2-κp+shift)^2*constm
     single_Ham[ja][1][3*(jb-1)+3,3*(jb-1)+3]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2-κm+shift)^2*constb
    end
  
    for jb in eachindex(wave)
      single_Ham[ja][2][3*(jb-1)+1,3*(jb-1)+1]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2+κp+shift)^2*constt+Eg
      single_Ham[ja][2][3*(jb-1)+2,3*(jb-1)+2]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2+κp+shift)^2*constm
      single_Ham[ja][2][3*(jb-1)+3,3*(jb-1)+3]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2+κm+shift)^2*constb
     end
    
  
     
    FFF=eigen(single_Ham[ja][1])
     uncoupled_eigenvector[ja][1]=FFF.vectors

    FFF=eigen(single_Ham[ja][2])
     uncoupled_eigenvector[ja][2]=FFF.vectors

  
  
    single_Ham[ja][1]+=single_MoirePo+tunnel
    single_Ham[ja][2]+=single_MoirePo+M_tunnel
   
    FFF=eigen(single_Ham[ja][1])
    single_eigenvalue[ja][1]=real(FFF.values)
    single_eigenvector[ja][1]=FFF.vectors
  
    FFF=eigen(single_Ham[ja][2])
    single_eigenvalue[ja][2]=real(FFF.values)
    single_eigenvector[ja][2]=FFF.vectors

  end
  


  println("values",single_eigenvalue[1][1][2*length(wave)-2:2*length(wave)+2])

   
  BG_DensityMatrix=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nx*Ny]
  for ja in 1:Nx*Ny, vi in 1:2, jb in 1:length(wave)*2
    BG_DensityMatrix[ja][vi]+=(uncoupled_eigenvector[ja][vi][:,jb]*(uncoupled_eigenvector[ja][vi][:,jb])') 
    
  end
  
  input_DensityMatrix=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nx*Ny]
  
  for ja in 1:Nx*Ny, vi in 1:1, statecount in 1:2
    A=10^(-1)*randn(ComplexF64,dimension)
    A=A/norm(A)
      input_DensityMatrix[ja][vi]+=A*A'
  end
   
  for ja in 1:Nx*Ny, vi in 2:2, statecount in 1:1
    A=10^(-1)*randn(ComplexF64,dimension)
    A=A/norm(A)
      input_DensityMatrix[ja][vi]+=A*A'
  end

   #=
  for ja in 1:Nx*Ny, vi in 2:2
    A=10^(-1)*randn(ComplexF64,dimension,dimension)
      input_DensityMatrix[ja][vi]+=A+A'
  end
  =#



  
  for ja in 1:Nx*Ny, vi in 1:2, jb in 1:length(wave)*2-2
   input_DensityMatrix[ja][vi]+=(single_eigenvector[ja][vi][:,jb]*(single_eigenvector[ja][vi][:,jb])') 
  end

 
 
    input_DensityMatrix-=BG_DensityMatrix

  
 
 single_chern=zeros(ComplexF64,4)
 single_chern[1]=calculate_chern(single_eigenvector,Nx,Ny,Minv,dimension,wave,allowedq,1,2*length(wave)-1)
 single_chern[2]=calculate_chern(single_eigenvector,Nx,Ny,Minv,dimension,wave,allowedq,1,2*length(wave))
 single_chern[3]=calculate_chern(single_eigenvector,Nx,Ny,Minv,dimension,wave,allowedq,2,2*length(wave)-1)
 single_chern[4]=calculate_chern(single_eigenvector,Nx,Ny,Minv,dimension,wave,allowedq,2,2*length(wave))
 println("single_chern",single_chern)
  
 return  wave, input_DensityMatrix, BG_DensityMatrix, single_Ham, single_eigenvalue,allowedq, T1, T2, a1m, a2m, Minv,constq,single_chern
    
end











function construct_loop_dic(wave::Vector{Vector{Int}})::Dict{Vector{Int},Any}
    g_dic=Dict{Vector{Int},Int}()
    for ja in eachindex(wave)
      g_dic[wave[ja]]=ja
    end

    loop_dic=Dict{Vector{Int},Any}()


    for gindex in eachindex(wave), g2index in eachindex(wave)
        if !haskey(loop_dic,wave[g2index]-wave[gindex])
           loop_dic[wave[g2index]-wave[gindex]]=Dict{Vector{Int},Vector{Vector{Int}}}()
        end
        
        loop_dic[wave[g2index]-wave[gindex]][[gindex,g2index,3*gindex-2,3*gindex,3*g2index-2,3*g2index]]=Vector{Int64}[]
        
         
    end
    
    for gindex in eachindex(wave), g1index in eachindex(wave),g2index in eachindex(wave)
        g3=wave[gindex]+wave[g1index]-wave[g2index]
        if haskey(g_dic,wave[gindex]+wave[g1index]-wave[g2index])
            push!(loop_dic[wave[g2index]-wave[gindex]][[gindex,g2index,3*gindex-2,3*gindex,3*g2index-2,3*g2index]],[g1index,g_dic[g3],3*g1index-2,3*g1index,3*g_dic[g3]-2,3*g_dic[g3]])
        end
    end
    
   return loop_dic

end






function Construct_DensityMatrix(loop_dic::Dict{Vector{Int},Any},allowedq::Vector{Vector{Int}},T1::Vector{Float64},T2::Vector{Float64},Nx::Int64,Ny::Int64,Minv::Matrix{Int64},wave::Vector{Vector{Int64}},input_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},BG_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},single_Ham::Vector{Vector{Matrix{ComplexF64}}},constq::Float64,Parnum::Int64)::Tuple{Float64,Vector{Vector{Matrix{ComplexF64}}},Vector{Vector{Matrix{ComplexF64}}},Vector{Vector{Vector{Float64}}},Vector{Vector{Matrix{ComplexF64}}},Float64,Float64}
  
 
    dimension=3*length(wave)
    HartreeMatrix=zeros(ComplexF64,dimension,dimension)
    FockMatrix=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nx*Ny]
  
    output_DensityMatrix=[Vector{Matrix{ComplexF64}}(undef,2) for _ in 1:Nx*Ny]
    DeltaMatrix=[Vector{Matrix{ComplexF64}}(undef,2) for _ in 1:Nx*Ny]
    NewDensityMatrix=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nx*Ny]
    HF_eigenvalue=[Vector{Vector{Float64}}(undef,2) for _ in 1:Nx*Ny]
    HF_eigenvector=[Vector{Matrix{ComplexF64}}(undef,2) for _ in 1:Nx*Ny]
  
   tic=time()
    Threads.@threads for jk in 1:Nx*Ny
      Fk = FockMatrix[jk]
      for jk1 in 1:Nx*Ny
          dmk = input_DensityMatrix[jk1]
          q=allowedq[jk1]-allowedq[jk]
         for (dg,loop_dic_dg) in loop_dic
             CoulF=CoulombMatrix(q+dg,T1,T2)     
             for (gg2,loop_dic_dg_gg2) in loop_dic_dg         
             for g1g3 in loop_dic_dg_gg2
                 Fk[1][g1g3[5]:g1g3[6],gg2[3]:gg2[4]]+=dmk[1][g1g3[3]:g1g3[4],gg2[5]:gg2[6]] .* CoulF
                 Fk[2][g1g3[5]:g1g3[6],gg2[3]:gg2[4]]+=dmk[2][g1g3[3]:g1g3[4],gg2[5]:gg2[6]] .* CoulF   
              end 
             end
     
         end    
      end
    end
    toc=time()
    println(toc-tic)
  
    
    Hartree_Density=zeros(ComplexF64,dimension,dimension)
    
    for ja in 1:Nx*Ny
     Hartree_Density+=input_DensityMatrix[ja][1]+input_DensityMatrix[ja][2]
    end
  
   
     for dg in keys(loop_dic)
            CoulH=CoulombMatrix(dg,T1,T2)
            for gg2 in keys(loop_dic[dg])                   
            for g1g3 in loop_dic[dg][gg2]           
                HartreeMatrix[gg2[5]:gg2[6],gg2[3]:gg2[4]]+=diagm(CoulH*diag(Hartree_Density[g1g3[3]:g1g3[4],g1g3[5]:g1g3[6]],0))            
            end 
            end
      end    
   
    
  
   
  
   for ja in 1:Nx*Ny, jb in 1:2
     FFF=eigen(single_Ham[ja][jb]+constq*HartreeMatrix-constq*FockMatrix[ja][jb])
     HF_eigenvalue[ja][jb]=real(FFF.values)
     HF_eigenvector[ja][jb]=FFF.vectors
   end


   
    #=
   htelement=sum(abs.(vec(HartreeMatrix)))
    println("htelement",htelement)



   fkelement=0
    for ja in 1:Nx*Ny, jb in 1:2
     fkelement+=sum(abs.(vec(FockMatrix[ja][jb])))/(Nx*Ny)
    end
  
  println("fkelement",fkelement)
   =#
   bound=(sort(reduce(vcat,reduce(vcat,HF_eigenvalue)))[Parnum+1]+sort(reduce(vcat,reduce(vcat,HF_eigenvalue)))[Parnum])/2
   
   
    for ja in 1:Nx*Ny,vi in 1:2   
         for jd in eachindex(HF_eigenvalue[ja][vi])
            if HF_eigenvalue[ja][vi][jd]<bound
               NewDensityMatrix[ja][vi]+=HF_eigenvector[ja][vi][:,jd]*(HF_eigenvector[ja][vi][:,jd])'
            end
         end
         NewDensityMatrix[ja][vi]-=BG_DensityMatrix[ja][vi]
    end
  


   
    DeltaMatrix=NewDensityMatrix-input_DensityMatrix
  
    output_DensityMatrix=0.0*input_DensityMatrix+1.0*NewDensityMatrix

    l3=0.0
    for ja in 1:length(wave), jb in 1:Nx*Ny,vi in 1:2
     l3+=NewDensityMatrix[jb][vi][3*(ja-1)+1,3*(ja-1)+1]
    end
   println("l3=",l3)
    
    e1=0.0
    
    #e2=0.0
    for ja in 1:Nx*Ny,vi in 1:2
      e1+=tr(DeltaMatrix[ja][vi]'*DeltaMatrix[ja][vi])
      #e2+=sum((abs.(diag(DeltaMatrix[ja][vi],0))).^2)
    end
    eout=real(e1)/(2*Nx*Ny)
    #println("e2=",sqrt(e2/(12*Nx*Ny*length(wave)^2)))
   
    energy=0.0
    for ja in 1:Nx*Ny, vi in 1:2
      Energy_Matrix=single_Ham[ja][vi]+(constq*HartreeMatrix-constq*FockMatrix[ja][vi])/2
        energy+=tr(Energy_Matrix*input_DensityMatrix[ja][vi])
    end
    println("energy=",energy)

    HF_chern=zeros(ComplexF64,4)
    dimension=3*length(wave)
    HF_chern[1]=calculate_chern(HF_eigenvector,Nx,Ny,Minv,dimension,wave,allowedq,1,2*length(wave)-1)
    HF_chern[2]=calculate_chern(HF_eigenvector,Nx,Ny,Minv,dimension,wave,allowedq,1,2*length(wave))
    HF_chern[3]=calculate_chern(HF_eigenvector,Nx,Ny,Minv,dimension,wave,allowedq,2,2*length(wave)-1)
    HF_chern[4]=calculate_chern(HF_eigenvector,Nx,Ny,Minv,dimension,wave,allowedq,2,2*length(wave))

    println("HFchern",HF_chern)
    


  
   return  eout,output_DensityMatrix,DeltaMatrix,HF_eigenvalue,HF_eigenvector,bound,real(energy)
end
  
  






function iteration_loop(initial_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},BG_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},allowedq::Vector{Vector{Int64}},T1::Vector{Float64},T2::Vector{Float64},geonum::Int64,Minv::Matrix{Int64},wave::Vector{Vector{Int64}},single_Ham::Vector{Vector{Matrix{ComplexF64}}},constq::Float64,holenum::Int64)::Tuple{Vector{Vector{Vector{Matrix{ComplexF64}}}},Vector{Vector{Vector{Matrix{ComplexF64}}}},Vector{Vector{Vector{Float64}}},Vector{Vector{Matrix{ComplexF64}}},Float64,Float64,Vector{ComplexF64},Vector{Vector{Matrix{ComplexF64}}}}
    
  Nx,Ny,l1,l2=Geometry(geonum)
    eout=1.0
    itcount=0
    loop_dic=construct_loop_dic(wave)
    HF_eigenvalue=[Vector{Vector{Float64}}(undef,2) for _ in 1:Nx*Ny]
    HF_eigenvector=[Vector{Matrix{ComplexF64}}(undef,2) for _ in 1:Nx*Ny]
    DIIS_input_DensityMatrix=Vector{Vector{Vector{Matrix{ComplexF64}}}}(undef,3)
    DIIS_input_DeltaMatrix=Vector{Vector{Vector{Matrix{ComplexF64}}}}(undef,3)
    input_DensityMatrix=initial_DensityMatrix
    bad_count=0
    Parnum=2*Nx*Ny*(2*length(wave))-holenum*Nx*Ny
    bound=0.0
    energy=0.0

    while (eout>1*10^(-13)) || (bad_count<3)
      if  eout<1*10^(-13)
        bad_count+=1
      end
      tic=time()
      eout,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue, HF_eigenvector,bound,energy=Construct_DensityMatrix(loop_dic,allowedq,T1,T2,Nx,Ny,Minv,wave,input_DensityMatrix,BG_DensityMatrix,single_Ham,constq,Parnum)
                                                                                                                                          
      DIIS_input_DensityMatrix[mod(itcount,3)+1]=input_DensityMatrix
      input_DensityMatrix=output_DensityMatrix
      itcount+=1
      toc=time()
      println(toc-tic,"eout=$eout")
      flush(stdout)
     
    
    end


#=
    println("startDIIS",itcount)
    
    bad_count=0
    while (eout>10^-13) || (bad_count<4)
        if  eout<1*10^-13 
            bad_count+=1
        end
        tic=time()
        Bmatrix=zeros(ComplexF64,4,4)
        for ja in 1:3
         Bmatrix[ja,4]=1
         Bmatrix[4,ja]=1
        end
    
        for ja in 1:3,jb in 1:3
            for jc in 1:Nx*Ny, vi in 1:2
               Bmatrix[ja,jb]+=tr((DIIS_input_DeltaMatrix[ja][jc][vi])'*(DIIS_input_DeltaMatrix[jb][jc][vi]))
            end
        end
        coeff=inv(Bmatrix)*[0;0;0;1]
       
        dmk=coeff[1]*(DIIS_input_DensityMatrix[1]+DIIS_input_DeltaMatrix[1])+coeff[2]*(DIIS_input_DensityMatrix[2]+DIIS_input_DeltaMatrix[2])+coeff[3]*(DIIS_input_DensityMatrix[3]+DIIS_input_DeltaMatrix[3])
        eout,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue, HF_eigenvector,bound,energy=Construct_DensityMatrix(loop_dic,allowedq,T1,T2,Nx,Ny,Minv,wave,input_DensityMatrix,BG_DensityMatrix,single_Ham,constq,Parnum)
        DIIS_input_DensityMatrix[mod(itcount,3)+1]=dmk
        itcount+=1

        toc=time()
        println(toc-tic,"eout=$eout")
        flush(stdout)
    end

    =# 






    
    HF_chern=zeros(ComplexF64,4)
    dimension=3*length(wave)
    HF_chern[1]=calculate_chern(HF_eigenvector,Nx,Ny,Minv,dimension,wave,allowedq,1,2*length(wave)-1)
    HF_chern[2]=calculate_chern(HF_eigenvector,Nx,Ny,Minv,dimension,wave,allowedq,1,2*length(wave))
    HF_chern[3]=calculate_chern(HF_eigenvector,Nx,Ny,Minv,dimension,wave,allowedq,2,2*length(wave)-1)
    HF_chern[4]=calculate_chern(HF_eigenvector,Nx,Ny,Minv,dimension,wave,allowedq,2,2*length(wave))


     
    dope_hole_DM=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nx*Ny]
    for ja in 1:Nx*Ny, vi in 1:2, jb in length(wave)*2:length(wave)*2
      dope_hole_DM[ja][vi]+=(HF_eigenvector[ja][vi][:,jb]*(HF_eigenvector[ja][vi][:,jb])') 
      
    end




  return DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue, HF_eigenvector,bound,energy,HF_chern,dope_hole_DM

end



function calculate_layerpolarization(wave,HF_eigenvector,layer_num,Nx,Ny)
  layer_pol=zeros(Float64,layer_num,2,Nx*Ny,layer_num*length(wave))

    for ja in 1:layer_num, vi in 1:2, jc in 1:Nx*Ny, jd in 1:layer_num*length(wave),je in 1:length(wave)
      layer_pol[ja,vi,jc,jd]+=abs(HF_eigenvector[jc][vi][layer_num*(je-1)+ja,jd])^2
    end
    return layer_pol
end

function get_holeband(wave,HF_eigenvector,Nx,Ny)
   vec_1=Vector{Any}(undef,Nx*Ny)
   vec_2=Vector{Any}(undef,Nx*Ny)

   for ja in 1:Nx*Ny
   vec_1[ja]=HF_eigenvector[ja][1][:,2*length(wave)]
   vec_2[ja]=HF_eigenvector[ja][2][:,2*length(wave)]
  end
  return vec_1, vec_2
end




function Densitymap(a1m::Vector{Float64},a2m::Vector{Float64},wave::Vector{Vector{Int}},input_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},T1::Vector{Float64},T2::Vector{Float64})::Tuple{Array{ComplexF64},Array{ComplexF64},Array{ComplexF64}}
  dimension=3*length(wave)  
     N3=50
    xgrid=zeros(Float64,N3,N3)
    ygrid=zeros(Float64,N3,N3)
    zgrid=zeros(ComplexF64,N3,N3,2,3)
    Total_Density=[zeros(ComplexF64,dimension,dimension) for _ in 1:2]
    for jk1 in eachindex(input_DensityMatrix)
       Total_Density+=input_DensityMatrix[jk1]
    end
    
    for ja in 1:N3, jb in 1:N3
      xgrid[ja,jb]=(ja/N3*a1m[1]+jb/N3*a2m[1])
      ygrid[ja,jb]=(ja/N3*a1m[2]+jb/N3*a2m[2])
      rvec=ja/N3*a1m+jb/N3*a2m
      for vi in 1:2, Li in 1:3, jc in eachindex(wave), jd in eachindex(wave)
        gvec=[T1 T2]*(wave[jc]-wave[jd])
       zgrid[ja,jb,vi,Li]+=Total_Density[vi][3*(jc-1)+Li,3*(jd-1)+Li]*exp(im*(gvec[1]*rvec[1]+gvec[2]*rvec[2]))
      end
    
    end

    return xgrid,ygrid,zgrid
end

function sendtomesh(Minv::Matrix{Int64},q1::Vector{Int})::Vector{Int}
  Qvec=q1'*inv(Minv)
  return q1 .-vec((Int.(floor.(round.(Qvec,digits=5)))*Minv)')
end


function calculate_chern(eigenvector_matrix,Nx,Ny,Minv,dimension,wave,allowedq,vi,bi)
  chern_eigenvector=zeros(ComplexF64,dimension,Nx+1,Ny+1)
  for ja in 1:Nx+1, jb in 1:Ny+1
    if (ja<Nx+1) && (jb<Ny+1)
      pos=findfirst(item->item==[ja-1,jb-1],allowedq)
      chern_eigenvector[:,ja,jb]=eigenvector_matrix[pos][vi][:,bi]
    else
      G1=[ja-1,jb-1]-sendtomesh(Minv,[ja-1,jb-1])
      k=sendtomesh(Minv,[ja-1,jb-1])
      pos=findfirst(item->item==k,allowedq)
      for jc in eachindex(wave)
        pos_2=findfirst(item->item==G1+wave[jc],wave)
        if pos_2≠nothing
          chern_eigenvector[3*(jc-1)+1:3*jc,ja,jb]=eigenvector_matrix[pos][vi][3*(pos_2-1)+1:3*pos_2,bi]
        end
      end
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
  chern=sum(Flink)/(2*π*im)



  return chern
end








function Geometry(geonum::Int64)
  if geonum==1
      Nx=3;
      Ny=3;
      l1=[3,0]
      l2=[0,3]
  end
  
  if geonum==2
   Nx=3;
   Ny=4;
   l1=[3,0]
   l2=[0,4]
  end

  if geonum==3
   Nx=3;
   Ny=5;
   l1=[3,0]
   l2=[0,5]
  end

  if geonum==4
   Nx=4;
   Ny=4;
   l1=[4,0]
   l2=[0,4]
  end

  if geonum==5
   Nx=3;
   Ny=6;
   l1=[3,0]
   l2=[0,6]
  end

  if geonum==6
   Nx=4;
   Ny=5;
   l1=[4,0]
   l2=[0,5]
  end

  if geonum==7
   Nx=4;
   Ny=6;
   l1=[4,0]
   l2=[0,6]
  end
   


  if geonum==8
   Nx=4
   Ny=7;
   l1=[4,0]
   l2=[0,7]
  end

  
  if geonum==9
      Nx=5
      Ny=6
      l1=[5,0]
      l2=[0,6]
  end

  if geonum==10
      Nx=6
      Ny=6
      l1=[6,0]
      l2=[0,6]
  end

  return Nx,Ny,l1,l2
end