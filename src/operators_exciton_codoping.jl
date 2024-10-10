using LinearAlgebra
using Arpack
using Combinatorics

using Random

#The unit cell length is different here


 
 
 
function Coulomb(qab::Float64,k::Vector{Int64},Ld::Float64)::Float64
  
   #return k==[0,0] ? 0.0 : (exp(-qab*Ld)/qab)
   #return qab<0.5 ? 0.0 : (exp(-qab*Ld)/qab)
   return exp(-qab*Ld)/(qab+1.0)
end


function CoulombMatrix(k::Vector{Int64},T1::Vector{Float64},T2::Vector{Float64})::Matrix{Float64}
  
  #vp=[2.73,0.73,0.0]
  #vp=[0.0 2.0 2.0;2.0 0.0 0.73;2.0 0.73 0.0]
  vp=[0.0 2.0 0.0;2.0 0.0 2.0;0.0 2.0 0.0]
  Cmatrix=Matrix{Float64}(undef,3,3)
  qab=norm(k[1]*T1+k[2]*T2)
  for L1 in 1:3, L2 in 1:3
    #Cmatrix[L1,L2]=Coulomb(qab,k,abs(vp[L1]-vp[L2]))
    Cmatrix[L1,L2]=Coulomb(qab,k,vp[L1,L2])
  end
  return Cmatrix
  
end


function triangle_initial_Densitymatrix_control(parameters::Vector{Float64},Nq::Int64)
 

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







  alattice=0.352# In units of NM
  blattice=4*π/(√3*alattice)

  
  
  bm=blattice*θ
  am=4*π/(√3*bm);
  constt=38.09981949*1/(mt)
  constm=38.09981949*1/(mm)
  constb=38.09981949*1/(mb) # the 38 is \hbar^2/(2me*nm^2) in units of meV
  constq=1/(Nq^2*3^(1/2)/2*am^2)*1/ϵr*9047.5636

  
  b1=bm*[1,0]
  b2=bm*[-1/2,√3/2]
  

  a1m=am*[√3/2,1/2]
  a2m=am*[0,1]

  
  T1=b1/(Nq)
  T2=b2/(Nq)
  
  κp=bm*[-1/2,1/(2√3)]
  κm=bm*[-1/2,-1/(2√3)]
  
  
  b1T=Int.(round.(inv([T1 T2])*b1))
  b2T=Int.(round.(inv([T1 T2])*b2))
  
  
  
  allowedq=Vector{Int64}[]
  for ja in 0:Nq-1,jb in 0:Nq-1
      push!(allowedq,[ja,jb])
  end

  
  
  
  wave=Vector{Int64}[]
  cutoff=18
  cutoffstandard=3.01*bm
  for ja in -cutoff:cutoff, jb in -cutoff:cutoff
      gtest=ja*b1+jb*b2;
      if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
          push!(wave,ja*b1T+jb*b2T)
      end
  end
  dimension=3*length(wave)
 
  shift=1/2*T1+1/2*T2




  single_Ham=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
  single_MoirePo=zeros(ComplexF64,dimension,dimension)
  tunnel=zeros(ComplexF64,dimension,dimension)
  M_tunnel=zeros(ComplexF64,dimension,dimension)
  #seed_tunnel=zeros(ComplexF64,dimension,dimension)
  
  single_eigenvalue=[[zeros(Float64,dimension) for _ in 1:2] for _ in 1:Nq^2]
  single_eigenvector=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
  uncoupled_eigenvector=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
  #seed_eigenvector=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
 
  
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
#=
  for jc in eachindex(wave)
   
    seed_tunnel[3*(jc-1)+1,3*(jc-1)+2]+=2.0
    seed_tunnel[3*(jc-1)+1,3*(jc-1)+3]+=2.0
  end

  for jc in eachindex(wave)
    pos=findfirst(item->item==wave[jc]-b1T,wave)
    if pos≠nothing
      seed_tunnel[3*(jc-1)+1:3*jc,3*(pos-1)+1:3*pos]+=diagm([1*exp(im*ϕt),1*exp(im*ϕm),1*exp(im*ϕb)])
    end
    
  
    pos=findfirst(item->item==wave[jc]+b2T+b1T,wave)
    if pos≠nothing
      seed_tunnel[3*(jc-1)+1:3*jc,3*(pos-1)+1:3*pos]+=diagm([1*exp(im*ϕt),1*exp(im*ϕm),1*exp(im*ϕb)])
    end
  
    pos=findfirst(item->item==wave[jc]-b2T,wave)
    if pos≠nothing
      seed_tunnel[3*(jc-1)+1:3*jc,3*(pos-1)+1:3*pos]+=diagm([1*exp(im*ϕt),1*exp(im*ϕm),1*exp(im*ϕb)])
    end
  end
  
  seed_tunnel=seed_tunnel+seed_tunnel'
  =#
  
  
  
  Threads.@threads for ja in 1:Nq^2
    
    
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
    #=
     FFF=eigen(single_Ham[ja][1]+seed_tunnel)
     seed_eigenvector[ja][1]=FFF.vectors

    FFF=eigen(single_Ham[ja][2]+seed_tunnel)
     seed_eigenvector[ja][2]=FFF.vectors
    =#
  
  
    single_Ham[ja][1]+=single_MoirePo+tunnel
    single_Ham[ja][2]+=single_MoirePo+M_tunnel
   
    FFF=eigen(single_Ham[ja][1])
    single_eigenvalue[ja][1]=real(FFF.values)
    single_eigenvector[ja][1]=FFF.vectors
  
    FFF=eigen(single_Ham[ja][2])
    single_eigenvalue[ja][2]=real(FFF.values)
    single_eigenvector[ja][2]=FFF.vectors

  end
  




   
  BG_DensityMatrix=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
  for ja in 1:Nq^2, vi in 1:2, jb in 1:length(wave)*2
    BG_DensityMatrix[ja][vi]+=(uncoupled_eigenvector[ja][vi][:,jb]*(uncoupled_eigenvector[ja][vi][:,jb])') 
    
  end
  
  input_DensityMatrix=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
  
  for ja in 1:Nq^2, vi in 1:1
    A=10^(-1)*randn(ComplexF64,dimension,dimension)
      input_DensityMatrix[ja][vi]=A+A'
  end

  for ja in 1:Nq^2, vi in 2:2
    A=10^(-1)*randn(ComplexF64,dimension,dimension)
      input_DensityMatrix[ja][vi]=A+A'
  end
  



#=
  for ja in 1:Nq^2, vi in 1:1, jb in 1:length(wave)*2
   input_DensityMatrix[ja][vi]+=(seed_eigenvector[ja][vi][:,jb]*(seed_eigenvector[ja][vi][:,jb])') 
  end
  for ja in 1:Nq^2, vi in 2:2, jb in 1:length(wave)*2-2
    input_DensityMatrix[ja][vi]+=(seed_eigenvector[ja][vi][:,jb]*(seed_eigenvector[ja][vi][:,jb])') 
  end
  for ja in 1:Nq^2, vi in 2:2, jb in length(wave)*2+1:length(wave)*2+1
    input_DensityMatrix[ja][vi]+=(seed_eigenvector[ja][vi][:,jb]*(seed_eigenvector[ja][vi][:,jb])') 
  end

  input_DensityMatrix-=BG_DensityMatrix
 =#
 single_chern=zeros(ComplexF64,4)
 single_chern[1]=calculate_chern(single_eigenvector,Nq,dimension,wave,allowedq,1,2*length(wave)-1)
 single_chern[2]=calculate_chern(single_eigenvector,Nq,dimension,wave,allowedq,1,2*length(wave))
 single_chern[3]=calculate_chern(single_eigenvector,Nq,dimension,wave,allowedq,2,2*length(wave)-1)
 single_chern[4]=calculate_chern(single_eigenvector,Nq,dimension,wave,allowedq,2,2*length(wave))
 println("single_chern",single_chern)
  
 return  wave, input_DensityMatrix, BG_DensityMatrix, single_Ham, single_eigenvalue,allowedq, T1, T2, a1m, a2m,constq,single_chern
    
end



#=
function plot_band(parameters::Vector{Float64},Nq::Int64,final_densitymatrix::Vector{Vector{Matrix{ComplexF64}}})
 

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







  alattice=0.352# In units of NM
  blattice=4*π/(√3*alattice)

  
  
  bm=2*blattice*sin(θ/2)
  am=4*π/(√3*bm);
  constt=38.09981949*1/(mt)
  constm=38.09981949*1/(mm)
  constb=38.09981949*1/(mb) # the 38 is \hbar^2/(2me*nm^2) in units of meV
  constq=1/(Nq^2*3^(1/2)/2*am^2)*1/ϵr*9047.5636
  
  
  b1=bm*[1,0]
  b2=bm*[-1/2,√3/2]
  

  a1m=am*[√3/2,1/2]
  a2m=am*[0,1]
  
  
  T1=b1/(Nq)
  T2=b2/(Nq)
  
  κp=bm*[-1/2,1/(2√3)]
  κm=bm*[-1/2,-1/(2√3)]
  γ=[0.0,0.0]


  Npath=90
  path=Vector{Vector{Float64}}(undef,Npath)
  for ja in 1:Int(Npath/3)
  path[ja]=ja/Int(Npath/3)*κp+(1-ja/Int(Npath/3))*γ
  path[ja+Int(Npath/3)]=ja/Int(Npath/3)*κm+(1-ja/Int(Npath/3))*κp
  path[ja+Int(Npath/3)*2]=ja/Int(Npath/3)*γ+(1-ja/Int(Npath/3))*κm
  end

  
  b1T=Int.(round.(inv([T1 T2])*b1))
  b2T=Int.(round.(inv([T1 T2])*b2))
  
  
  
  allowedq=Vector{Int64}[]
  for ja in 0:Nq-1,jb in 0:Nq-1
      push!(allowedq,[ja,jb])
  end

  
  
  
  wave=Vector{Int64}[]
  cutoff=18
  cutoffstandard=4.01*bm
  for ja in -cutoff:cutoff, jb in -cutoff:cutoff
      gtest=ja*b1+jb*b2;
      if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
          push!(wave,ja*b1T+jb*b2T)
      end
  end

  loop_dic=construct_loop_dic(wave)



  dimension=3*length(wave)
 
  
 
  single_MoirePo=zeros(ComplexF64,dimension,dimension)
  tunnel=zeros(ComplexF64,dimension,dimension)
  M_tunnel=zeros(ComplexF64,dimension,dimension)
  
  
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
  path_eigenvalue=zeros(ComplexF64,Npath)
  M_path_eigenvalue=zeros(ComplexF64,Npath)
  
  
  
  Threads.@threads for ja in 1:Npath
    
    
    k=path[ja]
     single_Ham=zeros(CompelxF64,dimension,dimension)
    for jb in eachindex(wave)
     single_Ham[3*(jb-1)+1,3*(jb-1)+1]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2)^2*constt+Eg
     single_Ham[3*(jb-1)+2,3*(jb-1)+2]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2-κp)^2*constm
     single_Ham[3*(jb-1)+3,3*(jb-1)+3]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2-κm)^2*constb
    end
  
  

    hfmatrix=Construct_HFmatrix(loop_dic,k,1,allowedq,T1,T2,Nq,wave,final_densitymatrix,constq)
  
   
    FFF=eigen(single_Ham+tunnel+single_MoirePo+hfmatrix)
    path_eigenvalue[ja]=real(FFF.values)
    
  end


   
  Threads.@threads for ja in 1:Npath
    
    
    k=path[ja]
    single_Ham=zeros(CompelxF64,dimension,dimension)
    for jb in eachindex(wave)
     single_Ham[3*(jb-1)+1,3*(jb-1)+1]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2)^2*constt+Eg
     single_Ham[3*(jb-1)+2,3*(jb-1)+2]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2+κp)^2*constm
     single_Ham[3*(jb-1)+3,3*(jb-1)+3]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2+κm)^2*constb
    end
  
  
    hfmatrix=Construct_HFmatrix(loop_dic,k,2,allowedq,T1,T2,Nq,wave,final_densitymatrix,constq)
    
   
    FFF=eigen(single_Ham+M_tunnel+single_MoirePo+hfmatrix)
    M_path_eigenvalue[ja]=real(FFF.values)
    
  end
  




   
  

 
 

  
 return  path_eigenvalue,M_path_eigenvalue
    
end
=#







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






function Construct_DensityMatrix(loop_dic::Dict{Vector{Int},Any},allowedq::Vector{Vector{Int}},T1::Vector{Float64},T2::Vector{Float64},Nq::Int64,wave::Vector{Vector{Int64}},input_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},BG_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},single_Ham::Vector{Vector{Matrix{ComplexF64}}},constq::Float64,Parnum::Int64)::Tuple{Float64,Vector{Vector{Matrix{ComplexF64}}},Vector{Vector{Matrix{ComplexF64}}},Vector{Vector{Vector{Float64}}},Vector{Vector{Matrix{ComplexF64}}},Float64,Float64}
  
 
    dimension=3*length(wave)
    HartreeMatrix=zeros(ComplexF64,dimension,dimension)
    FockMatrix=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
  
    output_DensityMatrix=[Vector{Matrix{ComplexF64}}(undef,2) for _ in 1:Nq^2]
    DeltaMatrix=[Vector{Matrix{ComplexF64}}(undef,2) for _ in 1:Nq^2]
    NewDensityMatrix=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
    HF_eigenvalue=[Vector{Vector{Float64}}(undef,2) for _ in 1:Nq^2]
    HF_eigenvector=[Vector{Matrix{ComplexF64}}(undef,2) for _ in 1:Nq^2]
  
   tic=time()
    Threads.@threads for jk in 1:Nq^2
      Fk = FockMatrix[jk]
      for jk1 in 1:Nq^2
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
    
    for ja in 1:Nq^2
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
   
    
  
   
  
   for ja in 1:Nq^2, jb in 1:2
     FFF=eigen(single_Ham[ja][jb]+constq*HartreeMatrix-constq*FockMatrix[ja][jb])
     HF_eigenvalue[ja][jb]=real(FFF.values)
     HF_eigenvector[ja][jb]=FFF.vectors
   end


   
    #=
   htelement=sum(abs.(vec(HartreeMatrix)))
 println("htelement",htelement)



   fkelement=0
    for ja in 1:Nq^2, jb in 1:2
     fkelement+=sum(abs.(vec(FockMatrix[ja][jb])))/Nq^2
    end
  
  println("fkelement",fkelement)
=#
   bound=(sort(reduce(vcat,reduce(vcat,HF_eigenvalue)))[Parnum+1]+sort(reduce(vcat,reduce(vcat,HF_eigenvalue)))[Parnum])/2
   
   
    for ja in 1:Nq^2,vi in 1:2   
         for jd in eachindex(HF_eigenvalue[ja][vi])
            if HF_eigenvalue[ja][vi][jd]<bound
               NewDensityMatrix[ja][vi]+=HF_eigenvector[ja][vi][:,jd]*(HF_eigenvector[ja][vi][:,jd])'
            end
         end
         NewDensityMatrix[ja][vi]-=BG_DensityMatrix[ja][vi]
    end
  


   
    DeltaMatrix=NewDensityMatrix-input_DensityMatrix
  
    output_DensityMatrix=0.7*input_DensityMatrix+0.3*NewDensityMatrix

    l3=0.0
    for ja in 1:length(wave), jb in 1:Nq^2,vi in 1:2
     l3+=NewDensityMatrix[jb][vi][3*(ja-1)+1,3*(ja-1)+1]
    end
   println("l3=",l3)
    
    e1=0.0
    for ja in 1:Nq^2,vi in 1:2
      e1+=tr(DeltaMatrix[ja][vi]'*DeltaMatrix[ja][vi])
    end
    eout=real(e1)/(2*Nq^2)
   
    energy=0.0
    for ja in 1:Nq^2, vi in 1:2
      Energy_Matrix=single_Ham[ja][vi]+(constq*HartreeMatrix-constq*FockMatrix[ja][vi])/2
        energy+=tr(Energy_Matrix*input_DensityMatrix[ja][vi])
    end

    HF_chern=zeros(ComplexF64,4)
    dimension=3*length(wave)
    HF_chern[1]=calculate_chern(HF_eigenvector,Nq,dimension,wave,allowedq,1,2*length(wave)-1)
    HF_chern[2]=calculate_chern(HF_eigenvector,Nq,dimension,wave,allowedq,1,2*length(wave))
    HF_chern[3]=calculate_chern(HF_eigenvector,Nq,dimension,wave,allowedq,2,2*length(wave)-1)
    HF_chern[4]=calculate_chern(HF_eigenvector,Nq,dimension,wave,allowedq,2,2*length(wave))

    println("HFchern",HF_chern)
    


  
   return  eout,output_DensityMatrix,DeltaMatrix,HF_eigenvalue,HF_eigenvector,bound,real(energy)
end
  
  
#=
function Construct_HFmatrix(loop_dic::Dict{Vector{Int},Any},pathpoint::Vector{Float64},vi::Int64,allowedq::Vector{Vector{Int}},T1::Vector{Float64},T2::Vector{Float64},Nq::Int64,wave::Vector{Vector{Int64}},input_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},constq::Float64)::Matrix{ComplexF64}
  
  
  dimension=3*length(wave)
  HartreeMatrix=zeros(ComplexF64,dimension,dimension) 
  FockMatrix=zeros(ComplexF64,dimension,dimension) 
 
  


    
    for jk1 in 1:Nq^2
        dmk = input_DensityMatrix[jk1]
        q=allowedq[jk1][1]*T1+allowedq[jk1][2]*T2-pathpoint
       for (dg,loop_dic_dg) in loop_dic
           CoulF=CoulombMatrix_arbitrary(q+dg[1]*T1+dg[2]*T2,norm(T1))     
           for (gg2,loop_dic_dg_gg2) in loop_dic_dg          
           for g1g3 in loop_dic_dg_gg2
               FockMatrix[g1g3[5]:g1g3[6],gg2[3]:gg2[4]]+=dmk[vi][g1g3[3]:g1g3[4],gg2[5]:gg2[6]] .* CoulF
           end 
           end
       end    
    end

  

  Hartree_Density=zeros(ComplexF64,dimension,dimension)
  for jk1 in 1:Nq^2
   Hartree_Density+=input_DensityMatrix[jk1][1]+input_DensityMatrix[jk1][2]
  end

 

    

 
    for dg in keys(loop_dic)
        CoulH=CoulombMatrix(dg,T1,T2)
        for gg2 in keys(loop_dic[dg])                    
        for g1g3 in loop_dic[dg][gg2]            
             HartreeMatrix[gg2[5]:gg2[6],gg2[3]:gg2[4]]+=diagm(CoulH*diag(Hartree_Density[g1g3[3]:g1g3[4],g1g3[5]:g1g3[6]],0))                         
        end 
        end

    end    
 


 return  constq*(HartreeMatrix-FockMatrix)
end

=#













function iteration_loop(initial_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},BG_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},allowedq::Vector{Vector{Int64}},T1::Vector{Float64},T2::Vector{Float64},Nq::Int,wave::Vector{Vector{Int64}},single_Ham::Vector{Vector{Matrix{ComplexF64}}},constq::Float64,holenum::Int64)::Tuple{Vector{Vector{Vector{Matrix{ComplexF64}}}},Vector{Vector{Vector{Matrix{ComplexF64}}}},Vector{Vector{Vector{Float64}}},Vector{Vector{Matrix{ComplexF64}}},Float64,Float64,Vector{ComplexF64}}
    eout=1.0
    itcount=0
    loop_dic=construct_loop_dic(wave)
    HF_eigenvalue=[Vector{Vector{Float64}}(undef,2) for _ in 1:Nq^2]
    HF_eigenvector=[Vector{Matrix{ComplexF64}}(undef,2) for _ in 1:Nq^2]
    DIIS_input_DensityMatrix=Vector{Vector{Vector{Matrix{ComplexF64}}}}(undef,3)
    DIIS_input_DeltaMatrix=Vector{Vector{Vector{Matrix{ComplexF64}}}}(undef,3)
    input_DensityMatrix=initial_DensityMatrix
    bad_count=0
    Parnum=2*Nq^2*(2*length(wave))-holenum*Nq^2
    bound=0.0
    energy=0.0

    while (eout>1*10^-13) || (bad_count<4)
      if  eout<1*10^-13
        bad_count+=1
      end
      tic=time()
      eout,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue, HF_eigenvector,bound,energy=Construct_DensityMatrix(loop_dic,allowedq,T1,T2,Nq,wave,input_DensityMatrix,BG_DensityMatrix,single_Ham,constq,Parnum)
      DIIS_input_DensityMatrix[mod(itcount,3)+1]=input_DensityMatrix
      input_DensityMatrix=output_DensityMatrix
      itcount+=1
      toc=time()
      println(toc-tic,"eout=$eout")
      flush(stdout)
     
    
    end
    
    HF_chern=zeros(ComplexF64,4)
    dimension=3*length(wave)
    HF_chern[1]=calculate_chern(HF_eigenvector,Nq,dimension,wave,allowedq,1,2*length(wave)-1)
    HF_chern[2]=calculate_chern(HF_eigenvector,Nq,dimension,wave,allowedq,1,2*length(wave))
    HF_chern[3]=calculate_chern(HF_eigenvector,Nq,dimension,wave,allowedq,2,2*length(wave)-1)
    HF_chern[4]=calculate_chern(HF_eigenvector,Nq,dimension,wave,allowedq,2,2*length(wave))


  return DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue, HF_eigenvector,bound,energy,HF_chern

end










function Densitymap(a1m::Vector{Float64},a2m::Vector{Float64},wave::Vector{Vector{Int}},input_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}})::Tuple{Array{ComplexF64},Array{ComplexF64},Array{ComplexF64}}
  dimension=3*length(wave)  
     N3=50
    xgrid=zeros(Float64,N3,N3)
    ygrid=zeros(Float64,N3,N3)
    zgrid=zeros(ComplexF64,N3,N3,2,3)
    Total_Density=[zeros(ComplexF64,dimension,dimension) for _ in 1:2]
    for jk1 in 1:Nq^2
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




function calculate_chern(eigenvector_matrix,Nq,dimension,wave,allowedq,vi,bi)
  chern_eigenvector=zeros(ComplexF64,dimension,Nq+1,Nq+1)
  for ja in 1:Nq+1, jb in 1:Nq+1
    if (ja<Nq+1) && (jb<Nq+1)
      pos=findfirst(item->item==[ja-1,jb-1],allowedq)
      chern_eigenvector[:,ja,jb]=eigenvector_matrix[pos][vi][:,bi]
    else
      G1=[ja-1,jb-1]-[mod(ja-1,Nq),mod(jb-1,Nq)]
      k=[mod(ja-1,Nq),mod(jb-1,Nq)]
      pos=findfirst(item->item==k,allowedq)
      for jc in eachindex(wave)
        pos_2=findfirst(item->item==G1+wave[jc],wave)
        if pos_2≠nothing
          chern_eigenvector[3*(jc-1)+1:3*jc,ja,jb]=eigenvector_matrix[pos][vi][3*(pos_2-1)+1:3*pos_2,bi]
        end
      end
    end

  end

  Uonelink=zeros(ComplexF64,Nq,Nq+1)
  Utwolink=zeros(ComplexF64,Nq+1,Nq)

  
  for ja in 1:Nq, jb in 1:Nq+1
   
     Uonelink[ja,jb]=dot(chern_eigenvector[:,ja,jb],chern_eigenvector[:,ja+1,jb])/abs(dot(chern_eigenvector[:,ja,jb],chern_eigenvector[:,ja+1,jb]))
  end

  
  
  for ja in 1:Nq+1, jb in 1:Nq

   Utwolink[ja,jb]=dot(chern_eigenvector[:,ja,jb],chern_eigenvector[:,ja,jb+1])/abs(dot(chern_eigenvector[:,ja,jb],chern_eigenvector[:,ja,jb+1]))
  end
  
  Flink=zeros(ComplexF64,Nq,Nq)
  for ja in 1:Nq, jb in 1:Nq
   Flink[ja,jb]=log(Uonelink[ja,jb]*Utwolink[ja+1,jb]/(Uonelink[ja,jb+1]*Utwolink[ja,jb]))
  end
  chern=sum(Flink)/(2*π*im)



  return chern
end








