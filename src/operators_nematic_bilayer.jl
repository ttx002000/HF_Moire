using LinearAlgebra
using Arpack
using Combinatorics

using Random

#The unit cell length is different here


 
 
 
function Coulomb(qab::Float64,k::Vector{Int64},L1::Int64,L2::Int64)::Float64
  
   D=10
   d=0.73
   Ld=abs(L1-L2)
   
   return k==[0,0] ? (D-Ld*d) : ((exp(-Ld*d*qab)-2/(1+exp(2*D*qab)))/qab)
   
end


function CoulombMatrix(k::Vector{Int64},T1::Vector{Float64},T2::Vector{Float64})::Matrix{Float64}
  
  Cmatrix=Matrix{Float64}(undef,2,2)
  qab=norm(k[1]*T1+k[2]*T2)
  for L1 in 1:2, L2 in 1:2
    Cmatrix[L1,L2]=Coulomb(qab,k,L1,L2)
  end
  return Cmatrix
  
end



function generate_seed(generate_num,wavenum)
  seed_tunnel=[zeros(Float64,2*wavenum,2*wavenum) for _ in 1:2]
  if generate_num==1
    for ja in 1:wavenum
    seed_tunnel[1][2*(ja-1)+1,2*(ja-1)+2]+=1
    end
  end

  if generate_num==2
    for ja in 1:wavenum
    seed_tunnel[2][2*(ja-1)+1,2*(ja-1)+2]+=1
    end
  end
  
  if generate_num==3
    for ja in 1:wavenum
    seed_tunnel[1][2*(ja-1)+1,2*(ja-1)+2]+=1
    seed_tunnel[2][2*(ja-1)+1,2*(ja-1)+2]+=1
    end
  end


  seed_tunnel[1]=seed_tunnel[1]+seed_tunnel[1]'
  seed_tunnel[2]=seed_tunnel[2]+seed_tunnel[2]'
  return 2*seed_tunnel
end





function triangle_initial_Densitymatrix_control(parameters::Vector{Float64},Nq::Int64,seednum::Int64)
 

  mt=parameters[1]
  mb=parameters[2]
  Vt=parameters[3]
  ϕt=parameters[4]/180*π
  Vb=parameters[5]
  ϕb=parameters[6]/180*π
  ϵr=parameters[7]
  Eg=parameters[8]
  θ=parameters[9]/180*π
  w=parameters[10]
  omega=exp(i*2*π/3*parametersp[11])






  alattice=0.347# In units of NM
  blattice=4*π/(√3*alattice)

  
  
  bm=2*blattice*sin(θ/2)
  am=4*π/(√3*bm);
  constt=-38.09981949*1/(mt)
  constb=38.09981949*1/(mb) # the 38 is \hbar^2/(2me*nm^2) in units of meV
  constq=1/(Nq^2*3^(1/2)/2*am^2)*1/ϵr*9047.5636
  
  
  b1=bm*[1,0]
  b2=bm*[-1/2,√3/2]
  

  a1m=am*[√3/2,1/2]
  a2m=am*[0,1]
  
  
  T1=b1/(Nq)
  T2=b2/(Nq)
  
  κt=bm*[-1/2,1/(2√3)]
 κb=bm*[-1/2,-1/(2√3)]

  
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
  dimension=2*length(wave)
 
  




  single_Ham=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for ja in 1:Nq^2]
  single_MoirePo=zeros(ComplexF64,dimension,dimension)
  tunnel=zeros(ComplexF64,dimension,dimension)
  M_tunnel=zeros(ComplexF64,dimension,dimension)
  
  single_eigenvalue=[[zeros(Float64,dimension) for _ in 1:2] for _ in 1:Nq^2]
  single_eigenvector=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
  seed_eigenvector=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
  uncoupled_eigenvector=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
  

  
  for jc in eachindex(wave)
    pos=findfirst(item->item==wave[jc]-b1T,wave)
    if pos≠nothing
      single_MoirePo[2*(jc-1)+1:2*jc,2*(pos-1)+1:2*pos]=diagm([Vt*exp(im*ϕt),Vb*exp(im*ϕb)])
    end
    
  
    pos=findfirst(item->item==wave[jc]+b2T+b1T,wave)
    if pos≠nothing
      single_MoirePo[2*(jc-1)+1:2*jc,2*(pos-1)+1:2*pos]=diagm([Vt*exp(im*ϕt),Vb*exp(im*ϕb)])
    end
  
    pos=findfirst(item->item==wave[jc]-b2T,wave)
    if pos≠nothing
      single_MoirePo[2*(jc-1)+1:2*jc,2*(pos-1)+1:2*pos]=diagm([Vt*exp(im*ϕt),Vb*exp(im*ϕb)])
    end
  end
  
  
  single_MoirePo=single_MoirePo+single_MoirePo'
  
  for jc in eachindex(wave)
    tunnel[2*(jc-1)+1,2*(jc-1)+2]=w
    
  
    pos=findfirst(item->item==wave[jc]-b2T-b1T,wave)
    if pos≠nothing
      tunnel[2*(jc-1)+1,2*(pos-1)+2]=w*omega
      
    end
  
    pos=findfirst(item->item==wave[jc]-b2T,wave)
    if pos≠nothing
      tunnel[2*(jc-1)+1,2*(pos-1)+2]=w*omega'
     
    end
  end
  
  tunnel=tunnel+tunnel'
  
  for jc in eachindex(wave)
    M_tunnel[2*(jc-1)+1,2*(jc-1)+2]=w
   
  
    pos=findfirst(item->item==wave[jc]+b2T+b1T,wave)
    if pos≠nothing
      M_tunnel[2*(jc-1)+1,2*(pos-1)+2]=w*omega'
     
    end
  
    pos=findfirst(item->item==wave[jc]+b2T,wave)
    if pos≠nothing
      M_tunnel[2*(jc-1)+1,2*(pos-1)+2]=w*omega
    
    end
  end
  
  M_tunnel=M_tunnel+M_tunnel'
  
  
  seed_tunnel=generate_seed(seednum,length(wave))
  
  
  Threads.@threads for ja in 1:Nq^2
    
    
    k=allowedq[ja][1]*T1+allowedq[ja][2]*T2
  
    for jb in eachindex(wave)
     single_Ham[ja][1][2*(jb-1)+1,2*(jb-1)+1]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2-κt)^2*constt
     single_Ham[ja][1][2*(jb-1)+2,2*(jb-1)+2]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2-κb)^2*constb+Eg
    end
  
    for jb in eachindex(wave)
      single_Ham[ja][2][2*(jb-1)+1,2*(jb-1)+1]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2+κt)^2*constt
      single_Ham[ja][2][2*(jb-1)+2,2*(jb-1)+2]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2+κb)^2*constb+Eg
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

      
    FFF=eigen(single_Ham[ja][1]+seed_tunnel[1])
     seed_eigenvector[ja][1]=FFF.vectors

    FFF=eigen(single_Ham[ja][2]+seed_tunnel[2])
    seed_eigenvector[ja][2]=FFF.vectors
  
  end
  




   
  BG_DensityMatrix=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
  for ja in 1:Nq^2, vi in 1:2, jb in 1:length(wave)
    BG_DensityMatrix[ja][vi]+=(uncoupled_eigenvector[ja][vi][:,jb]*(uncoupled_eigenvector[ja][vi][:,jb])') 
  end
  
  input_DensityMatrix=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
  

 
  for ja in 1:Nq^2, vi in 1:1, jb in 1:length(wave)
      input_DensityMatrix[ja][vi]+=(seed_eigenvector[ja][vi][:,jb]*(seed_eigenvector[ja][vi][:,jb])') 
  end
    
  for ja in 1:Nq^2,vi in 2:2, jb in 1:length(wave)-1
      input_DensityMatrix[ja][vi]+=(seed_eigenvector[ja][vi][:,jb]*(seed_eigenvector[ja][vi][:,jb])') 
  end


 
    
 
  input_DensityMatrix=input_DensityMatrix-BG_DensityMatrix
 

  
 return  wave, input_DensityMatrix, BG_DensityMatrix, single_Ham, single_eigenvalue,allowedq, T1, T2, a1m, a2m,constq
    
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
        
        loop_dic[wave[g2index]-wave[gindex]][[gindex,g2index,2*gindex-1,2*gindex,2*g2index-1,2*g2index]]=Vector{Int64}[]
        
         
    end
    
    for gindex in eachindex(wave), g1index in eachindex(wave),g2index in eachindex(wave)
        g3=wave[gindex]+wave[g1index]-wave[g2index]
        if haskey(g_dic,wave[gindex]+wave[g1index]-wave[g2index])
            push!(loop_dic[wave[g2index]-wave[gindex]][[gindex,g2index,2*gindex-1,2*gindex,2*g2index-1,2*g2index]],[g1index,g_dic[g3],2*g1index-1,2*g1index,2*g_dic[g3]-1,2*g_dic[g3]])
        end
    end
    
   return loop_dic

end






function Construct_DensityMatrix(loop_dic::Dict{Vector{Int},Any},allowedq::Vector{Vector{Int}},T1::Vector{Float64},T2::Vector{Float64},Nq::Int64,wave::Vector{Vector{Int64}},input_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},BG_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},single_Ham::Vector{Vector{Matrix{ComplexF64}}},constq::Float64,Parnum::Int64)::Tuple{Float64,Vector{Vector{Matrix{ComplexF64}}},Vector{Vector{Matrix{ComplexF64}}},Vector{Vector{Vector{Float64}}},Float64}
  
 
    dimension=2*length(wave)
    HartreeMatrix=zeros(ComplexF64,dimension,dimension)
    FockMatrix=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
  
    output_DensityMatrix=[Vector{Matrix{ComplexF64}}(undef,2) for _ in 1:Nq^2]
    DeltaMatrix=[Vector{Matrix{ComplexF64}}(undef,2) for _ in 1:Nq^2]
    NewDensityMatrix=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
    HF_eigenvalue=[Vector{Vector{Float64}}(undef,2) for _ in 1:Nq^2]
    HF_eigenvector=[Vector{Matrix{ComplexF64}}(undef,2) for _ in 1:Nq^2]
  
   
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
  
   bound=sort(reduce(vcat,reduce(vcat,HF_eigenvalue)))[Parnum+1]
   
  
    for ja in 1:Nq^2,vi in 1:2     
         for jd in eachindex(HF_eigenvalue[ja][vi])
            if HF_eigenvalue[ja][vi][jd]<bound
               NewDensityMatrix[ja][vi]+=HF_eigenvector[ja][vi][:,jd]*(HF_eigenvector[ja][vi][:,jd])'
            end
         end
    end
  
    NewDensityMatrix=NewDensityMatrix-BG_DensityMatrix
    DeltaMatrix=NewDensityMatrix-input_DensityMatrix
    output_DensityMatrix=0.4*input_DensityMatrix+0.6*NewDensityMatrix
    
    e1=0.0
    for ja in 1:Nq^2,vi in 1:2
      e1+=tr(DeltaMatrix[ja][vi]'*DeltaMatrix[ja][vi])
    end
    eout=real(e1)/(2*Nq^2)
   
   
  
   return  eout,output_DensityMatrix,DeltaMatrix,HF_eigenvalue,bound
end
  
  




function iteration_loop(initial_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},BG_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},allowedq::Vector{Vector{Int64}},T1::Vector{Float64},T2::Vector{Float64},Nq::Int,wave::Vector{Vector{Int64}},single_Ham::Vector{Vector{Matrix{ComplexF64}}},constq::Float64,holenum::Int64)::Tuple{Vector{Vector{Vector{Matrix{ComplexF64}}}},Vector{Vector{Vector{Matrix{ComplexF64}}}},Vector{Vector{Vector{Float64}}},Float64}
    eout=1.0
    itcount=0
    loop_dic=construct_loop_dic(wave)
    HF_eigenvalue=[Vector{Vector{Float64}}(undef,2) for _ in 1:Nq^2]
    DIIS_input_DensityMatrix=Vector{Vector{Vector{Matrix{ComplexF64}}}}(undef,3)
    DIIS_input_DeltaMatrix=Vector{Vector{Vector{Matrix{ComplexF64}}}}(undef,3)
    input_DensityMatrix=initial_DensityMatrix
    bad_count=0
    Parnum=2*Nq^2*(length(wave))-holenum*Nq^2
    bound=0.0


    while (eout>1*10^-14) || (bad_count<4)
      if  eout<1*10^-14 
        bad_count+=1
      end
      tic=time()
      eout,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue,bound=Construct_DensityMatrix(loop_dic,allowedq,T1,T2,Nq,wave,input_DensityMatrix,BG_DensityMatrix,single_Ham,constq,Parnum)
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
    while (eout>10^-14) || (bad_count<4)
        if  eout<1*10^-14 
            bad_count+=1
        end
        tic=time()
        Bmatrix=zeros(ComplexF64,4,4)
        for ja in 1:3
         Bmatrix[ja,4]=1
         Bmatrix[4,ja]=1
        end
    
        for ja in 1:3,jb in 1:3
            for jc in 1:Nq^2, vi in 1:2
               Bmatrix[ja,jb]+=tr((DIIS_input_DeltaMatrix[ja][jc][vi])'*(DIIS_input_DeltaMatrix[jb][jc][vi]))
            end
        end
        coeff=inv(Bmatrix)*[0;0;0;1]
       
        dmk=coeff[1]*(DIIS_input_DensityMatrix[1]+DIIS_input_DeltaMatrix[1])+coeff[2]*(DIIS_input_DensityMatrix[2]+DIIS_input_DeltaMatrix[2])+coeff[3]*(DIIS_input_DensityMatrix[3]+DIIS_input_DeltaMatrix[3])
        eout,_,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue,bound=Construct_DensityMatrix(loop_dic,allowedq,T1,T2,Nq,wave,dmk,BG_DensityMatrix,single_Ham,constq,Parnum)
        DIIS_input_DensityMatrix[mod(itcount,3)+1]=dmk
        itcount+=1

        toc=time()
        println(toc-tic,"eout=$eout")
        flush(stdout)
    end
    =#
 


  return DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,bound

end





function iteration_loop_control(initial_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},BG_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},allowedq::Vector{Vector{Int64}},T1::Vector{Float64},T2::Vector{Float64},Nq::Int,wave::Vector{Vector{Int64}},single_Ham::Vector{Vector{Matrix{ComplexF64}}},constq::Float64,holenum::Int64,seednum::Int64)::Tuple{Vector{Vector{Vector{Matrix{ComplexF64}}}},Vector{Vector{Vector{Matrix{ComplexF64}}}},Vector{Vector{Vector{Float64}}},Float64}
  eout=1.0
  itcount=0
  loop_dic=construct_loop_dic(wave)
  HF_eigenvalue=[Vector{Vector{Float64}}(undef,2) for _ in 1:Nq^2]
  DIIS_input_DensityMatrix=Vector{Vector{Vector{Matrix{ComplexF64}}}}(undef,3)
  DIIS_input_DeltaMatrix=Vector{Vector{Vector{Matrix{ComplexF64}}}}(undef,3)
  input_DensityMatrix=initial_DensityMatrix
  bad_count=0
  Parnum=2*Nq^2*(length(wave))-holenum*Nq^2
  bound=0.0
  seed_tlmatrix=[generate_seed(seednum,length(wave)) for _ in 1:Nq^2]
  

  while itcount<20
  
    tic=time()
    eout,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue,bound=Construct_DensityMatrix(loop_dic,allowedq,T1,T2,Nq,wave,input_DensityMatrix,BG_DensityMatrix,single_Ham+seed_tlmatrix,constq,Parnum)
    DIIS_input_DensityMatrix[mod(itcount,3)+1]=input_DensityMatrix
    input_DensityMatrix=output_DensityMatrix
    itcount+=1
    toc=time()
    println(toc-tic,"eout=$eout")
    flush(stdout)
   
  
  end






  while (eout>1*10^-14) || (bad_count<4)
    if  eout<1*10^-14 
      bad_count+=1
    end
    tic=time()
    eout,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue,bound=Construct_DensityMatrix(loop_dic,allowedq,T1,T2,Nq,wave,input_DensityMatrix,BG_DensityMatrix,single_Ham,constq,Parnum)
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
  while (eout>10^-14) || (bad_count<4)
      if  eout<1*10^-14 
          bad_count+=1
      end
      tic=time()
      Bmatrix=zeros(ComplexF64,4,4)
      for ja in 1:3
       Bmatrix[ja,4]=1
       Bmatrix[4,ja]=1
      end
  
      for ja in 1:3,jb in 1:3
          for jc in 1:Nq^2, vi in 1:2
             Bmatrix[ja,jb]+=tr((DIIS_input_DeltaMatrix[ja][jc][vi])'*(DIIS_input_DeltaMatrix[jb][jc][vi]))
          end
      end
      coeff=inv(Bmatrix)*[0;0;0;1]
     
      dmk=coeff[1]*(DIIS_input_DensityMatrix[1]+DIIS_input_DeltaMatrix[1])+coeff[2]*(DIIS_input_DensityMatrix[2]+DIIS_input_DeltaMatrix[2])+coeff[3]*(DIIS_input_DensityMatrix[3]+DIIS_input_DeltaMatrix[3])
      eout,_,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue,bound=Construct_DensityMatrix(loop_dic,allowedq,T1,T2,Nq,wave,dmk,BG_DensityMatrix,single_Ham,constq,Parnum)
      DIIS_input_DensityMatrix[mod(itcount,3)+1]=dmk
      itcount+=1

      toc=time()
      println(toc-tic,"eout=$eout")
      flush(stdout)
  end
  =#



return DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,bound

end








function Densitymap(a1m::Vector{Float64},a2m::Vector{Float64},wave::Vector{Vector{Int}},input_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}})::Tuple{Array{ComplexF64},Array{ComplexF64},Array{ComplexF64}}
  dimension=2*length(wave)  
     N3=50
    xgrid=zeros(Float64,N3,N3)
    ygrid=zeros(Float64,N3,N3)
    zgrid=zeros(ComplexF64,N3,N3,2,2)#first valley,second layer
    Total_Density=[zeros(ComplexF64,dimension,dimension) for _ in 1:2]
    for jk1 in 1:Nq^2
       Total_Density+=input_DensityMatrix[jk1]
    end
    
    for ja in 1:N3, jb in 1:N3
      xgrid[ja,jb]=(ja/N3*a1m[1]+jb/N3*a2m[1])
      ygrid[ja,jb]=(ja/N3*a1m[2]+jb/N3*a2m[2])
      rvec=ja/N3*a1m+jb/N3*a2m
      for vi in 1:2, Li in 1:2, jc in eachindex(wave), jd in eachindex(wave)
        gvec=[T1 T2]*(wave[jc]-wave[jd])
       zgrid[ja,jb,vi,Li]+=Total_Density[vi][2*(jc-1)+Li,2*(jd-1)+Li]*exp(im*(gvec[1]*rvec[1]+gvec[2]*rvec[2]))
      end
    
    end

    return xgrid,ygrid,zgrid
end





function Construct_HFmatrix(loop_dic::Dict{Vector{Int},Any},pathpoint::Vector{Int64},vi::Int64,allowedq::Vector{Vector{Int}},T1::Vector{Float64},T2::Vector{Float64},Nq::Int64,wave::Vector{Vector{Int64}},input_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},constq::Float64)::Matrix{ComplexF64}
  
  
  dimension=2*length(wave)
  HartreeMatrix=zeros(ComplexF64,dimension,dimension) 
  FockMatrix=zeros(ComplexF64,dimension,dimension) 
 
  


    
    for jk1 in 1:Nq^2
        dmk = input_DensityMatrix[jk1]
        q=allowedq[jk1]-pathpoint
       for (dg,loop_dic_dg) in loop_dic
           CoulF=CoulombMatrix(q+dg,T1,T2)     
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








function calculate_energy(Nq::Int,wave::Vector{Vector{Int}},input_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},constq::Float64,T1::Vector{Float64},T2::Vector{Float64},allowedq::Vector{Vector{Int64}},single_Ham::Vector{Vector{Matrix{ComplexF64}}})

  
   
    dimension=2*length(wave)
    
    
    loop_dic=construct_loop_dic(wave)
    
    Energy_Matrix=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
    HF_vectors=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
    
    Threads.@threads for ja in eachindex(allowedq)
        for vi in 1:2
      HFmatrix=Construct_HFmatrix(loop_dic,allowedq[ja],vi,allowedq,T1,T2,Nq,wave,input_DensityMatrix,constq)
      Energy_Matrix[ja][vi]=1/2*HFmatrix+single_Ham[ja][vi]
      HF_vectors[ja][vi]=eigen(HFmatrix+single_Ham[ja][vi]).vectors
    end
    end
    
    
   energy=0
   for ja in 1:Nq^2, vi in 1:2
       energy+=tr(Energy_Matrix[ja][vi]*input_DensityMatrix[ja][vi])
   end


   return energy,HF_vectors

end