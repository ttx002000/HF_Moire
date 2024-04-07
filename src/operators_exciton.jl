using LinearAlgebra
using Arpack
using Combinatorics
using LinearAlgebra
using Random




 
 
 
function Coulomb(k::Vector{Int64},T1::Vector{Float64},T2::Vector{Float64})::Float64
  
   return k==[0,0] ? 0.0 : 1/norm(k[1]*T1+k[2]*T2)
   
end







function triangle_initial_Densitymatrix(parameters::Vector{Float64},Nq::Int64)
 

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







    alattice=0.3# In units of NM
    blattice=4*π/(√3*alattice)
 
    
    
    bm=√3*2*blattice*sin(θ/2)
    am=4*π/(√3*bm);
    constt=-38.09981949*1/(mt)
    constm=38.09981949*1/(mm)
    constb=-38.09981949*1/(mb) # the 38 is \hbar^2/(2me*nm^2) in units of meV
    constq=1/(Nq^2*3^(1/2)/2*am^2)*1/ϵr*9047.5636
    
    
    b1=bm*[1,0]
    b2=bm*[-1/2,√3/2]
    

    a1m=am*[1/2,√3/2]
    a2m=am*[1,0]
    
    
    T1=b1/(Nq)
    T2=b2/(Nq)
    
    
    
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
    dimension=length(wave)
   
    




    single_Ham=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for ja in 1:Nq^2]
    single_MoirePo=zeros(ComplexF64,dimension,dimension)
    tunnel=zeros(ComplexF64,dimension,dimension)
    M_tunnel=zeros(ComplexF64,dimension,dimension)
    
    single_eigenvalue=[[zeros(Float64,dimension) for _ in 1:2] for _ in 1:Nq^2]
    single_eigenvector=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
    uncoupled_eigenvector=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
    

    
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
      tunnel[3*(jc-1)+1,3*(jc-1)+2]=w
      tunnel[3*(jc-1)+3,3*(jc-1)+2]=w
    
      pos=findfirst(item->item==wave[jc]-b2T-b1T,wave)
      if pos≠nothing
        tunnel[3*(jc-1)+1,3*(pos-1)+2]=w
        tunnel[3*(jc-1)+3,3*(pos-1)+2]=w
      end
    
      pos=findfirst(item->item==wave[jc]-b2T,wave)
      if pos≠nothing
        tunnel[3*(jc-1)+1,3*(pos-1)+2]=w
        tunnel[3*(jc-1)+3,3*(pos-1)+2]=w
      end
    end
    
    tunnel=tunnel+tunnel'
    
    for jc in eachindex(wave)
      M_tunnel[3*(jc-1)+1,3*(jc-1)+2]=w
      M_tunnel[3*(jc-1)+3,3*(jc-1)+2]=w
    
      pos=findfirst(item->item==wave[jc]+b2T+b1T,wave)
      if pos≠nothing
        M_tunnel[3*(jc-1)+1,3*(pos-1)+2]=w
        M_tunnel[3*(jc-1)+3,3*(pos-1)+2]=w
      end
    
      pos=findfirst(item->item==wave[jc]+b2T,wave)
      if pos≠nothing
        M_tunnel[3*(jc-1)+1,3*(pos-1)+2]=w
        M_tunnel[3*(jc-1)+3,3*(pos-1)+2]=w
      end
    end
    
    M_tunnel=M_tunnel+M_tunnel'
    
    
    
    
    
    Threads.@threads for ja in 1:Nq^2
      
      
      k=allowedq[ja][1]*T1+allowedq[ja][2]*T2
    
      for jb in eachindex(wave)
       single_Ham[ja][1][3*(jb-1)+1,3*(jb-1)+1]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2-κt)^2*constt
       single_Ham[ja][1][3*(jb-1)+2,3*(jb-1)+2]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2-κm)^2*constm+Eg
       single_Ham[ja][1][3*(jb-1)+3,3*(jb-1)+3]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2-κb)^2*constb
      end
    
      for jb in eachindex(wave)
        single_Ham[ja][2][3*(jb-1)+1,3*(jb-1)+1]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2+κt)^2*constt
        single_Ham[ja][2][3*(jb-1)+2,3*(jb-1)+2]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2+κm)^2*constm+Eg
        single_Ham[ja][2][3*(jb-1)+3,3*(jb-1)+3]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2+κb)^2*constb
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
    



    input_DensityMatrix=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nq^2]
    for ja in 1:Nq^2, vi in 1:2
        A=randn(ComplexF64,dimension,dimension)
       input_DensityMatrix[ja][vi]+=A+A'
    end
     
    BG_DensityMatrix=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
    for ja in 1:Nq^2, vi in 1:2, jb in 1:length(wave)*2
      BG_DensityMatrix[ja][vi]+=(uncoupled_eigenvector[ja][vi][:,jb]*(uncoupled_eigenvector[ja][vi][:,jb])') 
    end
    
    
    
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






function Construct_DensityMatrix(loop_dic::Dict{Vector{Int},Any},allowedq::Vector{Vector{Int}},T1::Vector{Float64},T2::Vector{Float64},Nq::Int64,wave::Vector{Vector{Int64}},input_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},BG_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},single_Ham::Vector{Vector{Matrix{ComplexF64}}},constq::Float64,Parnum::Int64)::Tuple{Float64,Vector{Vector{Matrix{ComplexF64}}},Vector{Vector{Matrix{ComplexF64}}},Vector{Vector{Vector{Float64}}},Float64}
  
 
    dimension=3*length(wave)
    HartreeMatrix=zeros(ComplexF64,dimension,dimension)
    FockMatrix=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
  
    output_DensityMatrix=[Vector{Matrix{ComplexF64}}(undef,2) for _ in 1:Nq^2]
    DeltaMatrix=[Vector{Matrix{ComplexF64}}(undef,2) for _ in 1:Nq^2]
    NewDensityMatrix=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
    HF_eigenvalue=[Vector{Vector{Float64}}(undef,2) for _ in 1:Nq^2]
    HF_eigenvector=[Vector{Matrix{ComplexF64}}(undef,2) for _ in 1:Nq^2]
  
   tic=time()
    for jk in 1:Nq^2
      Fk = FockMatrix[jk]
      for jk1 in 1:Nq^2
          dmk = input_DensityMatrix[jk1]
          q=allowedq[jk1]-allowedq[jk]
         for (dg,loop_dic_dg) in loop_dic
             CoulF=Coulomb(q+dg,T1,T2)     
             for (gg2,loop_dic_dg_gg2) in loop_dic_dg         
             for g1g3 in loop_dic_dg_gg2
                    Fk[1][g1g3[5]:g1g3[6],gg2[3]:gg2[4]]+=dmk[1][g1g3[3]:g1g3[4],gg2[5]:gg2[6]]*CoulF
                 Fk[2][g1g3[5]:g1g3[6],gg2[3]:gg2[4]]+=dmk[2][g1g3[3]:g1g3[4],gg2[5]:gg2[6]]*CoulF   
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
  
    
  
   Identity=Matrix{Float64}(I,3,3)   
   
     for dg in keys(loop_dic)
            CoulH=Coulomb(dg,T1,T2)
            for gg2 in keys(loop_dic[dg])                   
            for g1g3 in loop_dic[dg][gg2]           
                HartreeMatrix[gg2[5]:gg2[6],gg2[3]:gg2[4]]+=tr(Hartree_Density[g1g3[3]:g1g3[4],g1g3[5]:g1g3[6]])*CoulH*Identity               
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
    output_DensityMatrix=0.0*input_DensityMatrix+1.0*NewDensityMatrix
    
    e1=0.0
    for ja in 1:Nq^2,vi in 1:2
      e1+=tr(DeltaMatrix[ja][vi]'*DeltaMatrix[ja][vi])
    end
    eout=real(e1)
   
   
  
   return  eout,output_DensityMatrix,DeltaMatrix,HF_eigenvalue,bound
end
  
  




function iteration_loop(initial_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},BG_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},allowedq::Vector{Vector{Int64}},T1::Vector{Float64},T2::Vector{Float64},Nq::Int,wave::Vector{Vector{Int64}},single_Ham::Vector{Matrix{ComplexF64}},constq::Float64)::Tuple{Vector{Vector{Matrix{ComplexF64}}},Vector{Vector{Matrix{ComplexF64}}},Vector{Vector{Float64}},Float64}
    eout=1.0
    itcount=0
    loop_dic=construct_loop_dic(wave)
    HF_eigenvalue=[Vector{Vector{Float64}}(undef,2) for _ in 1:Nq^2]
    DIIS_input_DensityMatrix=Vector{Vector{Vector{Matrix{ComplexF64}}}}(undef,3)
    DIIS_input_DeltaMatrix=Vector{Vector{Vector{Matrix{ComplexF64}}}}(undef,3)
    input_DensityMatrix=initial_DensityMatrix
    bad_count=0

    while (eout>1*10^-9) || (bad_count<4)
      if  eout<1*10^-9 
        bad_count+=1
      end
      tic=time()
      eout,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue=Construct_DensityMatrix(loop_dic,allowedq,T1,T2,Nq,wave,input_DensityMatrix,BG_DensityMatrix,single_Ham,constq,Parnum)
      DIIS_input_DensityMatrix[mod(itcount,3)+1]=input_DensityMatrix
      input_DensityMatrix=output_DensityMatrix
      itcount+=1
      toc=time()
      println(toc-tic,"eout=$eout")
      flush(stdout)
     
    
    end
    
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
            for jc in 1:Nq^2
               Bmatrix[ja,jb]+=tr((DIIS_input_DeltaMatrix[ja][jc])'*(DIIS_input_DeltaMatrix[jb][jc]))
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

 


  return DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,bound

end


function Densitymap(a1m::Vector{Float64},a2m::Vector{Float64},wave::Vector{Vector{Int}},input_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}})::Tuple{Array{ComplexF64},Array{ComplexF64},Array{ComplexF64}}
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





function Construct_HFmatrix(loop_dic::Dict{Vector{Int},Any},pathpointindex::Int64,pathpoint::Vector{Int64},allowedq::Vector{Vector{Int}},T1::Vector{Float64},T2::Vector{Float64},Nq::Int64,wave::Vector{Vector{Int64}},input_DensityMatrix::Vector{Matrix{ComplexF64}},constq::Float64,chern_overlapmatrix::Array{ComplexF64,4})::Matrix{ComplexF64}
  
  
  dimension=length(wave)
  HartreeMatrix=zeros(ComplexF64,dimension,dimension) 
  FockMatrix=zeros(ComplexF64,dimension,dimension) 
 
  


    
    for jk1 in 1:Nq^2
        dmk = input_DensityMatrix[jk1]
        q=allowedq[jk1]-pathpoint
       for (dg,loop_dic_dg) in loop_dic
           CoulF1=Coulomb(q+dg,T1,T2)     
           for (gg2,loop_dic_dg_gg2) in loop_dic_dg          
               CoulF=CoulF1*chern_overlapmatrix[jk1,gg2[2],pathpointindex,gg2[1]]
           for g1g3 in loop_dic_dg_gg2
               FockMatrix[g1g3[2],gg2[1]]+=dmk[g1g3[1],gg2[2]]*CoulF*chern_overlapmatrix[pathpointindex,g1g3[2],jk1,g1g3[1]]
           end 
           end
   
       end    
    end

  

  Hartree_Density=zeros(ComplexF64,dimension,dimension)
  for jk1 in 1:Nq^2
   Hartree_Density+=input_DensityMatrix[jk1] .* transpose(chern_overlapmatrix[jk1,:,jk1,:])
  end

 

    
 
    for dg in keys(loop_dic)
        CoulH1=Coulomb(dg,T1,T2)
        for gg2 in keys(loop_dic[dg])          
            CoulH=CoulH1*chern_overlapmatrix[pathpointindex,gg2[2],pathpointindex,gg2[1]]            
        for g1g3 in loop_dic[dg][gg2]           
            HartreeMatrix[gg2[2],gg2[1]]+=Hartree_Density[g1g3[1],g1g3[2]]*CoulH               
        end 
        end

    end    
 


 return  constq*(HartreeMatrix-FockMatrix)
end






function triangle_chern(Nq::Int,wave::Vector{Vector{Int}},scale::Float64,ϕ::Float64,flux::Float64,input_DensityMatrix::Vector{Matrix{ComplexF64}},constq::Float64)

  
    β=4*flux/(√3*scale^2)
    mass=0.5;
    dimension=length(wave)
    
    b1=scale*[0,1]
    b2=scale*[√3/2,-1/2]

    
    
    T1=b1/(Nq)
    T2=b2/(Nq)
    b1T=Int.(round.(inv([T1 T2])*b1))
    b2T=Int.(round.(inv([T1 T2])*b2))
    
    
    chern_allowedq=Vector{Int64}[]
    for ja in 0:Nq,jb in 0:Nq
        push!(chern_allowedq,[ja,jb])
    end
    
    chern_overlapmatrix=zeros(ComplexF64,(Nq+1)^2,length(wave),(Nq+1)^2,length(wave))
    for ja in 1:(Nq+1)^2, jb in eachindex(wave), jc in 1:(Nq+1)^2, jd in eachindex(wave)
       chern_overlapmatrix[ja,jb,jc,jd]=overlap([T1 T2]*(chern_allowedq[ja]+wave[jb]),[T1 T2]*(chern_allowedq[jc]+wave[jd]-wave[jb]-chern_allowedq[ja]),β)
    end
    loop_dic=construct_loop_dic(wave)
    
    
    eigenvector_intermediate_bc=Vector{Vector{ComplexF64}}(undef,(Nq+1)^2)
    eigenvector_intermediate_single=Vector{Vector{ComplexF64}}(undef,(Nq+1)^2)

    Threads.@threads for ja in eachindex(chern_allowedq)
        k=[T1 T2]*chern_allowedq[ja]
        chern_Ham=zeros(ComplexF64,dimension,dimension)
        chern_MoirePo=zeros(ComplexF64,dimension,dimension)
        for jb in eachindex(wave)
            chern_Ham[jb,jb]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2)^2/(2*mass)
        end
    
        for jc in eachindex(wave)
            k1=k+wave[jc][1]*T1+wave[jc][2]*T2
            pos=findfirst(item->item==wave[jc]-b1T,wave)
            if pos≠nothing
           chern_MoirePo[jc,pos]=V0*exp(im*ϕ)*overlap(k1,-b1,β)
            end
            
          
            pos=findfirst(item->item==wave[jc]+b2T+b1T,wave)
            if pos≠nothing
                chern_MoirePo[jc,pos]=V0*exp(im*ϕ)*overlap(k1,b2+b1,β)
            end
        
            pos=findfirst(item->item==wave[jc]-b2T,wave)
            if pos≠nothing
                chern_MoirePo[jc,pos]=V0*exp(im*ϕ)*overlap(k1,-b2,β)
            end
         end


      chern_MoirePo=chern_MoirePo+chern_MoirePo'
      HFmatrix=Construct_HFmatrix(loop_dic,ja,chern_allowedq[ja],allowedq,T1,T2,Nq,wave,input_DensityMatrix,constq,chern_overlapmatrix)
      eigenvector_intermediate_bc[ja]=eigvecs(chern_Ham+chern_MoirePo+HFmatrix)[:,1] 
      eigenvector_intermediate_single[ja]=eigvecs(chern_Ham+chern_MoirePo)[:,1] 
    
    end

    eigenvector_bc=zeros(ComplexF64,dimension,Nq+1,Nq+1)
    eigenvector_bc_single=zeros(ComplexF64,dimension,Nq+1,Nq+1)
    for ja in eachindex(chern_allowedq)
       eigenvector_bc[:,chern_allowedq[ja][1]+1,chern_allowedq[ja][2]+1]=eigenvector_intermediate_bc[ja]
       eigenvector_bc_single[:,chern_allowedq[ja][1]+1,chern_allowedq[ja][2]+1]=eigenvector_intermediate_single[ja]
    end
    
    
    Uonelink=zeros(ComplexF64,Nq,Nq+1)
    Utwolink=zeros(ComplexF64,Nq+1,Nq)
    tra=zeros(ComplexF64,Nq,Nq)
    tra_single=zeros(ComplexF64,Nq,Nq)
    
    for ja in 1:Nq, jb in 1:Nq+1
        Amatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T1,T1,T2)
       Uonelink[ja,jb]=dot(eigenvector_bc[:,ja,jb],Amatrix*eigenvector_bc[:,ja+1,jb])/abs(dot(eigenvector_bc[:,ja,jb],Amatrix*eigenvector_bc[:,ja+1,jb]))
    end

    
    
    for ja in 1:Nq+1, jb in 1:Nq
        Amatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T2,T1,T2)
     Utwolink[ja,jb]=dot(eigenvector_bc[:,ja,jb],Amatrix*eigenvector_bc[:,ja,jb+1])/abs(dot(eigenvector_bc[:,ja,jb],Amatrix*eigenvector_bc[:,ja,jb+1]))
    end
    
    dG=norm(T2)
    for ja in 1:Nq, jb in 1:Nq
        Amatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T1,T1,T2)
        Bmatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T2,T1,T2)
        Cmatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T2+T1,T1,T2)
        A1=dot(eigenvector_bc[:,ja,jb],Amatrix*eigenvector_bc[:,ja+1,jb])
        B1=dot(eigenvector_bc[:,ja,jb],Bmatrix*eigenvector_bc[:,ja,jb+1])
        C1=dot(eigenvector_bc[:,ja,jb],Cmatrix*eigenvector_bc[:,ja+1,jb+1])
       gyy=(1-abs(A1)^2)/dG^2
       gxx=(2-1/2*gyy*dG^2-abs(B1)^2-abs(C1)^2)/(1.5*dG^2)
       tra[ja,jb]+=(gxx+gyy)*3^(1/2)/2*dG^2

        A1=dot(eigenvector_bc_single[:,ja,jb],Amatrix*eigenvector_bc_single[:,ja+1,jb])
        B1=dot(eigenvector_bc_single[:,ja,jb],Bmatrix*eigenvector_bc_single[:,ja,jb+1])
        C1=dot(eigenvector_bc_single[:,ja,jb],Cmatrix*eigenvector_bc_single[:,ja+1,jb+1])
       gyy=(1-abs(A1)^2)/dG^2
       gxx=(2-1/2*gyy*dG^2-abs(B1)^2-abs(C1)^2)/(1.5*dG^2)
       tra_single[ja,jb]+=(gxx+gyy)*3^(1/2)/2*dG^2
      
      
    end

   


    Flink=zeros(ComplexF64,Nq,Nq)
    for ja in 1:Nq, jb in 1:Nq
     Flink[ja,jb]=log(Uonelink[ja,jb]*Utwolink[ja+1,jb]/(Uonelink[ja,jb+1]*Utwolink[ja,jb]))
    end
    chern=sum(Flink)/(2*π*im)
    aveF=sum(Flink)/Nq^2
    uniform=0
    for ja in 1:Nq, jb in 1:Nq
        uniform+=(imag(Flink[ja,jb])-imag(aveF))^2*Nq^2/(2π)^2
    end


    Uonelink=zeros(ComplexF64,Nq,Nq+1)
    Utwolink=zeros(ComplexF64,Nq+1,Nq)
  
    
    for ja in 1:Nq, jb in 1:Nq+1
        Amatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T1,T1,T2)
       Uonelink[ja,jb]=dot(eigenvector_bc_single[:,ja,jb],Amatrix*eigenvector_bc_single[:,ja+1,jb])/abs(dot(eigenvector_bc_single[:,ja,jb],Amatrix*eigenvector_bc_single[:,ja+1,jb]))
    end
    
    for ja in 1:Nq+1, jb in 1:Nq
        Amatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T2,T1,T2)
     Utwolink[ja,jb]=dot(eigenvector_bc_single[:,ja,jb],Amatrix*eigenvector_bc_single[:,ja,jb+1])/abs(dot(eigenvector_bc_single[:,ja,jb],Amatrix*eigenvector_bc_single[:,ja,jb+1]))
    end
    

   

    Flink_single=zeros(ComplexF64,Nq,Nq)
    for ja in 1:Nq, jb in 1:Nq
     Flink_single[ja,jb]=log(Uonelink[ja,jb]*Utwolink[ja+1,jb]/(Uonelink[ja,jb+1]*Utwolink[ja,jb]))
    end
    chern_single=sum(Flink_single)/(2*π*im)
    aveF=sum(Flink_single)/Nq^2
    uniform_single=0
    for ja in 1:Nq, jb in 1:Nq
        uniform_single+=(imag(Flink_single[ja,jb])-imag(aveF))^2*Nq^2/(2π)^2
    end
    
    

    
  trace_condition=sum(tra)-sum(abs.(Flink))
  trace_condition_single=sum(tra_single)-sum(abs.(Flink_single))
      


    
      
   return chern,Flink,chern_single,Flink_single,trace_condition,trace_condition_single,uniform,uniform_single

end







function calculate_energy(Nq::Int,wave::Vector{Vector{Int}},scale::Float64,ϕ::Float64,flux::Float64,input_DensityMatrix::Vector{Matrix{ComplexF64}},constq::Float64,overlapmatrix::Array{ComplexF64,4})

  
    β=4*flux/(√3*scale^2)
    mass=0.5;
    dimension=length(wave)
    
    b1=scale*[0,1]
    b2=scale*[√3/2,-1/2]

    
    
    T1=b1/(Nq)
    T2=b2/(Nq)
    b1T=Int.(round.(inv([T1 T2])*b1))
    b2T=Int.(round.(inv([T1 T2])*b2))
    
    
    chern_allowedq=Vector{Int64}[]
    for ja in 0:Nq-1,jb in 0:Nq-1
        push!(chern_allowedq,[ja,jb])
    end
    

    loop_dic=construct_loop_dic(wave)
    
    Energy_Matrix=[zeros(ComplexF64,length(wave),length(wave)) for _ in 1:Nq^2]
    
    Threads.@threads for ja in eachindex(chern_allowedq)
        k=[T1 T2]*chern_allowedq[ja]
        chern_Ham=zeros(ComplexF64,dimension,dimension)
        chern_MoirePo=zeros(ComplexF64,dimension,dimension)
        for jb in eachindex(wave)
            chern_Ham[jb,jb]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2)^2/(2*mass)
        end
    
        for jc in eachindex(wave)
            k1=k+wave[jc][1]*T1+wave[jc][2]*T2
            pos=findfirst(item->item==wave[jc]-b1T,wave)
            if pos≠nothing
           chern_MoirePo[jc,pos]=V0*exp(im*ϕ)*overlap(k1,-b1,β)
            end
            
          
            pos=findfirst(item->item==wave[jc]+b2T+b1T,wave)
            if pos≠nothing
                chern_MoirePo[jc,pos]=V0*exp(im*ϕ)*overlap(k1,b2+b1,β)
            end
        
            pos=findfirst(item->item==wave[jc]-b2T,wave)
            if pos≠nothing
                chern_MoirePo[jc,pos]=V0*exp(im*ϕ)*overlap(k1,-b2,β)
            end
         end


      chern_MoirePo=chern_MoirePo+chern_MoirePo'
      HFmatrix=Construct_HFmatrix(loop_dic,ja,chern_allowedq[ja],allowedq,T1,T2,Nq,wave,input_DensityMatrix,constq,overlapmatrix)
      Energy_Matrix[ja]=1/2*HFmatrix+chern_MoirePo+chern_Ham
    
    end
    
    
   energy=0
   for ja in 1:Nq^2
       energy+=tr(Energy_Matrix[ja]*input_DensityMatrix[ja])
   end


   return energy

end