using LinearAlgebra
using Arpack
using Combinatorics

using Random



function overlap(k::Vector{Float64},q::Vector{Float64},β::Float64)::ComplexF64
    v=q[1]^2+q[2]^2+2*im*(k[1]*q[2]-k[2]*q[1])
    #v=2*im*(k[1]*q[2]-k[2]*q[1])
    return exp(-β/4*v)
 end
 
 
 
 
function Coulomb(k::Vector{Int64},T1::Vector{Float64},T2::Vector{Float64})::Float64
  
 
   return k==[0,0] ? 0.0 : 1/norm(k[1]*T1+k[2]*T2)
   
end







function triangle_initial_Densitymatrix(flux::Float64,V0::Float64,ϕ::Float64,scale::Float64,Nq::Int64)
    am=4*π/(√3*scale);
    β=4*flux/(√3*scale^2)
    mass=0.5;
    
    b1=4*π/(√3*am)*[0,1]
    b2=4*π/(√3*am)*[√3/2,-1/2]


    a1m=am*[1/2,√3/2]
    a2m=am*[1,0]
    
    
    T1=b1/(Nq)
    T2=b2/(Nq)

    Vseed=0.2;
    ϕseed=π/3
    M_ϕseed=π
    
    
    
    b1T=Int.(round.(inv([T1 T2])*b1))
    b2T=Int.(round.(inv([T1 T2])*b2))
    
    
    
    allowedq=Vector{Int64}[]
    for ja in 0:Nq-1,jb in 0:Nq-1
        push!(allowedq,[ja,jb])
    end

    
    
    
    wave=Vector{Int64}[]
    cutoff=18
    cutoffstandard=4.01*scale
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        gtest=ja*b1+jb*b2;
        if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
            push!(wave,ja*b1T+jb*b2T)
        end
    end
    dimension=length(wave)
   
    

    overlapmatrix=[zeros(ComplexF64,Nq^2,length(wave),Nq^2,length(wave)) for _ in 1:2]
    for ja in 1:Nq^2, jb in eachindex(wave), jc in 1:Nq^2, jd in eachindex(wave)
       overlapmatrix[1][ja,jb,jc,jd]=overlap([T1 T2]*(allowedq[ja]+wave[jb]),[T1 T2]*(allowedq[jc]+wave[jd]-wave[jb]-allowedq[ja]),β)
    end
    
    overlapmatrix[2]=conj(overlapmatrix[1])
    
 
    
    single_Ham=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
    single_MoirePo=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
    single_eigenvalue=[[zeros(Float64,dimension) for _ in 1:2] for _ in 1:Nq^2]
    single_eigenvector=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
    seed_MoirePo=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
    seed_eigenvector=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]

    Threads.@threads for ja in 1:Nq^2
      
      
        k=allowedq[ja][1]*T1+allowedq[ja][2]*T2

        for jb in eachindex(wave), vi in 1:2
         single_Ham[ja][vi][jb,jb]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2)^2/(2*mass)
        end
        
       for jc in eachindex(wave)
          k1=k+wave[jc][1]*T1+wave[jc][2]*T2
          pos=findfirst(item->item==wave[jc]-b1T,wave)
          if pos≠nothing
            single_MoirePo[ja][1][jc,pos]=V0*exp(im*ϕ)*overlap(k1,-b1,β)
            single_MoirePo[ja][2][jc,pos]=V0*exp(im*ϕ)*conj(overlap(k1,-b1,β))
            seed_MoirePo[ja][1][jc,pos]=Vseed*exp(im*ϕseed)*overlap(k1,-b1,β)
            seed_MoirePo[ja][2][jc,pos]=Vseed*exp(im*M_ϕseed)*conj(overlap(k1,-b1,β))
          end
          
        
          pos=findfirst(item->item==wave[jc]+b2T+b1T,wave)
          if pos≠nothing
              single_MoirePo[ja][1][jc,pos]=V0*exp(im*ϕ)*overlap(k1,b2+b1,β)
              single_MoirePo[ja][2][jc,pos]=V0*exp(im*ϕ)*conj(overlap(k1,b2+b1,β))
              seed_MoirePo[ja][1][jc,pos]=Vseed*exp(im*ϕseed)*overlap(k1,b2+b1,β)
              seed_MoirePo[ja][2][jc,pos]=Vseed*exp(im*M_ϕseed)*conj(overlap(k1,b2+b1,β))
          end
      
          pos=findfirst(item->item==wave[jc]-b2T,wave)
          if pos≠nothing
              single_MoirePo[ja][1][jc,pos]=V0*exp(im*ϕ)*overlap(k1,-b2,β)
              single_MoirePo[ja][2][jc,pos]=V0*exp(im*ϕ)*conj(overlap(k1,-b2,β))
              seed_MoirePo[ja][1][jc,pos]=Vseed*exp(im*ϕseed)*overlap(k1,-b2,β)
              seed_MoirePo[ja][2][jc,pos]=Vseed*exp(im*M_ϕseed)*conj(overlap(k1,-b2,β))
          end
       end
      
       for vi in 1:2
       single_MoirePo[ja][vi]=single_MoirePo[ja][vi]+single_MoirePo[ja][vi]'
       seed_MoirePo[ja][vi]=seed_MoirePo[ja][vi]+seed_MoirePo[ja][vi]'

       FFF=eigen(single_MoirePo[ja][vi]+single_Ham[ja][vi])
      
        single_eigenvalue[ja][vi]=real(FFF.values)
        single_eigenvector[ja][vi]=FFF.vectors

        FFF=eigen(single_MoirePo[ja][vi]+single_Ham[ja][vi]+seed_MoirePo[ja][vi])
      
       
        seed_eigenvector[ja][vi]=FFF.vectors

       end
    
    end


    input_DensityMatrix=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
    #for ja in 1:Nq^2, vi in 1:2
     
       #input_DensityMatrix[ja][vi]+=(single_eigenvector[ja][vi][:,1]*(single_eigenvector[ja][vi][:,1])')
        
    #end

    for ja in 1:Nq^2, vi in 1:2
     
       input_DensityMatrix[ja][vi]+=(seed_eigenvector[ja][vi][:,1]*(seed_eigenvector[ja][vi][:,1])')
        
    end

    

    
    return overlapmatrix, wave, input_DensityMatrix, single_MoirePo, single_Ham, seed_MoirePo,single_eigenvalue,allowedq, T1, T2, a1m, a2m
      

       
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
        
        loop_dic[wave[g2index]-wave[gindex]][[gindex,g2index]]=Vector{Int64}[]
        
         
    end
    
    for gindex in eachindex(wave), g1index in eachindex(wave),g2index in eachindex(wave)
       
        if haskey(g_dic,wave[gindex]+wave[g1index]-wave[g2index])
            push!(loop_dic[wave[g2index]-wave[gindex]][[gindex,g2index]],[g1index,g_dic[wave[gindex]+wave[g1index]-wave[g2index]]])
        end
        
    end
    
   return loop_dic

end



function Construct_DensityMatrix(loop_dic::Dict{Vector{Int},Any},allowedq::Vector{Vector{Int}},T1::Vector{Float64},T2::Vector{Float64},Nq::Int64,wave::Vector{Vector{Int64}},input_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},single_Ham::Vector{Vector{Matrix{ComplexF64}}},single_MoirePo::Vector{Vector{Matrix{ComplexF64}}},constq::Float64, ζ::Float64,overlapmatrix::Vector{Array{ComplexF64,4}})::Tuple{Float64,Vector{Vector{Matrix{ComplexF64}}},Vector{Vector{Matrix{ComplexF64}}},Vector{Vector{Vector{Float64}}}}
  
 
    dimension=length(wave)
   HartreeMatrix=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
   FockMatrix=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
    output_DensityMatrix=Vector{Any}(undef,Nq^2)
   DeltaMatrix=Vector{Any}(undef,Nq^2)
   NewDensityMatrix=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
   HF_eigenvalue=[[zeros(Float64,dimension) for _ in 1:2] for _ in 1:Nq^2]
   HF_eigenvector=[[zeros(ComplexF64,dimension,dimension) for _ in 1:2] for _ in 1:Nq^2]
 
  
   Threads.@threads for jk in 1:Nq^2
    for vi in 1:2
     Fk = FockMatrix[jk][vi]
     for jk1 in 1:Nq^2
         dmk = input_DensityMatrix[jk1][vi]
         q=allowedq[jk1]-allowedq[jk]
        for (dg,loop_dic_dg) in loop_dic
            CoulF1=Coulomb(q+dg,T1,T2)     
            for (gg2,loop_dic_dg_gg2) in loop_dic_dg          
                CoulF=CoulF1*overlapmatrix[vi][jk1,gg2[2],jk,gg2[1]]
            for g1g3 in loop_dic_dg_gg2
                Fk[g1g3[2],gg2[1]]+=dmk[g1g3[1],gg2[2]]*CoulF*overlapmatrix[vi][jk,g1g3[2],jk1,g1g3[1]]
            end 
            end
    
        end    
     end
   end
   end
  
   Hartree_Density=[zeros(ComplexF64,dimension,dimension) for _ in 1:2]
   for ja in 1:Nq^2, vi in 1:2
    Hartree_Density[vi]+=input_DensityMatrix[ja][vi] .* transpose(overlapmatrix[vi][ja,:,ja,:])
   end
 
   
 
  Threads.@threads for jk in 1:Nq^2
   for vi in 1:2
        if vi==1
         ovi=2
        else 
         ovi=1
        end
 
       for dg in keys(loop_dic)
           CoulH1=Coulomb(dg,T1,T2)
           for gg2 in keys(loop_dic[dg])          
               CoulH=CoulH1*overlapmatrix[vi][jk,gg2[2],jk,gg2[1]]            
           for g1g3 in loop_dic[dg][gg2]           
               HartreeMatrix[jk][vi][gg2[2],gg2[1]]+=Hartree_Density[vi][g1g3[1],g1g3[2]]*CoulH+ζ*Hartree_Density[ovi][g1g3[1],g1g3[2]]*CoulH                  
           end 
           end
   
       end    
  end
  end
 
 
  
 
  Threads.@threads for ja in 1:Nq^2
    for vi in 1:2
    FFF=eigen(single_MoirePo[ja][vi]+single_Ham[ja][vi]+constq*HartreeMatrix[ja][vi]-constq*FockMatrix[ja][vi])
    HF_eigenvalue[ja][vi]=real(FFF.values)
    HF_eigenvector[ja][vi]=FFF.vectors
    end
  end
 
  bound=sort(reduce(vcat,reduce(vcat,HF_eigenvalue)))[2*Nq^2+1]
 
   for ja in 1:Nq^2, vi in 1:2
     
        for jd in eachindex(HF_eigenvalue[ja][vi])
           if HF_eigenvalue[ja][vi][jd]<bound
              NewDensityMatrix[ja][vi]+=HF_eigenvector[ja][vi][:,jd]*(HF_eigenvector[ja][vi][:,jd])'
           end
        end
       
   end
   DeltaMatrix=NewDensityMatrix-input_DensityMatrix
   output_DensityMatrix=0.0*input_DensityMatrix+1.0*NewDensityMatrix
 
   
   e1=0.0
   for ja in 1:Nq^2, vi in 1:2
     e1+=tr(DeltaMatrix[ja][vi]'*DeltaMatrix[ja][vi])
   end
   eout=real(e1)/(2*Nq^2)
   
   
  
 
  return  eout,output_DensityMatrix,DeltaMatrix,HF_eigenvalue
end
 




function iteration_loop(initial_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},allowedq::Vector{Vector{Int64}},T1::Vector{Float64},T2::Vector{Float64},Nq::Int,wave::Vector{Vector{Int64}},single_Ham::Vector{Vector{Matrix{ComplexF64}}},single_MoirePo::Vector{Vector{Matrix{ComplexF64}}},seed_MoirePo::Vector{Vector{Matrix{ComplexF64}}},constq::Float64,ζ::Float64,overlapmatrix::Vector{Array{ComplexF64,4}})::Tuple{Vector{Vector{Vector{Matrix{ComplexF64}}}},Vector{Vector{Vector{Matrix{ComplexF64}}}},Vector{Vector{Vector{Float64}}}}
    eout=1.0
    itcount=0
    loop_dic=construct_loop_dic(wave)
    HF_eigenvalue=Vector{Any}(undef,Nq^2)
    DIIS_input_DensityMatrix=Vector{Vector{Vector{Matrix{ComplexF64}}}}(undef,3)
    DIIS_input_DeltaMatrix=Vector{Vector{Vector{Matrix{ComplexF64}}}}(undef,3)
    input_DensityMatrix=initial_DensityMatrix
    bad_count=0
     
    while itcount<10
       
        tic=time()
        eout,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue=Construct_DensityMatrix(loop_dic,allowedq,T1,T2,Nq,wave,input_DensityMatrix,single_Ham+seed_MoirePo,single_MoirePo,constq,ζ,overlapmatrix)
        DIIS_input_DensityMatrix[mod(itcount,3)+1]=input_DensityMatrix
        input_DensityMatrix=output_DensityMatrix
        itcount+=1
        toc=time()
        println(toc-tic,"eout=$eout")
        flush(stdout)
       
      
    end
    println("takeaway the fake potential")
    flush(stdout)

    while (eout>1*10^-13) || (bad_count<4)
      if  eout<1*10^-13 
        bad_count+=1
      end
      tic=time()
      eout,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue=Construct_DensityMatrix(loop_dic,allowedq,T1,T2,Nq,wave,input_DensityMatrix,single_Ham,single_MoirePo,constq,ζ,overlapmatrix)
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
            for jc in 1:Nq^2, vi in 1:2
               Bmatrix[ja,jb]+=tr((DIIS_input_DeltaMatrix[ja][vi][jc])'*(DIIS_input_DeltaMatrix[jb][vi][jc]))
            end
        end
        coeff=inv(Bmatrix)*[0;0;0;1]
       
        dmk=coeff[1]*(DIIS_input_DensityMatrix[1]+DIIS_input_DeltaMatrix[1])+coeff[2]*(DIIS_input_DensityMatrix[2]+DIIS_input_DeltaMatrix[2])+coeff[3]*(DIIS_input_DensityMatrix[3]+DIIS_input_DeltaMatrix[3])
         eout,_,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue=Construct_DensityMatrix(loop_dic,allowedq,T1,T2,Nq,wave,dmk,single_Ham,single_MoirePo,constq,ζ,overlapmatrix)
        DIIS_input_DensityMatrix[mod(itcount,3)+1]=dmk
        itcount+=1

        toc=time()
        println(toc-tic,"eout=$eout")
        flush(stdout)
    end

 
  =#

  return DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue

end


function Densitymap(a1m::Vector{Float64},a2m::Vector{Float64},overlapmatrix::Vector{Array{ComplexF64,4}},wave::Vector{Vector{Int}},input_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}})
    N3=50
    dimension=length(wave)
    zgrid=[zeros(Float64,N3,N3) for _ in 1:2]
    Hartree_Density=[zeros(ComplexF64,dimension,dimension) for _ in 1:2]


    for jk1 in 1:Nq^2, vi in 1:2
        Hartree_Density[vi]+=input_DensityMatrix[jk1][vi] .* transpose(overlapmatrix[vi][jk1,:,jk1,:])
    end
    
    for ja in 1:50, jb in 1:50
      rvec=ja/25*a1m+jb/25*a2m
      for jc in eachindex(wave), jd in eachindex(wave), vi in 1:2
        gvec=[T1 T2]*(wave[jc]-wave[jd])
       zgrid[vi][ja,jb]+=real(Hartree_Density[vi][jc,jd]*exp(im*(gvec[1]*rvec[1]+gvec[2]*rvec[2])))
      end
    
    end
    
    return zgrid
end






function metric(wavelist::Vector{Vector{Int64}},β::Float64,k::Vector{Float64},q::Vector{Float64},T1::Vector{Float64},T2::Vector{Float64},vi::Int)::Matrix{ComplexF64}
    Amatrix=zeros(ComplexF64,length(wavelist),length(wavelist))
    if vi==1
     for ja in 1:length(wavelist)
     Amatrix[ja,ja]=overlap(k+wavelist[ja][1]*T1+wavelist[ja][2]*T2,q,β)
     end
    else
        for ja in 1:length(wavelist)
            Amatrix[ja,ja]=conj(overlap(k+wavelist[ja][1]*T1+wavelist[ja][2]*T2,q,β))
        end
    end
    return Amatrix
end




function Construct_HFmatrix(loop_dic::Dict{Vector{Int},Any},vi::Int64,pathpointindex::Int64,pathpoint::Vector{Int64},allowedq::Vector{Vector{Int}},T1::Vector{Float64},T2::Vector{Float64},Nq::Int64,wave::Vector{Vector{Int64}},input_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},constq::Float64,ζ::Float64,chern_overlapmatrix::Vector{Array{ComplexF64,4}})::Matrix{ComplexF64}
  
  
    dimension=length(wave)
    HartreeMatrix=zeros(ComplexF64,dimension,dimension) 
    FockMatrix=zeros(ComplexF64,dimension,dimension) 
   
  
  
  
      
      for jk1 in 1:Nq^2
          dmk = input_DensityMatrix[jk1][vi]
          q=allowedq[jk1]-pathpoint
         for (dg,loop_dic_dg) in loop_dic
             CoulF1=Coulomb(q+dg,T1,T2)     
             for (gg2,loop_dic_dg_gg2) in loop_dic_dg          
                 CoulF=CoulF1*chern_overlapmatrix[vi][jk1,gg2[2],pathpointindex,gg2[1]]
             for g1g3 in loop_dic_dg_gg2
                 FockMatrix[g1g3[2],gg2[1]]+=dmk[g1g3[1],gg2[2]]*CoulF*chern_overlapmatrix[vi][pathpointindex,g1g3[2],jk1,g1g3[1]]
             end 
             end
     
         end    
      end
  
    
  
    Hartree_Density=[zeros(ComplexF64,dimension,dimension) for _ in 1:2]
    for jk1 in 1:Nq^2, vi in 1:2
     Hartree_Density[vi]+=input_DensityMatrix[jk1][vi] .* transpose(chern_overlapmatrix[vi][jk1,:,jk1,:])
    end
  
    if vi==1
      ovi=2
    else
      ovi=1
    end
   
  
      
   
      for dg in keys(loop_dic)
          CoulH1=Coulomb(dg,T1,T2)
          for gg2 in keys(loop_dic[dg])          
              CoulH=CoulH1*chern_overlapmatrix[vi][pathpointindex,gg2[2],pathpointindex,gg2[1]]            
          for g1g3 in loop_dic[dg][gg2]           
              HartreeMatrix[gg2[2],gg2[1]]+=Hartree_Density[vi][g1g3[1],g1g3[2]]*CoulH+ζ*Hartree_Density[ovi][g1g3[1],g1g3[2]]*CoulH             
          end 
          end
  
      end    
   
  
  
   
  
   return  constq*(HartreeMatrix-FockMatrix)
end






function triangle_chern(Nq::Int,wave::Vector{Vector{Int}},scale::Float64,ϕ::Float64,flux::Float64,input_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},constq::Float64,ζ::Float64)

  
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
    
    chern_overlapmatrix=[zeros(ComplexF64,(Nq+1)^2,length(wave),(Nq+1)^2,length(wave)) for _ in 1:2]
    for ja in 1:(Nq+1)^2, jb in eachindex(wave), jc in 1:(Nq+1)^2, jd in eachindex(wave)
       chern_overlapmatrix[1][ja,jb,jc,jd]=overlap([T1 T2]*(chern_allowedq[ja]+wave[jb]),[T1 T2]*(chern_allowedq[jc]+wave[jd]-wave[jb]-chern_allowedq[ja]),β)
    end
    chern_overlapmatrix[2]=conj(chern_overlapmatrix[1])

    loop_dic=construct_loop_dic(wave)
    
    
    eigenvector_intermediate_bc=[[zeros(ComplexF64,dimension) for _ in 1:2] for _ in 1:(Nq+1)^2]
    eigenvector_intermediate_single=[[zeros(ComplexF64,dimension) for _ in 1:2] for _ in 1:(Nq+1)^2]

   for ja in eachindex(chern_allowedq)
        for vi in 1:2
        k=[T1 T2]*chern_allowedq[ja]
        chern_Ham=zeros(ComplexF64,dimension,dimension)
        chern_MoirePo=zeros(ComplexF64,dimension,dimension)
        for jb in eachindex(wave)
            chern_Ham[jb,jb]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2)^2/(2*mass)
        end
        if vi==1
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

        else
            for jc in eachindex(wave)
                k1=k+wave[jc][1]*T1+wave[jc][2]*T2
                pos=findfirst(item->item==wave[jc]-b1T,wave)
                if pos≠nothing
               chern_MoirePo[jc,pos]=V0*exp(im*ϕ)*conj(overlap(k1,-b1,β))
                end
                
              
                pos=findfirst(item->item==wave[jc]+b2T+b1T,wave)
                if pos≠nothing
                    chern_MoirePo[jc,pos]=V0*exp(im*ϕ)*conj(overlap(k1,b2+b1,β))
                end
            
                pos=findfirst(item->item==wave[jc]-b2T,wave)
                if pos≠nothing
                    chern_MoirePo[jc,pos]=V0*exp(im*ϕ)*conj(overlap(k1,-b2,β))
                end
             end
        end


      chern_MoirePo=chern_MoirePo+chern_MoirePo'
      HFmatrix=Construct_HFmatrix(loop_dic,vi,ja,chern_allowedq[ja],allowedq,T1,T2,Nq,wave,input_DensityMatrix,constq,ζ,chern_overlapmatrix)
      eigenvector_intermediate_bc[ja][vi]=eigvecs(chern_Ham+chern_MoirePo+HFmatrix)[:,1] 
      eigenvector_intermediate_single[ja][vi]=eigvecs(chern_Ham+chern_MoirePo)[:,1] 
    
     end
  end

    eigenvector_bc=[zeros(ComplexF64,dimension,Nq+1,Nq+1) for _ in 1:2]
    eigenvector_bc_single=[zeros(ComplexF64,dimension,Nq+1,Nq+1) for _ in 1:2]
    for ja in eachindex(chern_allowedq), vi in 1:2
       eigenvector_bc[vi][:,chern_allowedq[ja][1]+1,chern_allowedq[ja][2]+1]=eigenvector_intermediate_bc[ja][vi]
       eigenvector_bc_single[vi][:,chern_allowedq[ja][1]+1,chern_allowedq[ja][2]+1]=eigenvector_intermediate_single[ja][vi]
    end
    
    
    Uonelink=[zeros(ComplexF64,Nq,Nq+1) for _ in 1:2]
    Utwolink=[zeros(ComplexF64,Nq+1,Nq) for _ in 1:2]

    Uonelink_single=[zeros(ComplexF64,Nq,Nq+1) for _ in 1:2]
    Utwolink_single=[zeros(ComplexF64,Nq+1,Nq) for _ in 1:2]
    tra=[zeros(ComplexF64,Nq,Nq) for _ in 1:2]
    tra_single=[zeros(ComplexF64,Nq,Nq) for _ in 1:2]
    Flink=[zeros(ComplexF64,Nq,Nq) for _ in 1:2]
    Flink_single=[zeros(ComplexF64,Nq,Nq) for _ in 1:2]
 for vi in 1:2
    
    for ja in 1:Nq, jb in 1:Nq+1
        Amatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T1,T1,T2,vi)
       Uonelink[vi][ja,jb]=dot(eigenvector_bc[vi][:,ja,jb],Amatrix*eigenvector_bc[vi][:,ja+1,jb])/abs(dot(eigenvector_bc[vi][:,ja,jb],Amatrix*eigenvector_bc[vi][:,ja+1,jb]))
    end

    
    
    for ja in 1:Nq+1, jb in 1:Nq
        Amatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T2,T1,T2,vi)
     Utwolink[vi][ja,jb]=dot(eigenvector_bc[vi][:,ja,jb],Amatrix*eigenvector_bc[vi][:,ja,jb+1])/abs(dot(eigenvector_bc[vi][:,ja,jb],Amatrix*eigenvector_bc[vi][:,ja,jb+1]))
    end
    
    dG=norm(T2)
    for ja in 1:Nq, jb in 1:Nq
        Amatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T1,T1,T2,vi)
        Bmatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T2,T1,T2,vi)
        Cmatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T2+T1,T1,T2,vi)
        A1=dot(eigenvector_bc[vi][:,ja,jb],Amatrix*eigenvector_bc[vi][:,ja+1,jb])
        B1=dot(eigenvector_bc[vi][:,ja,jb],Bmatrix*eigenvector_bc[vi][:,ja,jb+1])
        C1=dot(eigenvector_bc[vi][:,ja,jb],Cmatrix*eigenvector_bc[vi][:,ja+1,jb+1])
       gyy=(1-abs(A1)^2)/dG^2
       gxx=(2-1/2*gyy*dG^2-abs(B1)^2-abs(C1)^2)/(1.5*dG^2)
       tra[vi][ja,jb]+=(gxx+gyy)*3^(1/2)/2*dG^2

        A1=dot(eigenvector_bc_single[vi][:,ja,jb],Amatrix*eigenvector_bc_single[vi][:,ja+1,jb])
        B1=dot(eigenvector_bc_single[vi][:,ja,jb],Bmatrix*eigenvector_bc_single[vi][:,ja,jb+1])
        C1=dot(eigenvector_bc_single[vi][:,ja,jb],Cmatrix*eigenvector_bc_single[vi][:,ja+1,jb+1])
       gyy=(1-abs(A1)^2)/dG^2
       gxx=(2-1/2*gyy*dG^2-abs(B1)^2-abs(C1)^2)/(1.5*dG^2)
       tra_single[vi][ja,jb]+=(gxx+gyy)*3^(1/2)/2*dG^2
      
      
    end

    for ja in 1:Nq, jb in 1:Nq
        Flink[vi][ja,jb]=log(Uonelink[vi][ja,jb]*Utwolink[vi][ja+1,jb]/(Uonelink[vi][ja,jb+1]*Utwolink[vi][ja,jb]))
    end

  
    
    for ja in 1:Nq, jb in 1:Nq+1
        Amatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T1,T1,T2,vi)
       Uonelink_single[vi][ja,jb]=dot(eigenvector_bc_single[vi][:,ja,jb],Amatrix*eigenvector_bc_single[vi][:,ja+1,jb])/abs(dot(eigenvector_bc_single[vi][:,ja,jb],Amatrix*eigenvector_bc_single[vi][:,ja+1,jb]))
    end
    
    for ja in 1:Nq+1, jb in 1:Nq
        Amatrix=metric(wave,β,(ja-1)*T1+(jb-1)*T2,T2,T1,T2,vi)
     Utwolink_single[vi][ja,jb]=dot(eigenvector_bc_single[vi][:,ja,jb],Amatrix*eigenvector_bc_single[vi][:,ja,jb+1])/abs(dot(eigenvector_bc_single[vi][:,ja,jb],Amatrix*eigenvector_bc_single[vi][:,ja,jb+1]))
    end
  
    for ja in 1:Nq, jb in 1:Nq
     Flink_single[vi][ja,jb]=log(Uonelink[vi][ja,jb]*Utwolink[vi][ja+1,jb]/(Uonelink[vi][ja,jb+1]*Utwolink[vi][ja,jb]))
    end
 end

   
   chern=zeros(ComplexF64,2)
   uniform=zeros(ComplexF64,2)
   chern_single=zeros(ComplexF64,2)
   uniform_single=zeros(ComplexF64,2)
   trace_condition_single=zeros(ComplexF64,2)
   trace_condition=zeros(ComplexF64,2)
   
for vi in 1:2
    chern[vi]=sum(Flink[vi])/(2*π*im)
    aveF=sum(Flink[vi])/Nq^2
  
    for ja in 1:Nq, jb in 1:Nq
        uniform[vi]+=(imag(Flink[vi][ja,jb])-imag(aveF))^2*Nq^2/(2π)^2
    end


    chern_single[vi]=sum(Flink_single[vi])/(2*π*im)
    aveF=sum(Flink_single[vi])/Nq^2
    
    for ja in 1:Nq, jb in 1:Nq
        uniform_single[vi]+=(imag(Flink_single[vi][ja,jb])-imag(aveF))^2*Nq^2/(2π)^2
    end

    

    
  trace_condition[vi]=sum(tra[vi])-sum(abs.(Flink[vi]))
  trace_condition_single[vi]=sum(tra_single[vi])-sum(abs.(Flink_single[vi]))
      
 end

    
      
   return chern,Flink,chern_single,Flink_single,trace_condition,trace_condition_single,uniform,uniform_single

end







function calculate_energy(Nq::Int,wave::Vector{Vector{Int}},scale::Float64,ϕ::Float64,flux::Float64,input_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},constq::Float64,ζ::Float64,overlapmatrix::Vector{Array{ComplexF64,4}})

  
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
    
    Energy_Matrix=[[zeros(ComplexF64,length(wave),length(wave)) for _ in 1:2] for _ in 1:Nq^2]
    
    Threads.@threads for ja in eachindex(chern_allowedq)
        for vi in 1:2
        k=[T1 T2]*chern_allowedq[ja]
        chern_Ham=zeros(ComplexF64,dimension,dimension)
        chern_MoirePo=zeros(ComplexF64,dimension,dimension)
        for jb in eachindex(wave)
            chern_Ham[jb,jb]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2)^2/(2*mass)
        end
      if vi==1
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
      else
        for jc in eachindex(wave)
            k1=k+wave[jc][1]*T1+wave[jc][2]*T2
            pos=findfirst(item->item==wave[jc]-b1T,wave)
            if pos≠nothing
           chern_MoirePo[jc,pos]=V0*exp(im*ϕ)*conj(overlap(k1,-b1,β))
            end
            
          
            pos=findfirst(item->item==wave[jc]+b2T+b1T,wave)
            if pos≠nothing
                chern_MoirePo[jc,pos]=V0*exp(im*ϕ)*conj(overlap(k1,b2+b1,β))
            end
        
            pos=findfirst(item->item==wave[jc]-b2T,wave)
            if pos≠nothing
                chern_MoirePo[jc,pos]=V0*exp(im*ϕ)*conj(overlap(k1,-b2,β))
            end
         end
      end


      chern_MoirePo=chern_MoirePo+chern_MoirePo'
      HFmatrix=Construct_HFmatrix(loop_dic,vi,ja,chern_allowedq[ja],allowedq,T1,T2,Nq,wave,input_DensityMatrix,constq,ζ,overlapmatrix)
      Energy_Matrix[ja][vi]=1/2*HFmatrix+chern_MoirePo+chern_Ham
        end
    end
    
    
   energy=0
   for ja in 1:Nq^2, vi in 1:2
       energy+=tr(Energy_Matrix[ja][vi]*input_DensityMatrix[ja][vi])
   end


   return energy

end