using LinearAlgebra
using Arpack
using Combinatorics
using Random



 

function overlap(k::Vector{Float64},q::Vector{Float64},spin::Float64,M::Float64)::ComplexF64
    #v=(M^2+norm(k)^2+k[1]*q[1]+k[2]*q[2]-im*(k[1]*q[2]-k[2]*q[1]))^(Int(2*spin))/((M^2+norm(k)^2)^spin*(M^2+norm(k+q)^2)^spin)
    f1=-2*M^2*spin*(norm(q)^2)+(norm(k+q)^2+M^2)*(norm(k)^2+M^2)
    f2=(norm(k)^2+M^2+k[1]*q[1]+k[2]*q[2]-im*(k[1]*q[2]-k[2]*q[1]))^(Int(2*spin)-2)
    f3=(norm(k)^2+M^2)^spin
    f4=(norm(k+q)^2+M^2)^spin
    return f1*f2/(f3*f4)
end
 
 
function Coulomb(k::Vector{Int64},T1::Vector{Float64},T2::Vector{Float64})::Float64
  
 
   return k==[0,0] ? 0.0 : 1/norm(k[1]*T1+k[2]*T2)
   
end








function triangle_initial_Densitymatrix(spin::Float64,vf::Float64,V0::Float64,ϕ::Float64,scale::Float64,Geonum::Int64)
    am=4*π/(√3*scale);
    mass=0.5;
    
    Nx,Ny,l1,l2=Geometry(Geonum)
    b1=4*π/(√3*am)*[0,1]
    b2=4*π/(√3*am)*[√3/2,-1/2]

   
    a2m=am*[1/2,√3/2]
    a1m=am*[1,0]


    L1=l1[1]*a1m+l1[2]*a2m;
    L2=l2[1]*a1m+l2[2]*a2m;
    area=abs(L1[1]*L2[2]-L2[1]*L1[2]);
    Rotminus90=[0 1;-1 0]
    T1=2*π/area*Rotminus90*L2
    T2=-2*π/area*Rotminus90*L1
    
   
    
    
    
    b1T=Int.(round.(inv([T1 T2])*b1))
    b2T=Int.(round.(inv([T1 T2])*b2))
    
    
    
    allowedq=Vector{Int64}[]
    for ja in 0:Nx-1,jb in 0:Ny-1
        push!(allowedq,[ja,jb])
    end

     
   for ja in 1:Nx*Ny
    allowedq[ja]=sendtomesh([b1T';b2T'],allowedq[ja])
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
   
    
    
    overlapmatrix=zeros(ComplexF64,length(allowedq),length(wave),length(allowedq),length(wave))
    for ja in eachindex(allowedq), jb in eachindex(wave), jc in eachindex(allowedq), jd in eachindex(wave)
       overlapmatrix[ja,jb,jc,jd]=overlap([T1 T2]*(allowedq[ja]+wave[jb]),[T1 T2]*(allowedq[jc]+wave[jd]-wave[jb]-allowedq[ja]),spin,mass*vf/2)
    end


    
     
    single_Ham=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nx*Ny]
    single_MoirePo=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nx*Ny]
    single_eigenvalue=[zeros(Float64,dimension) for _ in 1:Nx*Ny]
    single_eigenvector=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nx*Ny]
    
    
    Threads.@threads for ja in eachindex(allowedq)
      
      
      k=allowedq[ja][1]*T1+allowedq[ja][2]*T2
    
      for jb in eachindex(wave)
       single_Ham[ja][jb,jb]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2)^2/(2*mass)
      end



     for jc in eachindex(wave)
        k1=k+wave[jc][1]*T1+wave[jc][2]*T2
        pos=findfirst(item->item==wave[jc]-b1T,wave)
        if pos≠nothing
          single_MoirePo[ja][jc,pos]=V0*exp(im*ϕ)*overlap(k1,-b1,spin,mass*vf/2)
        end
        
      
        pos=findfirst(item->item==wave[jc]+b2T+b1T,wave)
        if pos≠nothing
            single_MoirePo[ja][jc,pos]=V0*exp(im*ϕ)*overlap(k1,+b2+b1,spin,mass*vf/2)
        end
    
        pos=findfirst(item->item==wave[jc]-b2T,wave)
        if pos≠nothing
            single_MoirePo[ja][jc,pos]=V0*exp(im*ϕ)*overlap(k1,-b2,spin,mass*vf/2)
        end
     end
    
     single_MoirePo[ja]=single_MoirePo[ja]+single_MoirePo[ja]'
     FFF=eigen(single_MoirePo[ja]+single_Ham[ja])
    
      single_eigenvalue[ja]=real(FFF.values)
      single_eigenvector[ja]=FFF.vectors
    
    end


    input_DensityMatrix=[zeros(ComplexF64,dimension,dimension) for _ in 1:Nx*Ny]
    for ja in 1:Nx*Ny
     
       input_DensityMatrix[ja]+=(single_eigenvector[ja][:,1]*(single_eigenvector[ja][:,1])')
        
    end
     

    

    
    return overlapmatrix, wave, input_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,allowedq, T1, T2, a1m, a2m, b1T,b2T
      

       
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






function Construct_DensityMatrix(loop_dic::Dict{Vector{Int},Any},allowedq::Vector{Vector{Int}},T1::Vector{Float64},T2::Vector{Float64},wave::Vector{Vector{Int64}},input_DensityMatrix::Vector{Matrix{ComplexF64}},single_Ham::Vector{Matrix{ComplexF64}},single_MoirePo::Vector{Matrix{ComplexF64}},constq::Float64,overlapmatrix::Array{ComplexF64,4})::Tuple{Float64,ComplexF64,Vector{Matrix{ComplexF64}},Vector{Matrix{ComplexF64}},Vector{Vector{Float64}},Vector{Matrix{ComplexF64}}}
  
 
   dimension=length(wave)
  HartreeMatrix=[zeros(ComplexF64,dimension,dimension) for _ in eachindex(allowedq)]
  FockMatrix=[zeros(ComplexF64,dimension,dimension) for _ in eachindex(allowedq)]
  output_DensityMatrix=Vector{Matrix{ComplexF64}}(undef,length(allowedq))
  DeltaMatrix=Vector{Matrix{ComplexF64}}(undef,length(allowedq))
  NewDensityMatrix=[zeros(ComplexF64,dimension,dimension) for _ in eachindex(allowedq)]
  HF_eigenvalue=Vector{Vector{Float64}}(undef,length(allowedq))
  HF_eigenvector=Vector{Any}(undef,length(allowedq))

 
  Threads.@threads for jk in eachindex(allowedq)
    Fk = FockMatrix[jk]
    for jk1 in eachindex(allowedq)
        dmk = input_DensityMatrix[jk1]
        q=allowedq[jk1]-allowedq[jk]
       for (dg,loop_dic_dg) in loop_dic
           CoulF1=Coulomb(q+dg,T1,T2)     
           for (gg2,loop_dic_dg_gg2) in loop_dic_dg          
               CoulF=CoulF1*overlapmatrix[jk1,gg2[2],jk,gg2[1]]
           for g1g3 in loop_dic_dg_gg2
               Fk[g1g3[2],gg2[1]]+=dmk[g1g3[1],gg2[2]]*CoulF*overlapmatrix[jk,g1g3[2],jk1,g1g3[1]]
           end 
           end
   
       end    
    end
  end
 
  Hartree_Density=zeros(ComplexF64,dimension,dimension)
  for ja in eachindex(allowedq)
   Hartree_Density+=input_DensityMatrix[ja] .* transpose(overlapmatrix[ja,:,ja,:])
  end

  


  Threads.@threads for jk in eachindex(allowedq)
      for dg in keys(loop_dic)
          CoulH1=Coulomb(dg,T1,T2)
          for gg2 in keys(loop_dic[dg])          
              CoulH=CoulH1*overlapmatrix[jk,gg2[2],jk,gg2[1]]            
          for g1g3 in loop_dic[dg][gg2]           
              HartreeMatrix[jk][gg2[2],gg2[1]]+=Hartree_Density[g1g3[1],g1g3[2]]*CoulH               
          end 
          end
  
      end    
 end


 

 for ja in eachindex(allowedq)
   FFF=eigen(single_MoirePo[ja]+single_Ham[ja]+constq*HartreeMatrix[ja]-constq*FockMatrix[ja])
   HF_eigenvalue[ja]=real(FFF.values)
   HF_eigenvector[ja]=FFF.vectors
   
 end

 bound=sort(reduce(vcat,HF_eigenvalue))[length(allowedq)+1]

  for ja in eachindex(allowedq)
    
       for jd in eachindex(HF_eigenvalue[ja])
          if HF_eigenvalue[ja][jd]<bound
             NewDensityMatrix[ja]+=HF_eigenvector[ja][:,jd]*(HF_eigenvector[ja][:,jd])'
          end
       end
       DeltaMatrix[ja]=NewDensityMatrix[ja]-input_DensityMatrix[ja]
       output_DensityMatrix[ja]=0.0*input_DensityMatrix[ja]+1.0*NewDensityMatrix[ja]
  end


  
  e1=0.0
  energy=0.0
  for ja in eachindex(allowedq)
    e1+=tr(DeltaMatrix[ja]'*DeltaMatrix[ja])
    energy+=tr((single_MoirePo[ja]+single_Ham[ja]+constq/2*HartreeMatrix[ja]-constq/2*FockMatrix[ja])*output_DensityMatrix[ja])
  end
  eout=real(e1)/length(allowedq)
  


  
 

 return  eout,energy,output_DensityMatrix,DeltaMatrix,HF_eigenvalue,HF_eigenvector
end




function iteration_loop(initial_DensityMatrix::Vector{Matrix{ComplexF64}},allowedq::Vector{Vector{Int64}},T1::Vector{Float64},T2::Vector{Float64},wave::Vector{Vector{Int64}},single_Ham::Vector{Matrix{ComplexF64}},single_MoirePo::Vector{Matrix{ComplexF64}},constq::Float64,overlapmatrix::Array{ComplexF64,4})::Tuple{Vector{Vector{Matrix{ComplexF64}}},Vector{Vector{Matrix{ComplexF64}}},Vector{Vector{Float64}},Vector{Matrix{ComplexF64}},ComplexF64}
    eout=1.0
    energy=0.0

    itcount=0
    loop_dic=construct_loop_dic(wave)
    HF_eigenvalue=Vector{Any}(undef,length(allowedq))
    HF_eigenvector=Vector{Any}(undef,length(allowedq))
    DIIS_input_DensityMatrix=Vector{Vector{Matrix{ComplexF64}}}(undef,3)
    DIIS_input_DeltaMatrix=Vector{Vector{Matrix{ComplexF64}}}(undef,3)
    input_DensityMatrix=initial_DensityMatrix
    bad_count=0
    while (eout>1*10^-8) || (bad_count<4)
        if  eout<1*10^-8 
            bad_count+=1
        end
      tic=time()
      eout,energy,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue,HF_eigenvector=Construct_DensityMatrix(loop_dic,allowedq,T1,T2,wave,input_DensityMatrix,single_Ham,single_MoirePo,constq,overlapmatrix)
      DIIS_input_DensityMatrix[mod(itcount,3)+1]=input_DensityMatrix
      input_DensityMatrix=output_DensityMatrix
      itcount+=1
      toc=time()
      println(toc-tic,"eout=$eout")
      flush(stdout)
     
    
    end
    
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
            for jc in eachindex(allowedq)
               Bmatrix[ja,jb]+=tr((DIIS_input_DeltaMatrix[ja][jc])'*(DIIS_input_DeltaMatrix[jb][jc]))
            end
        end
        coeff=inv(Bmatrix)*[0;0;0;1]
       
        dmk=coeff[1]*(DIIS_input_DensityMatrix[1]+DIIS_input_DeltaMatrix[1])+coeff[2]*(DIIS_input_DensityMatrix[2]+DIIS_input_DeltaMatrix[2])+coeff[3]*(DIIS_input_DensityMatrix[3]+DIIS_input_DeltaMatrix[3])
         eout,energy,_,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue,HF_eigenvector=Construct_DensityMatrix(loop_dic,allowedq,T1,T2,wave,dmk,single_Ham,single_MoirePo,constq,overlapmatrix)
        DIIS_input_DensityMatrix[mod(itcount,3)+1]=dmk
        itcount+=1

        toc=time()
        println(toc-tic,"eout=$eout")
        flush(stdout)
    end

 


  return DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy

end


function Densitymap(a1m::Vector{Float64},a2m::Vector{Float64},overlapmatrix::Array{ComplexF64,4},wave::Vector{Vector{Int}},input_DensityMatrix::Vector{Matrix{ComplexF64}})
    N3=50
    dimension=length(wave)
    zgrid=zeros(Float64,N3,N3)
    Hartree_Density=zeros(ComplexF64,dimension,dimension)
    for jk1 in eachindex(input_DensityMatrix)
       Hartree_Density+=input_DensityMatrix[jk1] .* transpose(overlapmatrix[jk1,:,jk1,:])
    end
    
    for ja in 1:50, jb in 1:50
      rvec=ja/25*a1m+jb/25*a2m
      for jc in eachindex(wave), jd in eachindex(wave)
        gvec=[T1 T2]*(wave[jc]-wave[jd])
       zgrid[ja,jb]+=real(Hartree_Density[jc,jd]*exp(im*(gvec[1]*rvec[1]+gvec[2]*rvec[2])))
      end
    
    end
    
    return zgrid
end



  




function sendtomesh(Minv::Matrix{Int64},q1::Vector{Int})::Vector{Int}
    Qvec=q1'*inv(Minv)
    return q1 .-vec((Int.(floor.(round.(Qvec,digits=5)))*Minv)')
end






function Geometry(geonum::Int64)

    

 
    if geonum==1
     Nx=4;
     Ny=4;
     l1=[4,0]
     l2=[0,4]
    end
 

 
    if geonum==2
     Nx=4;
     Ny=5;
     l1=[4,0]
     l2=[0,5]
    end
 
    if geonum==3
     Nx=4;
     Ny=6;
     l1=[4,0]
     l2=[0,6]
    end
     
    if geonum==4
     Nx=1
     Ny=24;
     l1=[1,4]
     l2=[5,-4]
    end
    
 
    if geonum==5
     Nx=4
     Ny=7;
     l1=[4,0]
     l2=[0,7]
    end
    
    if geonum==6
        Nx=5
        Ny=6;
        l1=[5,0]
        l2=[0,6]
    end


    if geonum==8
        Nx=2
        Ny=48;
        l1=[2,8]
        l2=[10,-8]
    end

    if geonum==9
        Nx=8;
        Ny=12;
        l1=[8,0]
        l2=[0,12]
    end

    if geonum==10
        Nx=8;
        Ny=8;
        l1=[8,0]
        l2=[0,8]
    end
    
    
 
    return Nx,Ny,l1,l2
 end