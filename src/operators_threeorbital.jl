function Construct_DensityMatrix_Hartreeonly(allowedq::Vector{Vector{Int}},T1::Vector{Float64},T2::Vector{Float64},Nq::Int64,input_DensityMatrix::Vector{Matrix{ComplexF64}},single_Ham::Vector{Matrix{ComplexF64}},U::Float64,Uprime::Float64,nu::Int64)::Tuple{Float64,Vector{Matrix{ComplexF64}},Vector{Matrix{ComplexF64}},Vector{Vector{Float64}},Matrix{ComplexF64},ComplexF64,Vector{Matrix{ComplexF64}}}

 
    FockMatrix=[zeros(ComplexF64,3,3) for _ in 1:Nq^2]
    HartreeMatrix=zeros(ComplexF64,3,3)
    output_DensityMatrix=Vector{Matrix{ComplexF64}}(undef,Nq^2)
    DeltaMatrix=Vector{Matrix{ComplexF64}}(undef,Nq^2)
    NewDensityMatrix=[zeros(ComplexF64,3,3) for _ in 1:Nq^2]
    HF_eigenvalue=Vector{Vector{Float64}}(undef,Nq^2)
    HF_eigenvector=Vector{Any}(undef,Nq^2)
   
    HartreeDensity=zeros(ComplexF64,3,3)
    for ja in 1:Nq^2
      HartreeDensity+=input_DensityMatrix[ja]
    end
   
   HartreeMatrix[1,1]+=HartreeDensity[2,2]*U/Nq^2*3
   HartreeMatrix[2,2]+=HartreeDensity[1,1]*U/Nq^2*3
   HartreeMatrix[1,1]+=HartreeDensity[3,3]*U/Nq^2*3
   HartreeMatrix[3,3]+=HartreeDensity[1,1]*U/Nq^2*3
   HartreeMatrix[2,2]+=HartreeDensity[3,3]*Uprime/Nq^2*3
   HartreeMatrix[3,3]+=HartreeDensity[2,2]*Uprime/Nq^2*3
   
   
    

   BmA=[[0,1],[0,0],[√3/2,1/2]]#from A to it's neighbour B
   CmA=[[0,1],[0,0],[-√3/2,1/2]]#from A to it's neighbour C
   BmC=[[√3/2,-1/2],[0,0],[√3/2,1/2]]#from C to it's neighbour B

   
  
   
   
   for ja in 1:Nq^2
     FFF=eigen(single_Ham[ja]+HartreeMatrix)
     HF_eigenvalue[ja]=real(FFF.values)
     HF_eigenvector[ja]=FFF.vectors
     
   end
   
   
   
   bound=(sort(reduce(vcat,HF_eigenvalue))[nu*Nq^2+1]+sort(reduce(vcat,HF_eigenvalue))[nu*Nq^2])/2
   
    for ja in 1:Nq^2
      
         for jd in eachindex(HF_eigenvalue[ja])
            if HF_eigenvalue[ja][jd]<bound
               NewDensityMatrix[ja]+=HF_eigenvector[ja][:,jd]*(HF_eigenvector[ja][:,jd])'
            end
         end
         DeltaMatrix[ja]=NewDensityMatrix[ja]-input_DensityMatrix[ja]
         output_DensityMatrix[ja]=0.0*input_DensityMatrix[ja]+1.0*NewDensityMatrix[ja]
    end
   
   
    
    e1=0.0
    for ja in 1:Nq^2
      e1+=tr(DeltaMatrix[ja]'*DeltaMatrix[ja])
    end
    eout=real(e1)/Nq^2
    
    energy=0.0
    for ja in 1:Nq^2
        energy+=tr(output_DensityMatrix[ja]*(single_Ham[ja]+1/2*HartreeMatrix))
    end 

    if eout<10^(-22)


     for ja in 1:Nq^2
       k=[T1 T2]*allowedq[ja]
       for jb in 1:Nq^2
           kprime=[T1 T2]*allowedq[jb]
           for deltaBA in BmA
              FockMatrix[ja][1,2]+=input_DensityMatrix[jb][1,2]*U/Nq^2*exp(im*dot(k-kprime,deltaBA))
           end
  
           for deltaCA in CmA
               FockMatrix[ja][1,3]+=input_DensityMatrix[jb][1,3]*U/Nq^2*exp(im*dot(k-kprime,deltaCA))
           end
       
           for deltaBC in BmC
               FockMatrix[ja][3,2]+=input_DensityMatrix[jb][3,2]*Uprime/Nq^2*exp(im*dot(k-kprime,deltaBC))
            end
       end
  
        FockMatrix[ja]=FockMatrix[ja]+FockMatrix[ja]'
      end
    
      energy=0.0
      for ja in 1:Nq^2
          energy+=tr(output_DensityMatrix[ja]*(single_Ham[ja]+1/2*HartreeMatrix-1/2*FockMatrix[ja]))
      end 
    end


   
   
   return  eout,output_DensityMatrix,DeltaMatrix,HF_eigenvalue, HartreeMatrix, energy, HF_eigenvector
end
   


function Construct_DensityMatrix(allowedq::Vector{Vector{Int}},T1::Vector{Float64},T2::Vector{Float64},Nq::Int64,input_DensityMatrix::Vector{Matrix{ComplexF64}},single_Ham::Vector{Matrix{ComplexF64}},U::Float64,Uprime::Float64,nu::Int64)::Tuple{Float64,Vector{Matrix{ComplexF64}},Vector{Matrix{ComplexF64}},Vector{Vector{Float64}},Matrix{ComplexF64},ComplexF64,Vector{Matrix{ComplexF64}}}

 
  
    HartreeMatrix=zeros(ComplexF64,3,3)
    FockMatrix=[zeros(ComplexF64,3,3) for _ in 1:Nq^2]
    output_DensityMatrix=Vector{Matrix{ComplexF64}}(undef,Nq^2)
    DeltaMatrix=Vector{Matrix{ComplexF64}}(undef,Nq^2)
    NewDensityMatrix=[zeros(ComplexF64,3,3) for _ in 1:Nq^2]
    HF_eigenvalue=Vector{Vector{Float64}}(undef,Nq^2)
    HF_eigenvector=Vector{Any}(undef,Nq^2)
   
    HartreeDensity=zeros(ComplexF64,3,3)
    for ja in 1:Nq^2
      HartreeDensity+=input_DensityMatrix[ja]
    end
   
   HartreeMatrix[1,1]+=HartreeDensity[2,2]*U/Nq^2*3
   HartreeMatrix[2,2]+=HartreeDensity[1,1]*U/Nq^2*3
   HartreeMatrix[1,1]+=HartreeDensity[3,3]*U/Nq^2*3
   HartreeMatrix[3,3]+=HartreeDensity[1,1]*U/Nq^2*3
   HartreeMatrix[2,2]+=HartreeDensity[3,3]*Uprime/Nq^2*3
   HartreeMatrix[3,3]+=HartreeDensity[2,2]*Uprime/Nq^2*3
   
   
   BmA=[[0,1],[0,0],[√3/2,1/2]]#from A to it's neighbour B
   CmA=[[0,1],[0,0],[-√3/2,1/2]]#from A to it's neighbour C
   BmC=[[√3/2,-1/2],[0,0],[√3/2,1/2]]#from C to it's neighbour B

   
   for ja in 1:Nq^2
    k=[T1 T2]*allowedq[ja]
    for jb in 1:Nq^2
     kprime=[T1 T2]*allowedq[jb]
     for deltaBA in BmA
     FockMatrix[ja][1,2]+=input_DensityMatrix[jb][1,2]*U/Nq^2*exp(im*dot(k-kprime,deltaBA))
     end

     for deltaCA in CmA
     FockMatrix[ja][1,3]+=input_DensityMatrix[jb][1,3]*U/Nq^2*exp(im*dot(k-kprime,deltaCA))
     end
     
     for deltaBC in BmC
     FockMatrix[ja][3,2]+=input_DensityMatrix[jb][3,2]*Uprime/Nq^2*exp(im*dot(k-kprime,deltaBC))
     end
    end

    FockMatrix[ja]=FockMatrix[ja]+FockMatrix[ja]'
   end
   
   
   
   for ja in 1:Nq^2
     FFF=eigen(single_Ham[ja]+HartreeMatrix-FockMatrix[ja])
     HF_eigenvalue[ja]=real(FFF.values)
     HF_eigenvector[ja]=FFF.vectors
     
   end
   
   
   
   bound=(sort(reduce(vcat,HF_eigenvalue))[nu*Nq^2+1]+sort(reduce(vcat,HF_eigenvalue))[nu*Nq^2])/2
   
    for ja in 1:Nq^2
      
         for jd in eachindex(HF_eigenvalue[ja])
            if HF_eigenvalue[ja][jd]<bound
               NewDensityMatrix[ja]+=HF_eigenvector[ja][:,jd]*(HF_eigenvector[ja][:,jd])'
            end
         end
         DeltaMatrix[ja]=NewDensityMatrix[ja]-input_DensityMatrix[ja]
         output_DensityMatrix[ja]=0.0*input_DensityMatrix[ja]+1.0*NewDensityMatrix[ja]
    end
   
   
    
    e1=0.0
    for ja in 1:Nq^2
      e1+=tr(DeltaMatrix[ja]'*DeltaMatrix[ja])
    end
    eout=real(e1)/Nq^2
    
    energy=0.0
    for ja in 1:Nq^2
        energy+=tr(output_DensityMatrix[ja]*(single_Ham[ja]+1/2*HartreeMatrix-1/2*FockMatrix[ja]))
    end
    
    
   
   return  eout,output_DensityMatrix,DeltaMatrix,HF_eigenvalue, HartreeMatrix, energy, HF_eigenvector
end


function get_Ham(ϵA::Float64,ϵB::Float64,ϵC::Float64,t1::Float64,t2::Float64,t3::Float64,k::Vector{Float64})
   
    H=zeros(ComplexF64,3,3)
   
    H[3,1]+=t2*exp(im*π/3)+t2*exp(-im*π/3)*exp(im*dot(k,[0,-1]))+t2*exp(-im*π)*exp(im*dot(k,[√3/2,-1/2]))
    H[2,1]+=t1*exp(im*dot(k,[-√3/2,-1/2]))+t1*exp(-im*2*π/3)*exp(im*dot(k,[0,-1]))+t1*exp(-im*4*π/3)
    H[2,3]+=t3+t3*exp(im*dot(k,[-√3/2,-1/2]))+t3*exp(im*dot(k,[-√3/2,1/2]))

    H=H+H'
 
    H[1,1]+=ϵA
    H[2,2]+=ϵB
    H[3,3]+=ϵC
 
    return H
 
end


function get_input(ϵA::Float64,ϵB::Float64,ϵC::Float64,t1::Float64,t2::Float64,t3::Float64,nu::Int,Nq::Int)

    a1=[√3/2,-1/2]
    a2=[√3/2,1/2]
    

    b1=4*π/√3*[1/2,-√3/2]
    b2=4*π/√3*[1/2,√3/2]
    
    T1=b1/(Nq)
    T2=b2/(Nq)
    
    
    
    b1T=Int.(round.(inv([T1 T2])*b1))
    b2T=Int.(round.(inv([T1 T2])*b2))
    
    
    allowedq=Vector{Int64}[]
    for ja in 0:Nq-1,jb in 0:Nq-1
        push!(allowedq,[ja,jb])
    end
    
    
    single_Ham=[zeros(ComplexF64,3,3) for _ in 1:Nq^2]
    single_eigenvector=[zeros(ComplexF64,3,3) for _ in 1:Nq^2]
    single_eigenvalue=[zeros(Float64,3) for _ in 1:Nq^2]
    input_DensityMatrix=[zeros(ComplexF64,3,3) for _ in 1:Nq^2]
    single_DensityMatrix=[zeros(ComplexF64,3,3) for _ in 1:Nq^2]
    
    
    for ja in 1:Nq^2
        k=[T1 T2]*allowedq[ja]
        single_Ham[ja]=get_Ham(ϵA,ϵB,ϵC,t1,t2,t3,k)
        single_eigenvector[ja]=eigen(single_Ham[ja]).vectors
        single_eigenvalue[ja]=eigen(single_Ham[ja]).values
    
        if ja==1
        
            single_eigenvector[ja][:,2]=[0.0,1/√2,-1/√2]
            single_eigenvector[ja][:,1]=[0.0,1/√2,1/√2]
        end
        
    end
    
    for ja in 1:Nq^2,jb in 1:nu
        input_DensityMatrix[ja]+=(single_eigenvector[ja][:,jb]*(single_eigenvector[ja][:,jb])') 
        single_DensityMatrix[ja]+=(single_eigenvector[ja][:,jb]*(single_eigenvector[ja][:,jb])') 
    end
    
    
   
    return single_DensityMatrix,input_DensityMatrix,single_eigenvalue,single_Ham,allowedq,T1,T2
    

end




function iteration_loop(Nq::Int64,nu::Int64,allowedq::Vector{Vector{Int64}},T1::Vector{Float64},T2::Vector{Float64},single_Ham::Vector{Matrix{ComplexF64}},U::Float64,Uprime::Float64,input_DensityMatrix::Vector{Matrix{ComplexF64}})
    eout=1.0
    itcount=0
    energy=0
    HF_eigenvalue=[zeros(Float64,3) for _ in 1:Nq^2]
    HF_eigenvector=[zeros(ComplexF64,3,3) for _ in 1:Nq^2];
    DeltaMatrix=[zeros(ComplexF64,3,3) for _ in 1:Nq^2]
    DIIS_input_DensityMatrix=Vector{Vector{Matrix{ComplexF64}}}(undef,3)
    DIIS_input_DeltaMatrix=Vector{Vector{Matrix{ComplexF64}}}(undef,3)
    
    Hartree_Matrix=zeros(ComplexF64,3)
    bad_count=0
    bound=0.0

    while (eout>10^-22) || (bad_count<4)
        if eout<10^-22
            bad_count+=1 
        end
      #tic=time()
      eout,output_DensityMatrix,DeltaMatrix,HF_eigenvalue,Hartree_Matrix,energy,HF_eigenvector=Construct_DensityMatrix_Hartreeonly(allowedq,T1,T2,Nq,input_DensityMatrix,single_Ham,U,Uprime,nu)
      DIIS_input_DensityMatrix[mod(itcount,3)+1]=input_DensityMatrix
      DIIS_input_DeltaMatrix[mod(itcount,3)+1]=DeltaMatrix
      input_DensityMatrix=output_DensityMatrix
      itcount+=1
      #toc=time()
     #println(toc-tic,"eout=$eout")
     #flush(stdout)
     if itcount>30000
        break
     end
    
    end
    println(eout)
    flush(stdout)

    return DIIS_input_DensityMatrix[1],eout,energy,HF_eigenvalue,HF_eigenvector
   
end


function calculate_observable(Nq::Int64,input_DM::Vector{Matrix{ComplexF64}},T1::Vector{Float64},T2::Vector{Float64},allowedq::Vector{Vector{Int}})
     #Let's caluclate CdaggerC, BdaggerB, AdaggerA first
     BdB=0
   for ja in 1:Nq^2
      BdB+=input_DM[ja][2,2]/Nq^2
   end

    AdA=0
    for ja in 1:Nq^2
      AdA+=input_DM[ja][1,1]/Nq^2
    end

    CdC=0
   for ja in 1:Nq^2
     CdC+=input_DM[ja][3,3]/Nq^2
   end


   #We calculate BdA,CdB,CdA
   BmA=[[0,1],[0,0],[√3/2,1/2]]#from A to it's neighbour B
   CmA=[[0,1],[0,0],[-√3/2,1/2]]#from A to it's neighbour C
   BmC=[[√3/2,-1/2],[0,0],[√3/2,1/2]]#from C to it's neighbour B

   BdA=0
 for ja in 1:Nq^2
    kvec=[T1 T2]*allowedq[ja]
    #BdA+=input_DM[ja][1,2]*exp(-im*dot(kvec,BmA[3]))/Nq^2
    BdA+=input_DM[ja][1,2]/Nq^2
 end

  CdB=0
  for ja in 1:Nq^2
     kvec=[T1 T2]*allowedq[ja]
     #CdB+=input_DM[ja][2,3]*exp(im*dot(kvec,BmC[3]))/Nq^2
     CdB+=input_DM[ja][2,3]/Nq^2
  end

  AdC=0
  for ja in 1:Nq^2
     kvec=[T1 T2]*allowedq[ja]
     #AdC+=input_DM[ja][3,1]*exp(im*dot(kvec,CmA[3]))/Nq^2
     AdC+=input_DM[ja][3,1]/Nq^2
  end

 

   
  return AdA,BdB,CdC,BdA,CdB,AdC






end



function excecute_loop()
   # ϵAspace=collect(0.0:0.2:6.0)
    ϵAspace=collect(1.0:0.5:9.0)
    ϵBspace=[0.0]
    ϵCspace=[0.0]
    t1space=[1.0]
    t2space=[1.0]
    t3space=[0.0]
    Nq=15
    Uspace=collect(0.05:0.05:1.0)
    #Uspace=[0.1]
    nu=1
    allowedq=0

    single_eigenvalue=[[[zeros(ComplexF64,3) for _ in 1:Nq^2] for _ in eachindex(Uspace)] for _ in eachindex(ϵAspace)]
    HF_eigenvalue=[[[zeros(ComplexF64,3) for _ in 1:Nq^2] for _ in eachindex(Uspace)] for _ in eachindex(ϵAspace)]

    eoutmatrix=[zeros(ComplexF64,length(Uspace)) for _ in eachindex(ϵAspace)]
    energymatrix=[zeros(ComplexF64,length(Uspace)) for _ in eachindex(ϵAspace)]
    AdAmatrix=[zeros(ComplexF64,length(Uspace)) for _ in eachindex(ϵAspace)]
    BdBmatrix=[zeros(ComplexF64,length(Uspace)) for _ in eachindex(ϵAspace)]
    CdCmatrix=[zeros(ComplexF64,length(Uspace)) for _ in eachindex(ϵAspace)]
    AdCmatrix=[zeros(ComplexF64,length(Uspace)) for _ in eachindex(ϵAspace)]
    BdAmatrix=[zeros(ComplexF64,length(Uspace)) for _ in eachindex(ϵAspace)]
    CdBmatrix=[zeros(ComplexF64,length(Uspace)) for _ in eachindex(ϵAspace)]

    single_AdAmatrix=[zeros(ComplexF64,length(Uspace)) for _ in eachindex(ϵAspace)]
    single_BdBmatrix=[zeros(ComplexF64,length(Uspace)) for _ in eachindex(ϵAspace)]
    single_CdCmatrix=[zeros(ComplexF64,length(Uspace)) for _ in eachindex(ϵAspace)]
    single_AdCmatrix=[zeros(ComplexF64,length(Uspace)) for _ in eachindex(ϵAspace)]
    single_BdAmatrix=[zeros(ComplexF64,length(Uspace)) for _ in eachindex(ϵAspace)]
    single_CdBmatrix=[zeros(ComplexF64,length(Uspace)) for _ in eachindex(ϵAspace)]

  Threads.@threads for ja in eachindex(ϵAspace)
        for jb in eachindex(Uspace)
            ϵA=ϵAspace[ja]
            U=Uspace[jb]
            Uprime=0.9*U
            ϵB=ϵBspace[1]
            ϵC=ϵCspace[1]
            t1=t1space[1]
            t2=t2space[1]
            t3=t3space[1]

            single_DensityMatrix,input_DensityMatrix,single_eigenvalue[ja][jb],single_Ham,allowedq,T1,T2=get_input(ϵA,ϵB,ϵC,t1,t2,t3,nu,Nq)
            for jc in 1:Nq^2
                A=randn(3,3)+im*randn(3,3)
                input_DensityMatrix[jc]+=(A+A')*0.1
            end

            single_AdAmatrix[ja][jb],single_BdBmatrix[ja][jb],single_CdCmatrix[ja][jb],single_BdAmatrix[ja][jb],single_CdBmatrix[ja][jb],single_AdCmatrix[ja][jb]=calculate_observable(Nq,single_DensityMatrix,T1,T2,allowedq)
     
            HF_DensityMatrix,eoutmatrix[ja][jb],energymatrix[ja][jb],HF_eigenvalue[ja][jb],HF_eigenvector=iteration_loop(Nq,nu,allowedq,T1,T2,single_Ham,U,Uprime,input_DensityMatrix)
          
            AdAmatrix[ja][jb],BdBmatrix[ja][jb],CdCmatrix[ja][jb],BdAmatrix[ja][jb],CdBmatrix[ja][jb],AdCmatrix[ja][jb]=calculate_observable(Nq,HF_DensityMatrix,T1,T2,allowedq)
       
        end
    end




    observable=Dict{String,Any}()
    observable["AdAmatrix"]=AdAmatrix
    observable["BdBmatrix"]=BdBmatrix
    observable["CdCmatrix"]=CdCmatrix
    observable["AdCmatrix"]=AdCmatrix
    observable["BdAmatrix"]=BdAmatrix
    observable["CdBmatrix"]=CdBmatrix

    observable["single_AdAmatrix"]=single_AdAmatrix
    observable["single_BdBmatrix"]=single_BdBmatrix
    observable["single_CdCmatrix"]=single_CdCmatrix
    observable["single_AdCmatrix"]=single_AdCmatrix
    observable["single_BdAmatrix"]=single_BdAmatrix
    observable["single_CdBmatrix"]=single_CdBmatrix


    parameters=Dict{String,Any}()
    parameters["ϵAspace"]=ϵAspace
    parameters["ϵBspace"]=ϵBspace
    parameters["ϵCspace"]=ϵCspace
    parameters["t1space"]=t1space
    parameters["t2space"]=t2space
    parameters["t3space"]=t3space
    parameters["Nq"]=Nq
    parameters["Uspace"]=Uspace
    parameters["nu"]=nu


  return observable,parameters,energymatrix,eoutmatrix,single_eigenvalue,HF_eigenvalue
end

function test_func()
   A=[zeros(Float64,2) for _ in 1:2]
   #print(A)
   B=copy(A)
   for ja in 1:2, jb in 1:2
       (B[ja][jb],A[ja][jb])=(1.0,2.0)
       
   end
   print(B)
   print(A)
end