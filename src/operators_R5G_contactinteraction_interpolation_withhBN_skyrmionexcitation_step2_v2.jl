using LinearAlgebra
BLAS.set_num_threads(1)
using Arpack
using Combinatorics
using Random








function build_k_loopup(allowedq::Vector{Vector{Int}},allowedq_dic::Dict{Vector{Int}},Nq::Int)
  sum_idx=Matrix{Int}(undef,length(allowedq),length(allowedq))
  diff_idx=Matrix{Int}(undef,length(allowedq),length(allowedq))

  for i in eachindex(allowedq), j in eachindex(allowedq)
     sum_idx[i,j]=allowedq_dic[mod.(allowedq[i]+allowedq[j],Nq)]
     diff_idx[i,j]=allowedq_dic[mod.(allowedq[i]-allowedq[j],Nq)]
  end

  return sum_idx,diff_idx
end




function Coulomb(k::Vector{Int},T1::Vector{Float64},T2::Vector{Float64})::Float64
   D=25
   return k==[0,0] ? D*9047.5636 : tanh(norm([T1 T2]*k*D))/norm(k[1]*T1+k[2]*T2)*9047.5636
end

function Construct_DensityMatrix(allowedq::Vector{Vector{Int}},
                                 Nq::Int64,Nband::Int64,wave_diff::Vector{Vector{Int64}},
                               input_DensityMatrix::Matrix{ComplexF64},single_Ham::Matrix{ComplexF64},
                               single_MoirePo::Matrix{ComplexF64},Coulomb_element::Matrix{Float64},formf::Array{Matrix{ComplexF64}},
                               energy_input::Float64,filling::Int,Area::Float64,sum_idx::Matrix{Int},diff_idx::Matrix{Int})
  

  output_DensityMatrix=zeros(ComplexF64,Nq^2*Nband,Nq^2*Nband)
  DeltaMatrix=zeros(ComplexF64,Nq^2*Nband,Nq^2*Nband)
  NewDensityMatrix=zeros(ComplexF64,Nq^2*Nband,Nq^2*Nband)
  HF_eigenvalue=zeros(Float64,Nq^2*Nband)
  HF_eigenvector=zeros(ComplexF64,Nq^2*Nband,Nq^2*Nband)
 
   num_spin=2
  Nbs=Int(Nband/num_spin)
  input_ds_reshaped=reshape(input_DensityMatrix,Nq^2,Nbs,num_spin,Nq^2,Nbs,num_spin)
   
 
  HartreeMatrix=zeros(ComplexF64,Nq^2,Nbs,num_spin,Nq^2,Nbs,num_spin)
  FockMatrix=zeros(ComplexF64,Nq^2,Nbs,num_spin,Nq^2,Nbs,num_spin)
   nt = Threads.maxthreadid()
   Xbuf = [Matrix{ComplexF64}(undef, Nbs, Nbs) for _ in 1:nt]
  Ybuf = [Matrix{ComplexF64}(undef, Nbs, Nbs) for _ in 1:nt]
  ρbuf = [Matrix{ComplexF64}(undef, Nbs, Nbs) for _ in 1:nt]
  tic=time()
  Threads.@threads for k2 in eachindex(allowedq)
    tid=Threads.threadid()
    ρb = ρbuf[tid]
    X=Xbuf[tid]
    Y=Ybuf[tid]
    for k3 in eachindex(allowedq),s1 in 1:num_spin, s2 in 1:num_spin
    Fk_local=@view  FockMatrix[k2,:,s1,k3,:,s2]
    Fk_buf=zeros(ComplexF64,Nbs,Nbs) 
       for k1 in eachindex(allowedq)
          #qindex=allowedq_dic[mod.(allowedq[k1]-allowedq[k3],Nq)]
          #k4=allowedq_dic[mod.(allowedq[k1]+allowedq[k2]-allowedq[k3],Nq)]
           qindex=diff_idx[k1,k3]
           k4=sum_idx[qindex,k2]
           ρview = @view  input_ds_reshaped[k4, :, s1, k1, :, s2]
           copyto!(ρb, ρview)
           for qg in eachindex(wave_diff)
           
             α=Coulomb_element[qindex,qg]  
             A=formf[k2,qindex,qg,s1]   
             B=formf[k3,qindex,qg,s2]
             mul!(X,ρb,B')
             mul!(Y,A,X)
             Fk_buf.+=α.*Y
            end
         
       end
    copyto!(Fk_local, Fk_buf)
      end
  end
  toc=time()
  println(toc-tic,"finisheFock")
  
  tic=time()

  Threads.@threads for k2 in eachindex(allowedq)
    for k4 in eachindex(allowedq)
    HT_local=@view HartreeMatrix[k2,:,1,k4,:,1]
    #qindex=allowedq_dic[mod.(allowedq[k4]-allowedq[k2],Nq)]
     qindex=diff_idx[k4,k2]
     for qg in eachindex(wave_diff)
        htdensity=zero(ComplexF64)
        for k1 in eachindex(allowedq), si in 1:num_spin
          #k3=allowedq_dic[mod.(allowedq[k1]+allowedq[k2]-allowedq[k4],Nq)]
          k3=diff_idx[k1,qindex]
          #htdensity+=tr(input_ds_reshaped[k3,:,k1,:]*formf[k3,qindex,qg]')
          htdensity+=dot(formf[k3,qindex,qg,si],input_ds_reshaped[k3,:,si,k1,:,si])
        end
        β=Coulomb_element[qindex,qg]*htdensity
      
         HT_local.+=β.*formf[k2,qindex,qg,1] #I used the knowledge here that Hartree are the same for both spin
      
      end
   end
  end
  HartreeMatrix[:,:,2,:,:,2]=copy(HartreeMatrix[:,:,1,:,:,1])

  toc=time()
  println(toc-tic,"finishHartree")

   
   

  HartreeMatrix=reshape(HartreeMatrix,Nq^2*Nband,Nq^2*Nband)
  FockMatrix=reshape(FockMatrix,Nq^2*Nband,Nq^2*Nband)


  HartreeMatrix=(HartreeMatrix+HartreeMatrix')/(2*Area)
  FockMatrix=(FockMatrix+FockMatrix')/(2*Area)

 


   FFF=eigen(single_MoirePo+single_Ham+HartreeMatrix-FockMatrix)
   HF_eigenvalue=real(FFF.values)
   HF_eigenvector=FFF.vectors


 bound=(sort(reduce(vcat,HF_eigenvalue))[filling+1]+sort(reduce(vcat,HF_eigenvalue))[filling])/2


    
       for jd in eachindex(HF_eigenvalue)
          if HF_eigenvalue[jd]<bound
             NewDensityMatrix+=HF_eigenvector[:,jd]*(HF_eigenvector[:,jd])'
          end
       end
        NewDensityMatrix=0.5*(NewDensityMatrix'+ NewDensityMatrix)
       DeltaMatrix=NewDensityMatrix-input_DensityMatrix
       mix_ratio=0.5
       output_DensityMatrix=mix_ratio*input_DensityMatrix+(1-mix_ratio)*NewDensityMatrix


   


  
  eout = sum(abs2, DeltaMatrix)/Nq^2


  energy = real(dot(input_DensityMatrix, single_MoirePo)) +real(dot(input_DensityMatrix, single_Ham)) +  real(dot(input_DensityMatrix, pinning_po))+ 0.5*real(dot(input_DensityMatrix, HartreeMatrix)) -0.5*real(dot(input_DensityMatrix, FockMatrix))
    
   energy_change=real(energy-energy_input)



 return  eout,energy_change,output_DensityMatrix,DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,HartreeMatrix,FockMatrix
end


function get_ff(nop_HF_eigenvector::Vector{Matrix{ComplexF64}},NL::Int,Nband::Int,wave::Vector{Vector{Int}}
                ,wave_diff::Vector{Vector{Int}},spinor_set::Matrix{Vector{ComplexF64}},allowedq::Vector{Vector{Int}},allowedq_dic::Dict{Vector{Int}})
  num_spin=2
  Nbs=Int(Nband/num_spin)
  form_factors=[zeros(ComplexF64,Nbs,Nbs) for _ in eachindex(allowedq), _ in eachindex(allowedq), _ in eachindex(wave_diff), _ in 1:num_spin]

  eigenvector_PWbasis=[zeros(ComplexF64,2*NL,length(wave),Nbs) for _ in eachindex(allowedq)]
  diff_eigenvector_threads=[zeros(ComplexF64,2*NL,length(wave),Nbs,length(allowedq)) for _ in eachindex(wave_diff)]
  
  
 
  
 

  for ja in eachindex(allowedq)
   
    for jc in eachindex(wave), bi in 1:Nbs
      eigenvector_PWbasis[ja][:,jc,bi]=nop_HF_eigenvector[ja][jc,bi]*spinor_set[ja,jc]
    end
    
  end


   Threads.@threads for jd in eachindex(wave_diff)
      for ja in eachindex(wave)
            pos=findfirst(item->item==wave[ja]-wave_diff[jd],wave)
            if pos≠nothing
              for jk in eachindex(allowedq)
                 diff_eigenvector_threads[jd][:,ja,:,jk]=eigenvector_PWbasis[jk][:,pos,:]
              end
            end
      end
   end
   


   eigenvector_PWbasis=[reshape(eigenvector_PWbasis[ja],2*NL*length(wave),Nbs) for ja in eachindex(allowedq)]
   diff_eigenvector_threads=[reshape(diff_eigenvector_threads[ja],2*NL*length(wave),Nbs,length(allowedq)) for ja in eachindex(wave_diff)] 




   Threads.@threads for ja in 1:Nq^2
    ff_local=@view  form_factors[ja,:,:,:]
    for jb in 1:Nq^2
     
        meshk1plusq=mod.(allowedq[ja]+allowedq[jb],Nq)
        k2pos=allowedq_dic[meshk1plusq]
       for  jc in eachindex(wave_diff)
           gk1plusq=wave_diff[jc]+allowedq[ja]+allowedq[jb]-meshk1plusq
           gk1plusq_pos=findfirst(item->item==gk1plusq,wave_diff)
           
           if gk1plusq_pos≠nothing
           
             
                ff_local[jb,jc,1]=(diff_eigenvector_threads[gk1plusq_pos][:,:,ja])'*eigenvector_PWbasis[k2pos]
         
           end
       end
     end
   end
     form_factors[:,:,:,2]=form_factors[:,:,:,1]


   return form_factors
end

function initial_process(allowedq::Vector{Vector{Int}},nop_single_Ham::Vector{Matrix{ComplexF64}},
                         nop_single_MoirePo::Vector{Matrix{ComplexF64}},nop_HF_eigenvector::Vector{Matrix{ComplexF64}},
                         Nband::Int64,qcutoff::Float64,b1::Vector{Float64},b2::Vector{Float64},b1T::Vector{Int},b2T::Vector{Int})
  allowedq_dic=Dict{Vector{Int},Int}()
    for ja in eachindex(allowedq)
      allowedq_dic[allowedq[ja]]=ja
    end
    
    num_spin=2
    Nbs=Int(Nband/num_spin)
  single_Ham=zeros(ComplexF64,length(allowedq),Nbs,num_spin,length(allowedq),Nbs,num_spin)
  single_MoirePo=zeros(ComplexF64,length(allowedq),Nbs,num_spin,length(allowedq),Nbs,num_spin)
  for ja in eachindex(allowedq), si in 1:num_spin
    single_Ham[ja,:,si,ja,:,si]=(nop_HF_eigenvector[ja]'*nop_single_Ham[ja]*nop_HF_eigenvector[ja])[1:Nbs,1:Nbs]
    single_MoirePo[ja,:,si,ja,:,si]=(nop_HF_eigenvector[ja]'*nop_single_MoirePo[ja]*nop_HF_eigenvector[ja])[1:Nbs,1:Nbs]
  end

  single_Ham=reshape(single_Ham,length(allowedq)*Nband,length(allowedq)*Nband)
  single_MoirePo=reshape(single_MoirePo,length(allowedq)*Nband,length(allowedq)*Nband)

    
    wave_diff=Vector{Int}[]
     cutoff=Int(round.(max(qcutoff*norm(b1)/norm(b1),qcutoff*norm(b1)/norm(b2))+1))*8
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
      if norm(ja*b1+jb*b2)<qcutoff*norm(b1)
        push!(wave_diff,ja*b1T+jb*b2T)
      end
    end
    

    AA=randn(length(allowedq)*Nband,length(allowedq)*Nband)+im*randn(length(allowedq)*Nband,length(allowedq)*Nband)
    initial_DensityMatrix=AA+AA'


  return allowedq_dic,single_Ham,single_MoirePo,wave_diff,initial_DensityMatrix
        
end


function iteration_loop(initial_DensityMatrix::Matrix{ComplexF64},
                       allowedq::Vector{Vector{Int64}},allowedq_dic::Dict{Vector{Int}},T1::Vector{Float64},T2::Vector{Float64},
                       Nq::Int,wave::Vector{Vector{Int64}},wave_diff::Vector{Vector{Int64}},single_Ham::Matrix{ComplexF64},
                       single_MoirePo::Matrix{ComplexF64},constq::Float64,ϵr::Float64,formfactors::Array{Matrix{ComplexF64}},filling::Int,Area::Float64)
    eout=1.0
    itcount=0
 
     Coulomb_element=zeros(Float64,length(allowedq),length(wave_diff))
      for ja in eachindex(allowedq), jb in eachindex(wave_diff)
         Coulomb_element[ja,jb]=Coulomb(allowedq[ja]+wave_diff[jb],T1,T2)/ϵr+constq
      end
    
    HF_eigenvalue=zeros(Float64,Nq^2*Nband)
    HF_eigenvector=zeros(ComplexF64,Nq^2*Nband,Nq^2*Nband)
    HartreeMatrix=zeros(ComplexF64,Nq^2*Nband,Nq^2*Nband)
    FockMatrix=zeros(ComplexF64,Nq^2*Nband,Nq^2*Nband)
  
    sum_idx,diff_idx=build_k_loopup(allowedq,allowedq_dic,Nq)
    DIIS_size=5
    DIIS_input_DensityMatrix=Vector{Matrix{ComplexF64}}(undef,DIIS_size)
    DIIS_input_DeltaMatrix=Vector{Matrix{ComplexF64}}(undef,DIIS_size)
    input_DensityMatrix=initial_DensityMatrix
    bad_count=0
    energy=0.0
    energy_change=0.0
  

    while (eout>1*10^(-18)) || (bad_count<DIIS_size) || (abs(energy_change)>1*10^(-10))
 
      if eout<1*10^(-18)
       bad_count+=1
      else
        bad_count=0
      end
      
      tic=time()

      if (itcount>100 && abs(eout)>10^(-2)) || (itcount>30 && abs(eout)<10^(-10))
      
        dmk=implement_DIIS(DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,DIIS_size)
        if dmk==0
            itcount=0
             A=randn(Nq^2*Nband,Nq^2*Nband)+im*randn(Nq^2*Nband,Nq^2*Nband)
            dmk=(A+A')*0.01
          
        end
      

       eout,energy_change,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,DIIS_size)+1],HF_eigenvalue,HF_eigenvector,energy,HartreeMatrix,FockMatrix=Construct_DensityMatrix(allowedq,
                                                                                                                                                                                  Nq,Nband,wave_diff,
                                                                                                                                                                                dmk,single_Ham,
                                                                                                                                                                                single_MoirePo,Coulomb_element,formfactors,
                                                                                                                                                                                energy,filling,Area,sum_idx,diff_idx)
       
        DIIS_input_DensityMatrix[mod(itcount,DIIS_size)+1]=dmk
        input_DensityMatrix=output_DensityMatrix
        println("using DIIS")
       
      else
  
      

        eout,energy_change,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,DIIS_size)+1],HF_eigenvalue,HF_eigenvector,energy,HartreeMatrix,FockMatrix=Construct_DensityMatrix(allowedq,
                                                                                                                                                                                  Nq,Nband,wave_diff,
                                                                                                                                                                                input_DensityMatrix,single_Ham,
                                                                                                                                                                                single_MoirePo,Coulomb_element,formfactors,
                                                                                                                                                                                energy,filling,Area,sum_idx,diff_idx)
          DIIS_input_DensityMatrix[mod(itcount,DIIS_size)+1]=input_DensityMatrix
        input_DensityMatrix=output_DensityMatrix
        
       
         
     

      end

    



      itcount+=1
     
      toc=time()
      println(toc-tic,"eout=$eout","energy_change=$energy_change","itcount=$itcount")
      flush(stdout)
     
    
  end







  return DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,eout,HartreeMatrix,FockMatrix

end







function implement_DIIS(DIIS_input_projector::Vector{Matrix{ComplexF64}},DIIS_input_DeltaMatrix::Vector{Matrix{ComplexF64}},DIIS_size::Int64)


    Bmatrix=zeros(Float64,DIIS_size+1,DIIS_size+1)
      for ja in 1:DIIS_size
       Bmatrix[ja,DIIS_size+1]=1
       Bmatrix[DIIS_size+1,ja]=1
      end
  
       for ja in 1:DIIS_size,jb in 1:DIIS_size
   
             Bmatrix[ja,jb]+=real(dot(DIIS_input_DeltaMatrix[ja],DIIS_input_DeltaMatrix[jb]))
          
      end

      inB=safe_inverse(Bmatrix)
      if inB≠0
         onh=zeros(Float64,DIIS_size+1)
         onh[end]=1
         coeff=inB* onh
        
        dmk = zero(DIIS_input_projector[1])           # alloc once
        @inbounds for ja in 1:DIIS_size
            BLAS.axpy!(coeff[ja], DIIS_input_projector[ja], dmk)  # dmk += coeff[ja] * projector[ja]
       
            BLAS.axpy!(coeff[ja], DIIS_input_DeltaMatrix[ja], dmk)  # dmk += coeff[ja] * projector[ja]
        end
      return dmk
       
      else
        return 0
      end
end


function safe_inverse(A)
  try
      return inv(A)  # Attempt to compute inverse
  catch e
      if isa(e, SingularException)
          println("Matrix is singular, doing randomstart again.")
          return pinv(A, 0.1)  # Use pseudoinverse as an alternative
      else
          rethrow(e)  # If another error occurs, propagate it
      end
  end
end