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
 
  input_ds_reshaped=reshape(input_DensityMatrix,Nq^2,Nband,Nq^2,Nband)
   
  HartreeMatrix=zeros(ComplexF64,Nq^2,Nband,Nq^2,Nband)
  FockMatrix=zeros(ComplexF64,Nq^2,Nband,Nq^2,Nband)
   nt = Threads.maxthreadid()
  Xbuf = [Matrix{ComplexF64}(undef, Nband, Nband) for _ in 1:nt]
  Ybuf = [Matrix{ComplexF64}(undef, Nband, Nband) for _ in 1:nt]

  tic=time()
  Threads.@threads for k2 in eachindex(allowedq)
    tid=Threads.threadid()
    X=Xbuf[tid]
    Y=Ybuf[tid]
    for k3 in eachindex(allowedq)
    Fk_local=@view  FockMatrix[k2,:,k3,:] 
       for k1 in eachindex(allowedq)
          #qindex=allowedq_dic[mod.(allowedq[k1]-allowedq[k3],Nq)]
          #k4=allowedq_dic[mod.(allowedq[k1]+allowedq[k2]-allowedq[k3],Nq)]
           qindex=diff_idx[k1,k3]
           k4=sum_idx[qindex,k2]
           ρ=@view input_ds_reshaped[k4,:,k1,:]
           for qg in eachindex(wave_diff)
             #Fk_local.+=Coulomb_element[qindex,qg]*(formf[k2,qindex,qg]*input_ds_reshaped[k4,:,k1,:]*formf[k3,qindex,qg]')
             α=Coulomb_element[qindex,qg]  
             A=formf[k2,qindex,qg]   
             B=formf[k3,qindex,qg]
             mul!(X,ρ,B')
             mul!(Y,A,X)
             Fk_local.+=α.*Y
            end
         
       end
      end
  end
  toc=time()
  println(toc-tic,"finisheFock")
  
  tic=time()

  Threads.@threads for k2 in eachindex(allowedq)
    for k4 in eachindex(allowedq)
    HT_local=@view HartreeMatrix[k2,:,k4,:]
    #qindex=allowedq_dic[mod.(allowedq[k4]-allowedq[k2],Nq)]
     qindex=diff_idx[k4,k2]
     for qg in eachindex(wave_diff)
        htdensity=zero(ComplexF64)
        for k1 in eachindex(allowedq)
          #k3=allowedq_dic[mod.(allowedq[k1]+allowedq[k2]-allowedq[k4],Nq)]
          k3=diff_idx[k1,qindex]
          #htdensity+=tr(input_ds_reshaped[k3,:,k1,:]*formf[k3,qindex,qg]')
          htdensity+=dot(formf[k3,qindex,qg],input_ds_reshaped[k3,:,k1,:])
        end
        β=Coulomb_element[qindex,qg]*htdensity
        HT_local.+=β.*formf[k2,qindex,qg]
      end
   end
  end

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


   

  eout=real(tr(DeltaMatrix'*DeltaMatrix))/Nq^2
  
  


  ss=single_MoirePo+single_Ham+0.5*HartreeMatrix-0.5*FockMatrix
  energy=real(tr(ss*input_DensityMatrix))


   energy_change=real(energy-energy_input)



 return  eout,energy_change,output_DensityMatrix,DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,HartreeMatrix,FockMatrix
end


function get_ff(nop_HF_eigenvector::Vector{Matrix{ComplexF64}},NL::Int,Nband::Int,wave::Vector{Vector{Int}}
                ,wave_diff::Vector{Vector{Int}},spinor_set::Matrix{Vector{ComplexF64}},allowedq::Vector{Vector{Int}},allowedq_dic::Dict{Vector{Int}})
  num_spin=2
   
  form_factors=[zeros(ComplexF64,Nband,Nband) for _ in eachindex(allowedq), _ in eachindex(allowedq), _ in eachindex(wave_diff)]

  eigenvector_PWbasis=[zeros(ComplexF64,2*NL,num_spin,length(wave),Nband) for _ in eachindex(allowedq)]
  diff_eigenvector_threads=[zeros(ComplexF64,2*NL,num_spin,length(wave),Nband,length(allowedq)) for _ in eachindex(wave_diff)]
  
  
 
  
  nop_HF_eigenvector_reshaped=[reshape(nop_HF_eigenvector[ja],num_spin,length(wave),num_spin*length(wave)) for ja in eachindex(allowedq)]

  for ja in eachindex(allowedq)
   
    for jc in eachindex(wave), si in 1:num_spin,bi in 1:Nband
      eigenvector_PWbasis[ja][:,si,jc,bi]=nop_HF_eigenvector_reshaped[ja][si,jc,bi]*spinor_set[ja,jc]
    end
    
  end


   Threads.@threads for jd in eachindex(wave_diff)
      for ja in eachindex(wave)
            pos=findfirst(item->item==wave[ja]-wave_diff[jd],wave)
            if pos≠nothing
              for jk in eachindex(allowedq)
                 diff_eigenvector_threads[jd][:,:,ja,:,jk]=eigenvector_PWbasis[jk][:,:,pos,:]
              end
            end
      end
   end
   


   eigenvector_PWbasis=[reshape(eigenvector_PWbasis[ja],2*NL*num_spin*length(wave),Nband) for ja in eachindex(allowedq)]
   diff_eigenvector_threads=[reshape(diff_eigenvector_threads[ja],2*NL*num_spin*length(wave),Nband,length(allowedq)) for ja in eachindex(wave_diff)] 




   Threads.@threads for ja in 1:Nq^2
    ff_local=@view  form_factors[ja,:,:]
    for jb in 1:Nq^2
     
        meshk1plusq=mod.(allowedq[ja]+allowedq[jb],Nq)
        k2pos=allowedq_dic[meshk1plusq]
       for  jc in eachindex(wave_diff)
           gk1plusq=wave_diff[jc]+allowedq[ja]+allowedq[jb]-meshk1plusq
           gk1plusq_pos=findfirst(item->item==gk1plusq,wave_diff)
           
           if gk1plusq_pos≠nothing
                ff_local[jb,jc]=(diff_eigenvector_threads[gk1plusq_pos][:,:,ja])'*eigenvector_PWbasis[k2pos][:,:]
           end
       end
     end
   end



   return form_factors
end

function initial_process(allowedq::Vector{Vecotr{Int}},nop_single_Ham::Vector{Matrix{ComplexF64}},
                         nop_single_MoirePo::Vector{Matrix{ComplexF64}},nop_HF_eigenvector::Vector{Matrix{ComplexF64}},
                         Nband::Int64,qcutoff::Float64,b1::Vector{Float64},b2::Vector{Float64},b1T::Vector{Int},b2T::Vector{Int})
  allowedq_dic=Dict{Vector{Int},Int}()
    for ja in eachindex(allowedq)
      allowedq_dic[allowedq[ja]]=ja
    end
    
  single_Ham=zeros(ComplexF64,length(allowedq),Nband,length(allowedq),Nband)
  single_MoirePo=zeros(ComplexF64,length(allowedq),Nband,length(allowedq),Nband)
  for ja in eachindex(allowedq)
    single_Ham[ja,:,ja,:]=(nop_HF_eigenvector[ja]'*nop_single_Ham[ja]*nop_HF_eigenvector[ja])[1:Nband,1:Nband]
    single_MoirePo[ja,:,ja,:]=(nop_HF_eigenvector[ja]'*nop_single_MoirePo[ja]*nop_HF_eigenvector[ja])[1:Nband,1:Nband]
  end

  single_Ham=reshape(single_Ham,length(allowedq)*Nband,length(allowedq)*Nband)
  single_MoirePo=reshape(single_MoirePo,length(allowedq)*Nband,length(allowedq)*Nband)

    
    wave_diff=Vector{Int}[]
    cutoff=20
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
      if norm(ja*b1+jb*b2)<qcutoff*norm(b1)
        push!(wave_diff,ja*b1T+jb*b2T)
      end
    end


  return allowedq_dic,single_Ham,single_MoirePo,wave_diff
        
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
    DIIS_input_DensityMatrix=Vector{Matrix{ComplexF64}}(undef,3)
    DIIS_input_DeltaMatrix=Vector{Matrix{ComplexF64}}(undef,3)
    input_DensityMatrix=initial_DensityMatrix
    bad_count=0
    energy=0.0
    energy_change=0.0
  

    while (eout>1*10^(-22)) || (bad_count<4) || (energy_change>1*10^(-10))
 
      if eout<1*10^(-22)
       bad_count+=1
      end
      
      tic=time()

      if (itcount>100 && abs(eout)>10^(-2)) || (itcount>30 && abs(eout)<10^(-10))
      
        dmk=implement_DIIS(DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix)
        if dmk==0
            itcount=0
             A=randn(Nq^2*Nband,Nq^2*Nband)+im*randn(Nq^2*Nband,Nq^2*Nband)
            dmk=(A+A')*0.01
          
        end
      

       eout,energy_change,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue,HF_eigenvector,energy,HartreeMatrix,FockMatrix=Construct_DensityMatrix(allowedq,
                                                                                                                                                                                  Nq,Nband,wave_diff,
                                                                                                                                                                                dmk,single_Ham,
                                                                                                                                                                                single_MoirePo,Coulomb_element,formfactors,
                                                                                                                                                                                energy,filling,Area,sum_idx,diff_idx)
       
        DIIS_input_DensityMatrix[mod(itcount,3)+1]=dmk
        input_DensityMatrix=output_DensityMatrix
        println("using DIIS")
       
      else
  
      

        eout,energy_change,output_DensityMatrix,DIIS_input_DeltaMatrix[mod(itcount,3)+1],HF_eigenvalue,HF_eigenvector,energy,HartreeMatrix,FockMatrix=Construct_DensityMatrix(allowedq,
                                                                                                                                                                                  Nq,Nband,wave_diff,
                                                                                                                                                                                input_DensityMatrix,single_Ham,
                                                                                                                                                                                single_MoirePo,Coulomb_element,formfactors,
                                                                                                                                                                                energy,filling,Area,sum_idx,diff_idx)
        input_DensityMatrix=output_DensityMatrix
        
       
         
     

      end

    



      itcount+=1
     
      toc=time()
      println(toc-tic,"eout=$eout","energy_change=$energy_change","itcount=$itcount")
      flush(stdout)
     
    
  end







  return DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,eout,HartreeMatrix,FockMatrix

end







function implement_DIIS(DIIS_input_projector::Vector{Matrix{ComplexF64}},DIIS_input_DeltaMatrix::Vector{Matrix{ComplexF64}})



      Bmatrix=zeros(ComplexF64,4,4)
      for ja in 1:3
       Bmatrix[ja,4]=1
       Bmatrix[4,ja]=1
      end
  
      for ja in 1:3,jb in 1:3
          
             Bmatrix[ja,jb]+=real(tr(DIIS_input_DeltaMatrix[ja]'*DIIS_input_DeltaMatrix[jb]))
       
      end

      inB=safe_inverse(Bmatrix)
      if inB≠0
         coeff=inB*[0;0;0;1]
         dmk=coeff[1]*(DIIS_input_projector[1]+DIIS_input_DeltaMatrix[1])+coeff[2]*(DIIS_input_projector[2]+DIIS_input_DeltaMatrix[2])+coeff[3]*(DIIS_input_projector[3]+DIIS_input_DeltaMatrix[3])
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