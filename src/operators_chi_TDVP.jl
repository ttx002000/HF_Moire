using LinearAlgebra
using Arpack
using Combinatorics
using Random
using CSV,DataFrames

function big_func(powerindex::Int,bigQindex::Int,file_pos::Int)
    V0=0.0
    ϕ=0.0
    Nq=27
    scale=0.6
    if powerindex==1
      powerrange=collect(0:1:8)
    elseif powerindex==2
      powerrange=[0]
    elseif powerindex==3
      powerrange=[1]
    elseif powerindex==4
         powerrange=[3]
    elseif powerindex==5
         powerrange=[4]
    elseif powerindex==6
         powerrange=[5]
    elseif powerindex==7
         powerrange=[6]
    elseif powerindex==8
         powerrange=[7]
    elseif powerindex==9
         powerrange=[8]
    elseif powerindex==10
    powerrange=[2]


    end

    NL=5

    gcutoff=4.01

    am=4*π/(√3*scale);
    
    aGr=0.246
    t0=3100
    t1=380
    vf=√3/2*aGr*t0
    
    b1=scale*[0,1]
    b2=scale*[√3/2,-1/2]


    a1m=am*[1/2,√3/2]
    a2m=am*[1,0]


    
    lb=(√3/2*norm(a1m)^2/(2*π))^(1/2)
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
    cutoffstandard=gcutoff*scale
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        gtest=ja*b1+jb*b2;
        if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
            push!(wave,ja*b1T+jb*b2T)
        end
    end

    dimension=length(wave)


   
   
     scratch_dir = ENV["SCRATCH"]
    d1=Matrix(CSV.read(joinpath(scratch_dir, "chi_TDVP/data_output$(file_pos)/data_input/Re_LL_x8.3y7.5.csv"),DataFrame,header=:false))
    d2=Matrix(CSV.read(joinpath(scratch_dir, "chi_TDVP/data_output$(file_pos)/data_input/Im_LL_x8.3y7.5.csv"),DataFrame,header=:false))
     

    if size(d1)[1]≠Nq^2 || size(d2)[1]≠Nq^2
      error("This is wrong")
    end
    LL0_vals_p3=zeros(ComplexF64,length(allowedq))
    LL1_vals_p3=zeros(ComplexF64,length(allowedq))
    LL2_vals_p3=zeros(ComplexF64,length(allowedq))
    exp_p3=[exp(im*(dot([T1 T2]*(allowedq[ja]+wave[jb]),[8.3,7.5])))*get_spinor(NL,[T1 T2]*(allowedq[ja]+wave[jb]))[1] for ja in eachindex(allowedq), jb in eachindex(wave)]

    for ja in eachindex(d1[:,1])
      ppos=findfirst(item->item==Int.([d1[ja,1],d1[ja,2]]),allowedq)
      if Int.([d1[ja,1],d1[ja,2]])==Int.([d2[ja,1],d2[ja,2]])
        LL0_vals_p3[ppos]=d1[ja,3]+im*d2[ja,3]
        LL1_vals_p3[ppos]=d1[ja,4]+im*d2[ja,4]
        LL2_vals_p3[ppos]=d1[ja,5]+im*d2[ja,5]
      end

    end




   
    cnum=1
    Abz=√3/2*scale^2
    eigenvector_try=[zeros(ComplexF64,dimension) for _ in 1:Nq^2]
    k_dependent_const=zeros(ComplexF64,length(allowedq))
    for ja in 1:Nq^2, jb in eachindex(wave)
                k=allowedq[ja]
                g=wave[jb]
                kvec=k[1]*T1+k[2]*T2
                gvec=g[1]*T1+g[2]*T2
            
                phase1=exp(im*π/Abz*cnum*(kvec[1]*gvec[2]-kvec[2]*gvec[1]))
                nvec=Int.(round.(inv([b1T b2T])*g))
                phase2=exp(im*cnum*π*nvec[1]*nvec[2])
                factor1=(kvec[1]+gvec[1]+im*(kvec[2]+gvec[2]))^0
                eigenvector_try[ja][jb]=factor1/(exp(π*cnum/(2*Abz)*norm(kvec+gvec)^2))/get_spinor_norm(NL,kvec+gvec)*phase1*phase2
    end


    for ja in eachindex(allowedq)
      k_dependent_const[ja]=LL0_vals_p3[ja]*conj(LL0_vals_p3[1])/sum(eigenvector_try[ja].*exp_p3[ja,:])
    end

    for ja in eachindex(allowedq)
      eigenvector_try[ja]=eigenvector_try[ja]*k_dependent_const[ja]
    end


    


      
               println(bigQindex)
               flush(stdout)
                bigQ=allowedq[bigQindex]
                bigQvec=[T1 T2]*bigQ
                shiftvec=[-bigQvec[2],bigQvec[1]]*(norm(a1m)^2*√3/2)/(2*π)

            



                eigenvector_excited=[zeros(ComplexF64,dimension) for _ in eachindex(powerrange), _ in 1:Nq^2]
                eigenvector_excited_zeropower=[zeros(ComplexF64,dimension) for _ in 1:Nq^2]
                k_dependent_const_excited=zeros(ComplexF64,length(allowedq))

                for jp in eachindex(powerrange)
                  power=powerrange[jp]
                      for ja in 1:Nq^2, jb in eachindex(wave)
                                  k=allowedq[ja]
                                  g=wave[jb]
                                  kvec=k[1]*T1+k[2]*T2
                                  gvec=g[1]*T1+g[2]*T2
                              
                                  phase1=exp(im*π/Abz*cnum*(kvec[1]*gvec[2]-kvec[2]*gvec[1]))
                                  nvec=Int.(round.(inv([b1T b2T])*g))
                                  phase2=exp(im*cnum*π*nvec[1]*nvec[2])
                                  factor1=(-(kvec[1]+gvec[1]+im*(kvec[2]+gvec[2]))/√2)^power*lb^power/(√factorial(power))
                              
                                  phase3=exp(-im*dot(shiftvec,kvec+gvec))
                                  eigenvector_excited[jp,ja][jb]=factor1/(exp(π*cnum/(2*Abz)*norm(kvec+gvec)^2))/get_spinor_norm(NL,kvec+gvec)*phase1*phase2*phase3
                      end

                      
                end



                for ja in 1:Nq^2, jb in eachindex(wave)
                      k=allowedq[ja]
                      g=wave[jb]
                      kvec=k[1]*T1+k[2]*T2
                      gvec=g[1]*T1+g[2]*T2
                              
                      phase1=exp(im*π/Abz*cnum*(kvec[1]*gvec[2]-kvec[2]*gvec[1]))
                      nvec=Int.(round.(inv([b1T b2T])*g))
                      phase2=exp(im*cnum*π*nvec[1]*nvec[2])
                      factor1=(-(kvec[1]+gvec[1]+im*(kvec[2]+gvec[2]))/√2)^0
                              
                      phase3=exp(-im*dot(shiftvec,kvec+gvec))
                      eigenvector_excited_zeropower[ja][jb]=factor1/(exp(π*cnum/(2*Abz)*norm(kvec+gvec)^2))/get_spinor_norm(NL,kvec+gvec)*phase1*phase2*phase3
                end

                      

                for ja in eachindex(allowedq)
                  minusbigQpos=findfirst(item->item==mod.(-bigQ,Nq),allowedq)
                  k_minusbigQpos=findfirst(item->item==mod.(allowedq[ja]-bigQ,Nq),allowedq)

                  k_dependent_const_excited[ja]=LL0_vals_p3[k_minusbigQpos]*conj(LL0_vals_p3[minusbigQpos])/sum(eigenvector_excited_zeropower[ja].*exp_p3[ja,:])
                end

                for ja in eachindex(allowedq), jp in eachindex(powerrange)
                  eigenvector_excited[jp,ja]=eigenvector_excited[jp,ja]*k_dependent_const_excited[ja]
                end



                minus_bigQ=mod.(-allowedq[bigQindex],Nq)
                minus_bigQvec=[T1 T2]*minus_bigQ
                minus_shiftvec=[-minus_bigQvec[2],minus_bigQvec[1]]*(norm(a1m)^2*√3/2)/(2*π)


                eigenvector_excited_mq=[zeros(ComplexF64,dimension) for _ in eachindex(powerrange), _ in 1:Nq^2]
                eigenvector_excited_zeropower_mq=[zeros(ComplexF64,dimension) for _ in 1:Nq^2]
                k_dependent_const_excited_mq=zeros(ComplexF64,length(allowedq))

                for jp in eachindex(powerrange)
                  power=powerrange[jp]
                      for ja in 1:Nq^2, jb in eachindex(wave)
                                  k=allowedq[ja]
                                  g=wave[jb]
                                  kvec=k[1]*T1+k[2]*T2
                                  gvec=g[1]*T1+g[2]*T2
                              
                                  phase1=exp(im*π/Abz*cnum*(kvec[1]*gvec[2]-kvec[2]*gvec[1]))
                                  nvec=Int.(round.(inv([b1T b2T])*g))
                                  phase2=exp(im*cnum*π*nvec[1]*nvec[2])
                                  factor1=(-(kvec[1]+gvec[1]+im*(kvec[2]+gvec[2]))/√2)^power*lb^power/(√factorial(power))
                              
                                  phase3=exp(-im*dot(minus_shiftvec,kvec+gvec))
                                  eigenvector_excited_mq[jp,ja][jb]=factor1/(exp(π*cnum/(2*Abz)*norm(kvec+gvec)^2))/get_spinor_norm(NL,kvec+gvec)*phase1*phase2*phase3
                      end

                      
                end



                for ja in 1:Nq^2, jb in eachindex(wave)
                      k=allowedq[ja]
                      g=wave[jb]
                      kvec=k[1]*T1+k[2]*T2
                      gvec=g[1]*T1+g[2]*T2
                              
                      phase1=exp(im*π/Abz*cnum*(kvec[1]*gvec[2]-kvec[2]*gvec[1]))
                      nvec=Int.(round.(inv([b1T b2T])*g))
                      phase2=exp(im*cnum*π*nvec[1]*nvec[2])
                      factor1=(-(kvec[1]+gvec[1]+im*(kvec[2]+gvec[2]))/√2)^0
                              
                      phase3=exp(-im*dot(minus_shiftvec,kvec+gvec))
                      eigenvector_excited_zeropower_mq[ja][jb]=factor1/(exp(π*cnum/(2*Abz)*norm(kvec+gvec)^2))/get_spinor_norm(NL,kvec+gvec)*phase1*phase2*phase3
                end

                      

                for ja in eachindex(allowedq)
                  bigQpos=findfirst(item->item==mod.(bigQ,Nq),allowedq)
                  kplusq_pos=findfirst(item->item==mod.(allowedq[ja]+bigQ,Nq),allowedq)

                  k_dependent_const_excited_mq[ja]=LL0_vals_p3[kplusq_pos]*conj(LL0_vals_p3[bigQpos])/sum(eigenvector_excited_zeropower_mq[ja].*exp_p3[ja,:])
                end

                for ja in eachindex(allowedq), jp in eachindex(powerrange)
                  eigenvector_excited_mq[jp,ja]=eigenvector_excited_mq[jp,ja]*k_dependent_const_excited_mq[ja]
                end

                kinetic_mm=[zeros(Float64,length(wave)) for _ in eachindex(allowedq)]
                for ja in eachindex(allowedq), jb in eachindex(wave)
                kinetic_mm[ja][jb]=norm([T1 T2]*(allowedq[ja]+wave[jb]))^2*200
                end



                MMmatrix=zeros(ComplexF64,length(powerrange),length(powerrange)) #This is Gmatrix
               Threads.@threads for mindex in eachindex(powerrange)
                   for nindex in eachindex(powerrange)
                  
                      for emsite in eachindex(allowedq)
                              emk_plusq=mod.(allowedq[emsite]+bigQ,Nq)
                              emk_plusq_pos=findfirst(item->item==emk_plusq,allowedq)
                              v1_m=eigenvector_try[emk_plusq_pos]
                              v2_m=eigenvector_excited[mindex,emk_plusq_pos]
                              v1_n=eigenvector_try[emk_plusq_pos]
                              v2_n=eigenvector_excited[nindex,emk_plusq_pos]

                              dd1=det([v1_m'*v1_n v1_m'*v2_n; v2_m'*v1_n v2_m'*v2_n])
                              dd1=dd1/(norm(eigenvector_try[emsite]))^2
                              dd1=dd1/(norm(eigenvector_try[emk_plusq_pos]))^2
                              MMmatrix[mindex,nindex]+=dd1
                      end
                    end
                end
                   println("finishMM")
                   flush(stdout)


            MHmatrix=zeros(ComplexF64,length(powerrange),length(powerrange))# This is Amatrix
            Threads.@threads for mindex in eachindex(powerrange)
              for nindex in eachindex(powerrange) 

                      for emsite in eachindex(allowedq)
                            emk_plusq=mod.(allowedq[emsite]+bigQ,Nq)
                            emk_plusq_pos=findfirst(item->item==emk_plusq,allowedq)
                            

                              v1_m=eigenvector_try[emk_plusq_pos]
                              v2_m=eigenvector_excited[mindex,emk_plusq_pos]
                              v1_n=eigenvector_try[emk_plusq_pos]
                              v2_n=eigenvector_excited[nindex,emk_plusq_pos]
                              overlap_determinant=det([v1_m'*v1_n v1_m'*v2_n; v2_m'*v1_n v2_m'*v2_n])
                              nn1=(norm(eigenvector_try[emsite])^2)*(norm(eigenvector_try[emk_plusq_pos])^2)
                              overlap_determinant=overlap_determinant/nn1

                              overlap_determinant_revise1=det([v1_m'*(kinetic_mm[emk_plusq_pos].*v1_n) v1_m'*v2_n; v2_m'*(kinetic_mm[emk_plusq_pos].*v1_n) v2_m'*v2_n])
                              overlap_determinant_revise1=overlap_determinant_revise1/nn1
                              
                              overlap_determinant_revise2=det([v1_m'*v1_n v1_m'*(kinetic_mm[emk_plusq_pos].*v2_n); v2_m'*v1_n v2_m'*(kinetic_mm[emk_plusq_pos].*v2_n)])
                              overlap_determinant_revise2=overlap_determinant_revise2/nn1

                              overlap_determinant_revise3=eigenvector_try[emsite]'*(kinetic_mm[emsite].*eigenvector_try[emsite])/norm(eigenvector_try[emsite])^2
                              overlap_determinant_revise3=overlap_determinant_revise3*overlap_determinant

                              overlap_determinant_revise4=eigenvector_try[emk_plusq_pos]'*(kinetic_mm[emk_plusq_pos].*eigenvector_try[emk_plusq_pos])/norm(eigenvector_try[emk_plusq_pos])^2
                              overlap_determinant_revise4=overlap_determinant_revise4*overlap_determinant
                              


                              v3=overlap_determinant_revise1+overlap_determinant_revise2-overlap_determinant_revise3-overlap_determinant_revise4
                        



                            MHmatrix[mindex,nindex]+=v3

                            

                      end
                    end
                 end
                   println("finishMH")
                   flush(stdout)




            Bmatrix=zeros(ComplexF64,length(powerrange),length(powerrange))# This is Bmatrix

            Threads.@threads for mindex in eachindex(powerrange)
              for nindex in eachindex(powerrange) 

                  for emsite in eachindex(allowedq)
                        emk_plusq=mod.(allowedq[emsite]+bigQ,Nq)
                        emk_plusq_pos=findfirst(item->item==emk_plusq,allowedq)

                        v3=(eigenvector_excited[mindex,emk_plusq_pos]'*eigenvector_try[emk_plusq_pos])*(eigenvector_excited_mq[nindex,emsite]'*eigenvector_try[emsite])
                        v3=v3/(norm(eigenvector_try[emsite])^2*norm(eigenvector_try[emk_plusq_pos])^2)
                        v3=v3*(eigenvector_try[emsite]'*(kinetic_mm[emsite].*eigenvector_try[emsite])/norm(eigenvector_try[emsite])^2+eigenvector_try[emk_plusq_pos]'*(kinetic_mm[emk_plusq_pos].*eigenvector_try[emk_plusq_pos])/norm(eigenvector_try[emk_plusq_pos])^2)

                        v4=(eigenvector_excited[mindex,emk_plusq_pos]'*(kinetic_mm[emk_plusq_pos].*eigenvector_try[emk_plusq_pos]))*(eigenvector_excited_mq[nindex,emsite]'*eigenvector_try[emsite])/(norm(eigenvector_try[emsite])^2*norm(eigenvector_try[emk_plusq_pos])^2)
                        v5=(eigenvector_excited[mindex,emk_plusq_pos]'*eigenvector_try[emk_plusq_pos])*(eigenvector_excited_mq[nindex,emsite]'*(kinetic_mm[emsite].*eigenvector_try[emsite]))/(norm(eigenvector_try[emsite])^2*norm(eigenvector_try[emk_plusq_pos])^2)
                        

                      Bmatrix[mindex,nindex]+=v3-v4-v5
                  end
                end
            end
    

             println("finishB")
                   flush(stdout)
            
  


      


        


   return MHmatrix, Bmatrix,MMmatrix


end



function get_spinor(NL::Int,q::Vector{Float64})
    qh=q[1]+im*q[2]
    aGr=0.246
    t0=3100
    t1=380
    vf=√3/2*aGr*t0
    spinor=zeros(ComplexF64,NL)
    for ja in 1:NL
      spinor[ja]=qh^(ja-1)*(vf/t1)^(ja-1)
    end
  return spinor/norm(spinor)
end

function get_spinor_norm(NL::Int,q::Vector{Float64})
    qh=q[1]+im*q[2]
    aGr=0.246
    t0=3100
    t1=380
    vf=√3/2*aGr*t0
    spinor=zeros(ComplexF64,NL)
    for ja in 1:NL
      spinor[ja]=qh^(ja-1)*(vf/t1)^(ja-1)
    end
  return 1/norm(spinor)
end