
using LinearAlgebra,Plots,StatsBase

function get_f(k)
    delta1=1/√3*0.246*[0,1]
    delta2=1/√3*0.246*[√3/2,-1/2]
    delta3=1/√3*0.246*[-√3/2,-1/2]
 
    return exp(im*dot(k,delta1))+exp(im*dot(k,delta2))+exp(im*dot(k,delta3))
end



function get_RNGham(k::Vector{Float64},NL::Int,uD::Float64)
    Ham=zeros(ComplexF64,2*NL,2*NL)
    t0=3100
    t1=380
    t2=-21
    t3=290
    t4=141
    for layer in 1:NL-1
       Ham[2*layer-1:2*layer,2*layer+1:2*layer+2]=[t4*get_f(k) t3*conj(get_f(k));t1 t4*get_f(k)]
    end

    for layer in 1:NL-2
        Ham[2*layer-1:2*layer,2*layer+3:2*layer+4]=[0.0 t2/2;0.0 0.0]
    end

    Ham=Ham+Ham'

    for layer in 1:NL
        #Ham[2*layer-1:2*layer,2*layer-1:2*layer]=[uD*(layer+1-(NL-1)/2) -t0*get_f(k);-t0*conj(get_f(k)) uD*(layer+1-(NL-1)/2)]
        Ham[2*layer-1:2*layer,2*layer-1:2*layer]=[uD*(layer-(NL+1)/2) -t0*get_f(k);-t0*conj(get_f(k)) uD*(layer-(NL+1)/2)]
 
      end
   return Ham
end





function get_MoireHam(k::Vector{Float64},wave::Vector{Vector{Int}},NL::Int,uD::Float64,T1::Vector{Float64},T2::Vector{Float64},V1::Float64,V0::Float64,ψ::Float64,g1T::Vector{Int},g2T::Vector{Int})
    
  ω=exp(im*2π/3)
    Hamiltonian=zeros(ComplexF64,length(wave)*2*NL,length(wave)*2*NL)
    for jb in eachindex(wave)
        kvec=k+[T1 T2]*wave[jb]
        Hamiltonian[2*NL*(jb-1)+1:2*NL*jb,2*NL*(jb-1)+1:2*NL*jb]=get_RNGham(kvec,NL,uD)
    end
   
    Moire=zeros(ComplexF64,length(wave)*2*NL,length(wave)*2*NL)

    for jc in eachindex(wave)
       pos=findfirst(item->item==wave[jc]+g1T,wave)
       if pos≠nothing
          Moire[2*NL*(pos-1)+1:2*NL*(pos-1)+2,2*NL*(jc-1)+1:2*NL*(jc-1)+2]+=V1*exp(-im*ψ)*[1 1;ω ω]
       end

       pos=findfirst(item->item==wave[jc]+g2T,wave)
       if pos≠nothing
          Moire[2*NL*(pos-1)+1:2*NL*(pos-1)+2,2*NL*(jc-1)+1:2*NL*(jc-1)+2]+=V1*exp(-im*ψ)*[1 ω^2;ω^2 ω]

       end

       pos=findfirst(item->item==wave[jc]-g1T-g2T,wave)
       if pos≠nothing
        Moire[2*NL*(pos-1)+1:2*NL*(pos-1)+2,2*NL*(jc-1)+1:2*NL*(jc-1)+2]+=V1*exp(-im*ψ)*[1 ω;1 ω]
       end
   
       Moire[2*NL*(jc-1)+1:2*NL*(jc-1)+2,2*NL*(jc-1)+1:2*NL*(jc-1)+2]+=V0/2*[1 0; 0 1]
    end



    return Moire+Moire'+Hamiltonian




end





function sample_value(uD::Float64, numsample::Int,θ::Float64,Ecut::Float64)
   



    ac=0.246
    R1=ac*[1,0]
    R2=ac*[1/2,√3/2]
    G1=2π/ac*[1,-1/√3]
    G2=2π/ac*[0,2/√3]


    ϵ=0.2504/ac-1
    Rθ=[cos(θ) -sin(θ);sin(θ) cos(θ)]
    g1=G1-(1+ϵ)^(-1)*Rθ*G1
    g2=G2-(1+ϵ)^(-1)*Rθ*G2

    T1=g1 
    T2=g2

    V0=28.9
    V1=21.0
   

    ψ=-0.29
    
    g1T=Int.(round.(inv([T1 T2])*g1))
    g2T=Int.(round.(inv([T1 T2])*g2))
  
    NL=5
    KGr=4π/(3*ac)*[1,0]
    gm=norm(g1)
    am=4π/(gm*√3)
    ns=1/(√3/2*am)^2
    
     
    
   
    wave=Vector{Int64}[]
    cutoff=18
    cutoffstandard=3.1*norm(g1)
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        gtest=ja*g1+jb*g2;
        if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
             push!(wave,ja*g1T+jb*g2T)
         end
    end




    gross_valuesset=[Vector{Vector{Float64}}() for _ in 1:Threads.nthreads()]
    
    
    
    conduction_band_record=[Vector{Float64}() for _ in 1:Threads.nthreads()]
    valence_band_record=[Vector{Float64}() for _ in 1:Threads.nthreads()]
    
    Threads.@threads for ja in 1:2000
      kvec=KGr+(rand()-1/2)*g1+(rand()-1/2)*g2
      Hamiltonian=get_MoireHam(kvec,wave,NL,uD,T1,T2,V1,V0,ψ,g1T,g2T)
      FFF=real.(eigen(Hamiltonian).values)
     
      push!(gross_valuesset[Threads.threadid()],FFF)
    
      push!(conduction_band_record[Threads.threadid()],FFF[length(wave)*NL+1])
      push!(valence_band_record[Threads.threadid()],FFF[length(wave)*NL])
    end
    
    conduction_band_record=reduce(vcat,conduction_band_record)
    valence_band_record=reduce(vcat,valence_band_record)

    conduction_bandmin=sort(conduction_band_record)[1]
    conduction_bandmax=sort(conduction_band_record)[end]

    valence_bandmin=sort(valence_band_record)[1]
    valence_bandmax=sort(valence_band_record)[end]
    
   #I start over
    conduction_band_record=[Vector{Float64}() for _ in 1:Threads.nthreads()]
    valence_band_record=[Vector{Float64}() for _ in 1:Threads.nthreads()]
    
    
    Threads.@threads for ja in 1:numsample
      kvec=KGr+(rand()-1/2)*g1+(rand()-1/2)*g2
      Hamiltonian=get_MoireHam(kvec,wave,NL,uD,T1,T2,V1,V0,ψ,g1T,g2T)
      FFF=sort(real.(eigen(Hamiltonian).values))

       # Range to search
        lower_bound = valence_bandmax-Ecut-10
        upper_bound = conduction_bandmin+Ecut+10


        start_idx = searchsortedfirst(FFF, lower_bound)
        end_idx = searchsortedlast(FFF, upper_bound)
        elements_in_range = FFF[start_idx:end_idx]


      
        push!(gross_valuesset[Threads.threadid()],elements_in_range)
        push!(conduction_band_record[Threads.threadid()],FFF[length(wave)*NL+1])
        push!(valence_band_record[Threads.threadid()],FFF[length(wave)*NL])
    end
    
    gross_valuesset=reduce(vcat,reduce(vcat,gross_valuesset))
    conduction_band_record=reduce(vcat,conduction_band_record)
    valence_band_record=reduce(vcat,valence_band_record)

    valuesset=Float64[]
    conduction_bandmin=sort(conduction_band_record)[1]
    conduction_bandmax=sort(conduction_band_record)[end]

    valence_bandmin=sort(valence_band_record)[1]
    valence_bandmax=sort(valence_band_record)[end]


    for ja in eachindex(gross_valuesset)
       if valence_bandmax-Ecut<gross_valuesset[ja]<conduction_bandmin+Ecut
          push!(valuesset,gross_valuesset[ja])
       end
    end
    println("finisheddiagonalization")

    return valuesset,ns,gm,conduction_bandmin,conduction_bandmax,valence_bandmin,valence_bandmax


end

function process_data(valuesset::Vector{Float64},ns::Float64,numsample::Int,DOS_n_binnum::Int64,DOS_E_binnum::Int64,Density_start::Float64,Density_end::Float64,gm::Float64,CNP_point::Float64)
    h=fit(Histogram, valuesset, nbins=DOS_E_binnum) 
    bin_edges = collect(h.edges[1])
    Nstates = h.weights
    bin_centers =collect(0.5* (bin_edges[1:end-1] + bin_edges[2:end]))
    NN=numsample*(bin_centers[2]-bin_centers[1])
    Nstates=Nstates/(NN)*gm^2*√3/2*1/(4π^2)
  
  
    num_below_CNP=0
    for ja in eachindex(valuesset)
      if valuesset[ja]<CNP_point
         num_below_CNP+=1
      end 
    end
   n_background=num_below_CNP/numsample*gm^2*√3/2*1/(4π^2)

    nE=zeros(Float64,length(bin_centers))
 for ja in eachindex(nE)
   nE[ja]=sum(Nstates[1:ja])*(bin_centers[2]-bin_centers[1])-n_background
 end


  

    

 nE_edges = collect(range(start=Density_start*ns, stop=Density_end*ns, length=DOS_n_binnum+1))
      
 nE_bin_means = Float64[] #This will hold a list of density of states
 nE_bin_centers =(nE_edges[1:end-1] + nE_edges[2:end]) / 2 #This is the list of density of carriers

 for i in 1:length(nE_edges)-1
 
    in_bin = []
    for ja in eachindex(nE)
      if nE[ja]>=nE_edges[i]&& nE[ja]<nE_edges[i+1]
        push!(in_bin,ja)
      end
    end
  
    if length(in_bin)>=1
        push!(nE_bin_means,mean(Nstates[in_bin]))
    else
        push!(nE_bin_means, 0.0)
    end
 end

 return nE_bin_centers,nE_bin_means,bin_centers,nE,Nstates


end