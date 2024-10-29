
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

function sample_value(uD::Float64, numsample::Int,θ::Float64,rad::Float64)
   



    ac=0.246
    R1=ac*[1,0]
    R2=ac*[1/2,√3/2]
    G1=2π/ac*[1,-1/√3]
    G2=2π/ac*[0,2/√3]


    ϵ=0.2504/ac-1
    Rθ=[cos(θ) -sin(θ);sin(θ) cos(θ)]
    g1=G1-(1+ϵ)^(-1)*Rθ*G1
    g2=G2-(1+ϵ)^(-1)*Rθ*G2
    

  
    NL=5
    KGr=4π/(3*ac)*[1,0]
    gm=norm(g1)
    am=4π/(gm*√3)
    ns=1/(√3/2*am)^2
    
  
    gross_valuesset=[]
    for ja in 1:numsample
      ρ=rand()*gm*rad
      ang=rand()*2*π
      kvec=KGr+ρ*[cos(ang),sin(ang)]
      Hamiltonian=get_RNGham(kvec,NL,uD)
      FFF=real.(eigen(Hamiltonian).values)
      
        push!(gross_valuesset,FFF[NL+1])
    
    end
    
    
    valuesset=Float64[]
    bandmin=sort(gross_valuesset)[1]
    for ja in eachindex(gross_valuesset)
       if gross_valuesset[ja]<bandmin+15
          push!(valuesset,gross_valuesset[ja])
       end
    end
    

    return valuesset,ns,gm


end

function process_data(valuesset::Vector{Float64},ns::Float64,numsample::Int,rad::Float64,DOS_n_binnum::Int64,DOS_E_binnum::Int64,Density_start::Float64,Density_end::Float64,gm::Float64)
    h=fit(Histogram, valuesset, nbins=DOS_E_binnum) 
    bin_edges = collect(h.edges[1])
    Nstates = h.weights
    bin_centers =collect(0.5* (bin_edges[1:end-1] + bin_edges[2:end]))
    NN=numsample*(bin_centers[2]-bin_centers[1])
    Nstates=Nstates/(NN)*gm^2/(4π)*rad^2
  
    nE=zeros(Float64,length(bin_centers))
 for ja in eachindex(nE)
   nE[ja]=sum(Nstates[1:ja])*(bin_centers[2]-bin_centers[1])
 end

    

 nE_edges = collect(range(start=Density_start*ns, stop=Density_end*ns, length=DOS_n_binnum+1))
      
 nE_bin_means = Float64[]
 nE_bin_centers =(nE_edges[1:end-1] + nE_edges[2:end]) / 2

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