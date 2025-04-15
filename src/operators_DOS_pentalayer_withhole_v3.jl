
using LinearAlgebra,Plots,StatsBase

function get_f(k)
  delta1=1/√3*0.246*[0,1]
  delta2=1/√3*0.246*[√3/2,-1/2]
  delta3=1/√3*0.246*[-√3/2,-1/2]

  return exp(im*dot(k,delta1))+exp(im*dot(k,delta2))+exp(im*dot(k,delta3))
end



function get_ABCBA(k::Vector{Float64},NL::Int,uD::Float64)

  Ham=zeros(ComplexF64,2*NL,2*NL)
  Kac=4π/(3*0.246)*[1,0]*valley
  t0=3100
  t1=380

  for layer in 1:Int((NL/2-0.5))
     Ham[2*layer-1:2*layer,2*layer+1:2*layer+2]=[0 0;t1 0]
  end

  for layer in Int((NL/2+0.5)):NL-1
      Ham[2*layer-1:2*layer,2*layer+1:2*layer+2]=[0 t1;0 0]
  end

 

  Ham=Ham+Ham'

  for layer in 1:NL
      Ham[2*layer-1:2*layer,2*layer-1:2*layer]=[uD*(layer-(NL+1)/2) -t0*get_f((k+Kac)*stacking);-t0*conj(get_f((k+Kac)*stacking)) uD*(layer-(NL+1)/2)]
  end


 
 return Ham
end





function get_perturb_Ham(perturb::Int)
    
  NL=5
 
  perturb_Ham=zeros(Float64,2*NL,2*NL)
  if perturb==1
    for layer in 1:1
      perturb_Ham[2*layer-1:2*layer,2*layer-1:2*layer]=15.0*[1 0;0 -1]
    end
    for layer in NL:NL
      perturb_Ham[2*layer-1:2*layer,2*layer-1:2*layer]=15.0*[1 0;0 -1]
    end
  elseif perturb==2
    for layer in 1:1
      perturb_Ham[2*layer-1:2*layer,2*layer-1:2*layer]=15.0*[1 0;0 -1]
    end
    for layer in NL:NL
      perturb_Ham[2*layer-1:2*layer,2*layer-1:2*layer]=-15.0*[1 0;0 -1]
    end

  elseif perturb==3
    for layer in 1:1
      perturb_Ham[2*layer-1:2*layer,2*layer-1:2*layer]=4.0*[1 0;0 -1]+28.9*[1 0;0 1]
    end
    for layer in NL:NL
      perturb_Ham[2*layer-1:2*layer,2*layer-1:2*layer]=4.0*[1 0;0 -1]+28.9*[1 0;0 1]
    end
  
  elseif perturb==4
    for layer in 1:1
      perturb_Ham[2*layer-1:2*layer,2*layer-1:2*layer]=28.9*[1 0;0 1]
    end
 
  elseif perturb==5
    for layer in 1:1
      perturb_Ham[2*layer-1:2*layer,2*layer-1:2*layer]=15.0*[1 0;0 -1]
    end
  elseif perturb==6
    for layer in 1:1
      perturb_Ham[2*layer-1:2*layer,2*layer-1:2*layer]=100.0*[1 0;0 1]
    end
  elseif perturb==7
    for layer in 1:1
      perturb_Ham[2*layer-1:2*layer,2*layer-1:2*layer]=50.0*[1 0;0 1]
    end
  elseif perturb==8
    for layer in 1:1
      perturb_Ham[2*layer-1:2*layer,2*layer-1:2*layer]=15.0*[1 0;0 -1]+40*[1 0;0 1]
    end
    for layer in NL:NL
      perturb_Ham[2*layer-1:2*layer,2*layer-1:2*layer]=15.0*[1 0;0 -1]+40*[1 0;0 1]
    end
  elseif perturb==9
    for layer in 1:1
      perturb_Ham[2*layer-1:2*layer,2*layer-1:2*layer]=15.0*[1 0;0 -1]+30*[1 0;0 1]
    end
    for layer in NL:NL
      perturb_Ham[2*layer-1:2*layer,2*layer-1:2*layer]=15.0*[1 0;0 -1]+30*[1 0;0 1]
    end
 
  
  end

  return perturb_Ham


end

function sample_value(uD::Float64, numsample::Int,E_lower::Float64,E_upper::Float64,perturb_Ham::Matrix{Float64},kradius::Float64,Einterval::Float64)
 



  ac=0.246
  R1=ac*[1,0]
  R2=ac*[1/2,√3/2]
  G1=2π/ac*[1,-1/√3]
  G2=2π/ac*[0,2/√3]




  NL=5
  KGr=4π/(3*ac)*[1,0]
 
 

  

  valuesset=[Float64[] for _ in 1:Threads.nthreads()]

  N_valence_tt=[0 for _ in 1:Threads.nthreads()]
  Threads.@threads for ja in 1:numsample
    kvec=KGr+[rand()-0.5,rand()-0.5]*2*kradius
    Hamiltonian=get_RNGham(kvec,NL,uD)+perturb_Ham
    FFF=real.(eigen(Hamiltonian).values)
    for jb in eachindex(FFF)
      if FFF[jb]>E_lower && FFF[jb]<E_upper
        push!(valuesset[Threads.threadid()],FFF[jb])
        if jb<NL+1
          N_valence_tt[Threads.threadid()]+=1
        end

      end
    end
  end
  N_valence=sum(N_valence_tt)
  
  valuesset=reduce(vcat,valuesset)
  DOS_E_binnum=Int(round((E_upper-E_lower)/Einterval))

  println("finisheddiagonalization")

  return valuesset, N_valence, DOS_E_binnum


end

function process_data(valuesset::Vector{Float64},numsample::Int,kradius::Float64,DOS_n_binnum::Int64,DOS_E_binnum::Int64,nstart::Float64,nend::Float64)
  


  h=fit(Histogram, valuesset, nbins=DOS_E_binnum) 
  bin_edges = collect(h.edges[1])
  Nstates = h.weights
  bin_centers =collect(0.5* (bin_edges[1:end-1] + bin_edges[2:end]))
  NN=numsample*(bin_centers[2]-bin_centers[1])
  Nstates=Nstates/(NN)*kradius^2*4/(4π^2)


 n_background=N_valence/numsample*kradius^2*4/(4π^2)

  nE=zeros(Float64,length(bin_centers))
 for ja in eachindex(nE)
  nE[ja]=sum(Nstates[1:ja])*(bin_centers[2]-bin_centers[1])-n_background # this is the density of carriers as a function of energy
 end




  

 nE_edges = collect(range(start=nstart, stop=nend, length=DOS_n_binnum+1))
    
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

   #nE_bin_centers,nE_bin_means is a pair, one is the carrier density, the other is the corresponding DOS 
   # bin_centers,nE is a pair, one is the energy, the other is the corresponding DOS

 return nE_bin_centers,nE_bin_means,bin_centers,nE,Nstates


end