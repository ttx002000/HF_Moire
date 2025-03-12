
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


function find_E_cut(nstart::Float64,nend::Float64,uD::Float64,perturb_Ham::Matrix{Float64})
ac=0.246
KGr=4π/(3*ac)*[1,0]
NL=5


kx_grid=collect(-1.5:0.01:1.5)
ky_grid=collect(-1.5:0.01:1.5)
conduction_edge_set=[]
valence_edge_set=[]
for ja in [1,length(kx_grid)],jb in eachindex(ky_grid)
    kvec=[kx_grid[ja],ky_grid[jb]]
    Hamiltonian=get_RNGham(kvec+KGr,NL,uD)+perturb_Ham
    FFF=eigen(Hamiltonian)
    push!(conduction_edge_set,FFF.values[NL+1])
    push!(valence_edge_set,FFF.values[NL])
end

conduction_edge=sort(conduction_edge_set)[1]
valence_edge=sort(valence_edge_set)[end]

println(conduction_edge)
println(valence_edge)






Area=4π^2/((kx_grid[2]-kx_grid[1])*(ky_grid[2]-ky_grid[1]))
knorm=[]
Eval=[]
N_BG=0
for ja in eachindex(kx_grid),jb in eachindex(ky_grid)
    kvec=[kx_grid[ja],ky_grid[jb]]
    Hamiltonian=get_RNGham(KGr+kvec,NL,uD)
    FFF=eigen(Hamiltonian)
    for jc in eachindex(FFF.values)
        if valence_edge<FFF.values[jc]<conduction_edge
          push!(knorm,norm(kvec))
          push!(Eval,FFF.values[jc])
          if jc<NL+1
             N_BG+=1
          end
        end
    end


end
ord=sortperm(Eval)
Eval=Eval[ord]
knorm=knorm[ord]

nrange=collect(1:1:length(Eval))/Area.-N_BG/Area

if nend>nrange[end] || nstart<nrange[1]
  error("boundary wrong, try larger sampling area")
end

E_lower=Eval[searchsortedfirst(nrange,nstart)]
E_upper=Eval[searchsortedfirst(nrange,nend)]

k1=[]
for ja in eachindex(Eval)
   if E_upper-15<Eval[ja]<E_upper+15
     push!(k1,knorm[ja])
   end

end

for ja in eachindex(Eval)
  if E_lower-15<Eval[ja]<E_lower+15
    push!(k1,knorm[ja])
  end
end
kradius=sort(k1)[end]


return E_lower-5,E_upper+5, kradius*1.1
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

  N_valence=0
  Threads.@threads for ja in 1:numsample
    kvec=KGr+[rand()-0.5,rand()-0.5]*2*kradius
    Hamiltonian=get_RNGham(kvec,NL,uD)+perturb_Ham
    FFF=real.(eigen(Hamiltonian).values)
    for jb in eachindex(FFF)
      if FFF[jb]>E_lower && FFF[jb]<E_upper
        push!(valuesset[Threads.threadid()],FFF[jb])
        if jb<NL+1
          N_valence+=1
        end

      end
    end
  end
  
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