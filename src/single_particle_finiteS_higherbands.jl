using Combinatorics


using Combinatorics
using LinearAlgebra
using Plots



function overlapFS(k::Vector{Float64},q::Vector{Float64},spin::Float64,M::Float64)::ComplexF64
   
    v=(M^2+norm(k)^2+k[1]*q[1]+k[2]*q[2]-im*(k[1]*q[2]-k[2]*q[1]))^(Int(2*spin))/((M^2+norm(k)^2)^spin*(M^2+norm(k+q)^2)^spin)
    
    return v
 end
  
 
 
 function metricFS(wavelist::Vector{Vector{Int64}},spin::Float64,vf::Float64,k::Vector{Float64},q::Vector{Float64},T1::Vector{Float64},T2::Vector{Float64})::Matrix{ComplexF64}
    Amatrix=zeros(ComplexF64,length(wavelist),length(wavelist))
    mass=0.5
    for ja in 1:length(wavelist)
    Amatrix[ja,ja]=overlapFS(k+wavelist[ja][1]*T1+wavelist[ja][2]*T2,q,spin,mass*vf/2)
    end
    return Amatrix
 end
 






function getwave(cutoffstandard,g1,g3,g1T,g3T)
    wave=Vector{Int64}[]
    cutoff=18

    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        gtest=ja*g1+jb*g3;
        if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
            push!(wave,ja*g1T+jb*g3T)
        end
    end

    return wave
end

function getMoireglist(wave,g1T,g3T)
    Moireglist=[Int64[] for _ in 1:length(wave)]
    for ja in eachindex(wave)
        pos=findfirst(item->item==wave[ja]-g3T,wave)
        if pos≠nothing
        push!(Moireglist[ja],pos)
        end
        if isnothing(pos)
            push!(Moireglist[ja],0)
        end
    
        pos=findfirst(item->item==wave[ja]-g1T,wave)
        if pos≠nothing
        push!(Moireglist[ja],pos)
        end
        if isnothing(pos)
            push!(Moireglist[ja],0)
        end
    
        pos=findfirst(item->item==wave[ja]+g3T+g1T,wave)
        if pos≠nothing
        push!(Moireglist[ja],pos)
        end
        if isnothing(pos)
            push!(Moireglist[ja],0)
        end
    end
    return Moireglist
end


function getHamiltonian(k::Vector{Float64},wave::Vector{Vector{Int64}},Moireglist::Vector{Vector{Int64}},T1::Vector{Float64},T2::Vector{Float64},g1::Vector{Float64},g3::Vector{Float64},spin::Float64,vf::Float64,V0::Float64,ϕ::Float64)
    mass=0.5
    dimension=length(wave)
    Ham=zeros(ComplexF64,dimension,dimension)
    MoirePo=zeros(ComplexF64,dimension,dimension)
    
    for jb in 1:length(wave)
    Ham[jb,jb]=norm(k+wave[jb][1]*T1+wave[jb][2]*T2)^2/(2*mass)
    end
   
    for jc in 1:length(wave)
        k1=k+wave[jc][1]*T1+wave[jc][2]*T2
        
       if Moireglist[jc][1]≠0 
        MoirePo[jc,Moireglist[jc][1]]=V0*exp(im*ϕ)*overlapFS(k1,-g3,spin,mass*vf/2)
       end
       if Moireglist[jc][2]≠0 
        MoirePo[jc,Moireglist[jc][2]]=V0*exp(im*ϕ)*overlapFS(k1,-g1,spin,mass*vf/2)
       end
        if Moireglist[jc][3]≠0 
        MoirePo[jc,Moireglist[jc][3]]=V0*exp(im*ϕ)*overlapFS(k1,g1+g3,spin,mass*vf/2)
        end
    
      end
  

    MoirePo=MoirePo+MoirePo'

    return Ham+MoirePo
end




function calculatechern(wave::Vector{Vector{Int64}},Moireglist::Vector{Vector{Int64}},T1::Vector{Float64},T2::Vector{Float64},g1::Vector{Float64},g3::Vector{Float64},spin::Float64,vf::Float64,V0::Float64,ϕ::Float64,bandindex::Int64)::Tuple{ComplexF64,Float64,Matrix{ComplexF64},Float64,Float64}

   
    Nchern=15
    dimension=length(wave)
    chern_allowedq=Vector{Int64}[]
    for jb in 1:Nchern+1, ja in 1:Nchern+1
    push!(chern_allowedq,[ja,jb])
    end
    
    
    single_energy=Array{Float64}(undef,2,(Nchern+1)^2)
    chern_eigenvector_single=Array{ComplexF64}(undef,dimension,Nchern+1,Nchern+1)
    for ja in 1:(Nchern+1)*(Nchern+1)
    
        Ham=getHamiltonian(chern_allowedq[ja][1]/Nchern*g1+chern_allowedq[ja][2]/Nchern*g3,wave,Moireglist,T1,T2,g1,g3,spin,vf,V0,ϕ)
        F1=eigen(Ham)
        chern_eigenvector_single[:,chern_allowedq[ja][1],chern_allowedq[ja][2]]=F1.vectors[:,bandindex]
        single_energy[1,ja]=real(F1.values[bandindex])
        single_energy[2,ja]=real(F1.values[bandindex+1])
    end
    
   bandwidth=max(single_energy[1,:]...)-min(single_energy[1,:]...)
   gap=min((single_energy[2,:]-single_energy[1,:])...)

    
    Uonelink=zeros(ComplexF64,Nchern,Nchern+1)
    Utwolink=zeros(ComplexF64,Nchern+1,Nchern)
     
       
       for ja in 1:Nchern, jb in 1:Nchern+1
        Amatrix=metricFS(wave,spin,vf,ja/Nchern*g1+jb/Nchern*g3,1/Nchern*g1,T1,T2)
          Uonelink[ja,jb]=dot(chern_eigenvector_single[:,ja,jb],Amatrix*chern_eigenvector_single[:,ja+1,jb])/abs(dot(chern_eigenvector_single[:,ja,jb],Amatrix*chern_eigenvector_single[:,ja+1,jb]))
       end
    
       
       
       for ja in 1:Nchern+1, jb in 1:Nchern
        Amatrix=metricFS(wave,spin,vf,ja/Nchern*g1+jb/Nchern*g3,1/Nchern*g3,T1,T2)
        Utwolink[ja,jb]=dot(chern_eigenvector_single[:,ja,jb],Amatrix*chern_eigenvector_single[:,ja,jb+1])/abs(dot(chern_eigenvector_single[:,ja,jb],Amatrix*chern_eigenvector_single[:,ja,jb+1]))
       end
       
     
    
      
    
    
       Flink=zeros(ComplexF64,Nchern,Nchern)
       for ja in 1:Nchern, jb in 1:Nchern
        Flink[ja,jb]=log(Uonelink[ja,jb]*Utwolink[ja+1,jb]/(Uonelink[ja,jb+1]*Utwolink[ja,jb]))
       end
       
       aveF=sum(Flink)/Nchern^2
       uniform=0
       for ja in 1:Nchern, jb in 1:Nchern
           uniform+=(imag(Flink[ja,jb])-imag(aveF))^2*Nchern^2/(2π)^2
       end
    
       chern=sum(Flink)/(2*π*im)

  return chern,uniform,Flink,bandwidth,gap
end



function calculateformq(k::Vector{Float64},Trqpath::Vector{Vector{Float64}},wave::Vector{Vector{Int64}},Moireglist::Vector{Vector{Int64}},T1::Vector{Float64},T2::Vector{Float64},g1::Vector{Float64},g3::Vector{Float64},spin::Float64,vf::Float64,V0::Float64,ϕ::Float64,bandindex::Int64)::Vector{Float64}
  
    
    dimension=length(wave)
    Nqpath=length(Trqpath)
    eigenvector_single=Vector{ComplexF64}(undef,dimension)
    Ham=getHamiltonian(k,wave,Moireglist,T1,T2,g1,g3,spin,vf,V0,ϕ)
    
    F1=eigen(Ham)
    eigenvector_single=F1.vectors[:,bandindex]
        
    
    eigenvector_single_q=Matrix{ComplexF64}(undef,dimension,Nqpath)
    for jqpath=1:Nqpath
    
       
        Hamq=getHamiltonian(k+Trqpath[jqpath],wave,Moireglist,T1,T2,g1,g3,spin,vf,V0,ϕ)
        F1=eigen(Hamq)
        eigenvector_single_q[:,jqpath]=F1.vectors[:,bandindex]
        
    end
    
    
    formq=zeros(Float64,Nqpath)
 
   for ja in 1:Nqpath
      Amatrix=metricFS(wave,spin,vf,k,Trqpath[ja],T1,T2)
      formq[ja]+=abs(dot(eigenvector_single[:],Amatrix*eigenvector_single_q[:,ja]))^2
   end
    
    

  return formq
end