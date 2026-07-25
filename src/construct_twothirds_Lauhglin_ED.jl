using EllipticFunctions
using LinearAlgebra
using Plots
using JLD2,StaticArrays
using Combinatorics
using SparseArrays


function Cdag(op_index::Int64,state_index::Int64,sign::Int64)::Tuple{Int64,Int64}
    if (sign==0) | ((2^(op_index-1)& state_index)==2^(op_index-1))
    return 0,0
    end
 
    new_sign=sign*(-1)^count_ones((2^(op_index-1)-1) & (state_index))
    new_state_index=state_index+2^(op_index-1)
   
    return new_sign,new_state_index
end





function Cann(op_index::Int64,state_index::Int64,sign::Int64)::Tuple{Int64,Int64}
    if (sign==0) | ((2^(op_index-1)& state_index)≠2^(op_index-1))
    return 0,0
    end
 
    new_sign=sign*(-1)^count_ones((2^(op_index-1)-1) & (state_index))
    new_state_index=state_index-2^(op_index-1)
   
    return new_sign,new_state_index
end



function sendtomesh(Minv::Matrix{Int64},q1::Vector{Int})::Vector{Int}
  Qvec=q1'*inv(Minv)
  return q1 .-vec((Int.(floor.(round.(Qvec,digits=5)))*Minv)')
end

function v_cross(v1::Vector{Float64},v2::Vector{Float64})
   return v1[1]*v2[2]-v1[2]*v2[1]
end



function BuildVmatrix(Nx::Int64,Ny::Int64,wave::Vector{Vector{Int64}},allowedq::Vector{Vector{Int64}},T1::Vector{Float64},T2::Vector{Float64},Deltatheta::Vector{Float64},g1T,g2T,dimension,Lb,Minv,meshpos)
  Fmatrix=Array{ComplexF64}(undef,Nx*Ny,Nx*Ny,dimension)
  invmatrix=inv([g1T g2T])
  
  for  jk1 in 1:Nx*Ny, jq in 1:Nx*Ny, jqg in 1:dimension
    
    k1vec=[T1 T2]*allowedq[jk1]
    k2meshT=sendtomesh(Minv,allowedq[jk1]+allowedq[jq])
    k2vec=[T1 T2]*allowedq[meshpos[k2meshT]]
    qvec=[T1 T2]*(allowedq[jq]+wave[jqg])
    gT=(allowedq[jk1]+allowedq[jq]+wave[jqg]-allowedq[meshpos[k2meshT]])
    gn=invmatrix*gT
    @assert norm(gn-round.(gn))<10^(-9)
    gvec=[T1 T2]*gT

  
    sgnf=(-1)^(round(gn[1])*round(gn[2]))

    Fmatrix[jk1,jq,jqg]=exp(-Lb^2/4*norm(qvec)^2)*sgnf*exp(im*Lb^2/2*v_cross(k1vec-Deltatheta,gvec))*exp(im*Lb^2/2*v_cross(k2vec-Deltatheta,qvec))
   
  end
  
  
  Vmatrix=zeros(ComplexF64,Nx*Ny,Nx*Ny,Nx*Ny,Nx*Ny)


 
  
  
  for jk1 in 1:Nx*Ny,jk2 in 1:(jk1-1),qmesh in 1:Nx*Ny
     
   
        k3mesh_pos=meshpos[sendtomesh(Minv,allowedq[jk1]+allowedq[qmesh])]
        k4mesh_pos=meshpos[sendtomesh(Minv,allowedq[jk2]-allowedq[qmesh])]
        if k3mesh_pos≠nothing && k4mesh_pos≠nothing 
          for qg in 1:length(wave)
              mqT=-1*allowedq[qmesh]-wave[qg]
              
              mqmesh=sendtomesh(Minv,mqT)

              mqgT=mqT-mqmesh
              mqmesh_pos=findfirst(item->item==mqmesh,allowedq)
              mqg_pos=findfirst(item->item==mqgT,wave)
              qvec=(allowedq[qmesh][1]+wave[qg][1])*T1+(allowedq[qmesh][2]+wave[qg][2])*T2
              
              if mqg_pos≠nothing && mqmesh_pos≠nothing && norm(qvec)≠0.0
                Vmatrix[jk1,jk2,k3mesh_pos,k4mesh_pos]+=-(norm(qvec))^2*Fmatrix[jk1,qmesh,qg]*Fmatrix[jk2,mqmesh_pos,mqg_pos]
              end
            
          end
       end
     
   end

   return Vmatrix

end


function reducedV(Vmatrix::Array{ComplexF64},Nx::Int64,Ny::Int64)
  reduced_Vcol=ComplexF64[]
  reduced_Vcoor=Vector{Int64}[]
  
  for i in 1:Nx*Ny, j in 1:i-1, k in 1:Nx*Ny, p in 1:k-1
  if abs(Vmatrix[i,j,k,p]-Vmatrix[i,j,p,k])>10^(-10)
  push!(reduced_Vcol,2*Vmatrix[i,j,k,p]-2*Vmatrix[i,j,p,k])
  push!(reduced_Vcoor,[i,j,k,p])
  end
  end
  
  return reduced_Vcol, reduced_Vcoor
end

function do_ED(Nx::Int64,Ny::Int64,Nparticle::Int,
              flux1::Float64,flux2::Float64,q1::Float64,q2::Float64,
              num_vecs::Int)

        a1m=[1,0]
        a2m=[1/2,√3/2]


        g1=(2*π*[1 0; 0 1]*inv([a1m a2m]))[1,:]
        g2=(2*π*[1 0; 0 1]*inv([a1m a2m]))[2,:]


       
        L1=Nx*a1m
        L2=Ny*a2m;
     
        Lb=((a1m[1]*a2m[2]-a2m[1]*a1m[2])/(2π))^(1/2)


     
        area=abs(L1[1]*L2[2]-L2[1]*L1[2]);
        Rotminus90=[0 1;-1 0]

        T1=2*π/area*Rotminus90*L2
        T2=-2*π/area*Rotminus90*L1


        Deltatheta=flux1*T1+flux2*T2
        qq=q1*g1+q2*g2


        Deltatheta_LL=Deltatheta+qq





        g1T=Int.(round.(inv([T1 T2])*g1))
        g2T=Int.(round.(inv([T1 T2])*g2))

        Minv=[g1T';g2T']



        wave=Vector{Int64}[]
        cutoff=18
        cutoffstandard=8.01*norm(g1)
        for ja in -cutoff:cutoff, jb in -cutoff:cutoff
                gtest=ja*g1+jb*g2;
                if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
                    push!(wave,ja*g1T+jb*g2T)
                end
        end

        dimension=length(wave)
        allowedq=Vector{Int64}[]
            
        allowedq=Vector{Int64}[]
        for jb in 1:Ny, ja in 1:Nx
            push!(allowedq,[ja-1,jb-1])
        end


        for ja in 1:Ny*Nx
            allowedq[ja]=sendtomesh([g1T';g2T'],allowedq[ja])
        end

        meshpos=Dict{Vector{Int},Int}()
        for ja in eachindex(allowedq)
            meshpos[allowedq[ja]]=ja
        end


        Vmatrix=BuildVmatrix(Nx,Ny,wave,allowedq,T1,T2,Deltatheta_LL,g1T,g2T,dimension,Lb,Minv,meshpos)
        (reduced_Vcol,reduced_Vcoor)=reducedV(Vmatrix,Nx,Ny)


        values_record=Vector{Vector{ComplexF64}}(undef,Nx*Ny)
        eig_vec_record=Vector{Matrix{ComplexF64}}(undef,Nx*Ny)
        state_can_record=Vector{Vector{Vector{Int}}}(undef,Nx*Ny)
        state_integer_record=Vector{Vector{Int}}(undef,Nx*Ny)
        for sector in 1:Nx*Ny
            (state_can_record[sector], state_integer_record[sector])=Construct_MBstate(Nx,Ny,Nparticle,allowedq,Minv,meshpos,sector);
            values_record[sector],eig_vec_record[sector]=Construct_Manybodymatrix(reduced_Vcol,reduced_Vcoor,state_can_record[sector],state_integer_record[sector])
        end



        v1=[]
        v2=[]
        v3=[]
        for ja in 1:Nx*Ny, jb in eachindex(values_record[ja])
            push!(v1,ja)
            push!(v2,values_record[ja][jb])
            push!(v3,[ja,jb])

        end
        ord=sortperm(real.(v2))
        v1=v1[ord]
        v2=v2[ord]
        v3=v3[ord]


        
        Laughlin_state_vector=Vector{Vector{ComplexF64}}(undef,num_vecs)
        Laughlin_state_sector=Vector{Int}(undef,num_vecs)
        Laughlin_state_can=Vector{Vector{Vector{Int}}}(undef,num_vecs)
        Laughlin_state_integer=Vector{Vector{Int}}(undef,num_vecs)

        for ja in 1:num_vecs
            sector_in=v3[ja][1]
            vc_in=v3[ja][2]
            Laughlin_state_vector[ja]=eig_vec_record[sector_in][:,vc_in]
            Laughlin_state_sector[ja]=sector_in
            Laughlin_state_can[ja]=state_can_record[sector_in]
            Laughlin_state_integer[ja]=state_integer_record[sector_in]
        end

        PH_eigvector=deepcopy(Laughlin_state_vector)
        PH_state_can=Vector{Vector{Vector{Int}}}(undef,num_vecs) # I didn't sort this in order. But that's fine. We don't care about it anyway

        for ja in 1:num_vecs
            PH_state_can[ja]=[take_complement(Laughlin_state_can[ja][jb],Nx*Ny) for jb in eachindex(Laughlin_state_can[ja])]
                for jb in eachindex(Laughlin_state_vector[ja])
                    sgn=(-1)^(sum(Laughlin_state_can[ja][jb])-Nparticle)
                    PH_eigvector[ja][jb]=conj(Laughlin_state_vector[ja][jb])*sgn   
                end
         
        end

  

        return PH_eigvector,PH_state_can,allowedq,T1,T2, 
               values_record,eig_vec_record, state_can_record, state_integer_record,
                Laughlin_state_vector,Laughlin_state_sector,Laughlin_state_can, Laughlin_state_integer,
                a1m,a2m,g1,g2

     




end


function take_complement(vv::Vector{Int},Nk::Int)
    return  sort(setdiff(collect(1:1:Nk),vv))
end


function Construct_MBstate(Nx::Int64,Ny::Int64,Nparticle::Int64,allowedq::Vector{Vector{Int64}},Minv::Matrix{Int},meshpos::Dict{Vector{Int},Int},sector::Int64)



  MB_state=collect(combinations(1:Nx*Ny,Nparticle))
  MB_state_can=Vector{Int64}[]
  MB_state_integer=Int64[]
  target=allowedq[sector]
   
   for ja in eachindex(MB_state)
       QN=sum(allowedq[MB_state[ja]])
       QN_M=sendtomesh(Minv,QN)

       if QN_M==target
        push!(MB_state_can,MB_state[ja])
        push!(MB_state_integer,sum(2 .^ (MB_state[ja].-1)))
       end
     end
   

   
   
     sortindex=sortperm(MB_state_integer)
     MB_state_integer=MB_state_integer[sortindex]
     MB_state_can=MB_state_can[sortindex]


  return MB_state_can, MB_state_integer



end



function Construct_Manybodymatrix(reduced_Vcol::Vector{ComplexF64},reduced_Vcoor::Vector{Vector{Int64}},state_can,state_integer)

    

  
   
  Rows = [Vector{Int64}() for _ in 1:Threads.nthreads()]
  Cols = [Vector{Int64}() for _ in 1:Threads.nthreads()]
  Vals = [Vector{ComplexF64}() for _ in 1:Threads.nthreads()]

  Threads.@threads for jc in eachindex(reduced_Vcol)
    for jb in eachindex(state_can)
   
      (sign,state)=Cann(reduced_Vcoor[jc][3],state_integer[jb],1)
      (sign,state)=Cann(reduced_Vcoor[jc][4],state,sign)
      (sign,state)=Cdag(reduced_Vcoor[jc][2],state,sign)
      (sign,state)=Cdag(reduced_Vcoor[jc][1],state,sign)
 
       state_index=searchsortedfirst(state_integer,state)
       
      if (state*sign)≠0 &&  (state_integer[state_index]==state) 
       push!(Rows[Threads.threadid()],state_index)
       push!(Cols[Threads.threadid()],jb)
       push!(Vals[Threads.threadid()],sign*reduced_Vcol[jc]/2)
      end
   end

  end




  matrix_index1=reduce(vcat,Rows)
  matrix_index2=reduce(vcat,Cols)
  matrix_value=reduce(vcat,Vals)
  
  Rows=nothing;
  Cols=nothing;
  Vals=nothing;
  GC.gc()

 

    Hsp = sparse(
        matrix_index1,
        matrix_index2,
        matrix_value,
        length(state_can),
        length(state_can),
    )

    H = Hermitian(Matrix(Hsp))

    F = eigen(H)

    MB_spectrum = F.values
    ζ = F.vectors





return MB_spectrum,ζ

end