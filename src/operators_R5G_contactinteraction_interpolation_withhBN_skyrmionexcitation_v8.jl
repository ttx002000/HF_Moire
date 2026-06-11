using LinearAlgebra
using Arpack
using Combinatorics
using Random
BLAS.set_num_threads(1)


function get_Ham(k::Vector{Float64},uD::Float64,valley::Int64,stacking::Int,NL::Int)
 
   Ham=zeros(ComplexF64,2*NL,2*NL)
  Kac=4π/(3*0.246)*[1,0]*valley
  t0=3100
  t1=380
  t2=-21
  t3=290
  t4=141
  for layer in 1:NL-1
     Ham[2*layer-1:2*layer,2*layer+1:2*layer+2]=[t4*get_f((k+Kac)*stacking) t3*conj(get_f((k+Kac)*stacking));t1 t4*get_f((k+Kac)*stacking)]
  end

  if NL>2
   for layer in 1:NL-2
      Ham[2*layer-1:2*layer,2*layer+3:2*layer+4]=[0.0 t2/2;0.0 0.0]
   end
 end

  Ham=Ham+Ham'

  for layer in 1:NL
      Ham[2*layer-1:2*layer,2*layer-1:2*layer]=[uD*(layer-(NL+1)/2) -t0*get_f((k+Kac)*stacking);-t0*conj(get_f((k+Kac)*stacking)) uD*(layer-(NL+1)/2)]
  end

  
  
 return Ham
end




function get_Ham_Holomorphic(k::Vector{Float64},uD::Float64,valley::Int64,stacking::Int,NL::Int)
 
  Ham=zeros(ComplexF64,2*NL,2*NL)

  t0=3100
  t1=380
  fk=-√3/2*0.246*k[1]*valley+stacking*im*√3/2*0.246*k[2]

  for layer in 1:NL-1
     Ham[2*layer-1:2*layer,2*layer+1:2*layer+2]+=[0.0 0.0;t1 0.0]
  end



  Ham=Ham+Ham'

  for layer in 1:NL
      Ham[2*layer-1:2*layer,2*layer-1:2*layer]+=[0.0 -t0*fk;-t0*conj(fk) 0.0]
  end
  
   Ham[1:2,1:2]=[0.0 0.0; 0.0 0.0]

 for layer in 1:NL
   Ham[2*layer,2*layer]+=eigen(get_Ham(k,uD,valley,stacking,NL)).values[NL+1] 
 end

  for layer in 1:NL
   Ham[2*layer-1,2*layer-1]-=eigen(get_Ham(k,uD,valley,stacking,NL)).values[NL+1]  
 end
  
 return Ham
end




 

function get_f(k::Vector{Float64})
  delta1=1/√3*0.246*[0,1]
  delta2=1/√3*0.246*[√3/2,-1/2]
  delta3=1/√3*0.246*[-√3/2,-1/2]

  return exp(im*dot(k,delta1))+exp(im*dot(k,delta2))+exp(im*dot(k,delta3))
end





function sendtomesh(Minv::Matrix{Int64},q1::Vector{Int})::Vector{Int}
  Qvec=q1'*inv(Minv)
  return q1 .-vec((Int.(floor.(round.(Qvec,digits=5)))*Minv)')
end




function Geometry(geonum::Int64)
  if geonum==1
      Nx=3;
      Ny=3;
      l1=[3,0]
      l2=[0,3]
      Tr1=[1,0]
      Tr2=[0,1]
  end

  if geonum==2
      Nx=3;
      Ny=3;
      l1=[1,1]*3
      l2=[-1,2]*3
      Tr1=[1,1]
      Tr2=[-1,2]
  end

   if geonum==3
      Nx=6;
      Ny=6;
      l1=[6,0]
      l2=[0,6]
      Tr1=[1,0]
      Tr2=[0,1]
  end


  if geonum==4
      Nx=6;
      Ny=6;
      l1=[1,1]*6
      l2=[-1,2]*6
      Tr1=[1,1]
      Tr2=[-1,2]
  end
  
  

  return Nx,Ny,l1,l2,Tr1,Tr2
end



function triangle_initial_Densitymatrix(NL::Int,θ::Float64,geonum::Int64,gcutoff::Float64,
                                      uD::Float64,λ::Float64,
                                        V0_hBN::Float64,V1_hBN::Float64,ψ_hBN::Float64, 
                                        V2_scalar::Float64,ϕ::Float64,
                                        V3_scalar::Float64,ϕ3::Float64,r3::Int)
   
    aGr=0.246
    ϵ=0.2504/aGr-1 #This is the normal one
  
    G1=2π/aGr*[1,-1/√3]
    G2=2π/aGr*[0,2/√3]
    Nx,Ny,l1,l2,Tr1,Tr2=Geometry(geonum)

    Rθ=[cos(θ) -sin(θ);sin(θ) cos(θ)]
    b1_ps=(G1-(1+ϵ)^(-1)*Rθ*G1)
    b2_ps=(G2-(1+ϵ)^(-1)*Rθ*G2)
    a1m_ps=inv([b1_ps';b2_ps'])*[2π,0] #This is the pristine a
    a2m_ps=inv([b1_ps';b2_ps'])*[0,2π]
    
    a1m=Tr1[1]*a1m_ps+Tr1[2]*a2m_ps;  #This is the unit cell a
    a2m=Tr2[1]*a1m_ps+Tr2[2]*a2m_ps;
    b1=inv([a1m';a2m'])*[2π,0]
    b2=inv([a1m';a2m'])*[0,2π]
      
    L1=l1[1]*a1m_ps+l1[2]*a2m_ps; #This is the torus
    L2=l2[1]*a1m_ps+l2[2]*a2m_ps;

    Area=abs(L1[1]*L2[2]-L2[1]*L1[2]);
    
    Rotminus90=[0 1;-1 0]
    T1=2*π/Area*Rotminus90*L2
    T2=-2*π/Area*Rotminus90*L1


    
    
    am=norm(a1m);
    

    b1T=Int.(round.(inv([T1 T2])*b1))
    b2T=Int.(round.(inv([T1 T2])*b2))

    b1_ps_T=Int.(round.(inv([T1 T2])*b1_ps))
    b2_ps_T=Int.(round.(inv([T1 T2])*b2_ps))


    
    Rotation120=[cos(2π/3) -sin(2π/3);sin(2π/3) cos(2π/3)]

    R120_b1T=Int.(round.(inv([T1 T2])*Rotation120*b1))
   

   
 
    
 
    
    

    
    

    
    
    
   allowedq=Vector{Int64}[]
  for ja in 0:Nx-1,jb in 0:Ny-1
      push!(allowedq,[ja,jb])
  end
  
  for ja in 1:Nx*Ny
    allowedq[ja]=sendtomesh([b1T';b2T'],allowedq[ja])
  end
    
  allowedq=unique(allowedq)
  if length(allowedq)<Nx*Ny
   throw("there is an error in allowedq")
  end
    
    
    wave=Vector{Int64}[]
    
    cutoffstandard=gcutoff*norm(b1)
    cutoff=Int(round(gcutoff))*5
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        gtest=ja*b1+jb*b2;
        if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
            push!(wave,ja*b1T+jb*b2T)
        end
    end
    dimension=length(wave)
   
    
    
    spinor_set=Matrix{Vector{ComplexF64}}(undef,length(allowedq),length(wave))
    for ja in eachindex(allowedq), jb in eachindex(wave)
      hh=(1-λ)*get_Ham_Holomorphic([T1 T2]*(allowedq[ja]+wave[jb]),uD,1,1,NL)+(λ)*get_Ham([T1 T2]*(allowedq[ja]+wave[jb]),uD,1,1,NL)
      v1=eigen(hh).vectors[:,NL+1]
      if uD>0.0
        v1angle=angle(v1[2*NL])
      else
        throw("there is an error")
      end
      spinor_set[ja,jb]=v1*exp(-im*v1angle)/norm(v1)
    end



     overlapmatrix=zeros(ComplexF64,length(allowedq),length(wave),length(allowedq),length(wave))
    for ja in eachindex(allowedq)
      for jb in eachindex(wave), jc in eachindex(allowedq), jd in eachindex(wave)
       overlapmatrix[ja,jb,jc,jd]=spinor_set[ja,jb]'*spinor_set[jc,jd]
      end
    end



    
     
    single_Ham=[zeros(ComplexF64,dimension,dimension) for _ in eachindex(allowedq)]
    single_MoirePo=[zeros(ComplexF64,dimension,dimension) for _ in eachindex(allowedq)]
    single_eigenvalue=[zeros(Float64,dimension) for _ in eachindex(allowedq)]
    single_eigenvector=[zeros(ComplexF64,dimension,dimension) for _ in eachindex(allowedq)]

         ω=exp(im*2π/3)
       op_1=zeros(ComplexF64,2*NL,2*NL)
       op_1[1:2,1:2]=[1 1;ω ω]

       op_2=zeros(ComplexF64,2*NL,2*NL)
       op_2[1:2,1:2]=[1 ω^2;ω^2 ω]

       op_3=zeros(ComplexF64,2*NL,2*NL)
       op_3[1:2,1:2]=[1 ω;1 ω]

       
       op_4=zeros(ComplexF64,2*NL,2*NL)
       op_4[1:2,1:2]=[1 0; 0 1]

       op_5=Matrix{Float64}(I,2*NL,2*NL)


       if r3==1
        rcenter=a1m/3
       else
        error("not implemented r3")
       end
    
    Threads.@threads for ja in eachindex(allowedq)
      
      
      k=allowedq[ja][1]*T1+allowedq[ja][2]*T2
    
      for jb in eachindex(wave)
          hh=(1-λ)*get_Ham_Holomorphic([T1 T2]*(allowedq[ja]+wave[jb]),uD,1,1,NL)+(λ)*get_Ham([T1 T2]*(allowedq[ja]+wave[jb]),uD,1,1,NL)
   
         
          single_Ham[ja][jb,jb]=real(eigen(hh).values[NL+1])

         
      end


     






     # adding ordinary potential
     for jc in eachindex(wave)
    
        pos=findfirst(item->item==wave[jc]-b1_ps_T,wave)
        if pos≠nothing
       
          single_MoirePo[ja][jc,pos]+=V1_hBN*exp(-im*ψ_hBN)*(spinor_set[ja,jc]'*op_1*spinor_set[ja,pos])
          single_MoirePo[ja][jc,pos]+=V2_scalar*exp(-im*ϕ)*(spinor_set[ja,jc]'*op_5*spinor_set[ja,pos])
        end
        
        pos=findfirst(item->item==wave[jc]-b2_ps_T,wave)
        if pos≠nothing

            single_MoirePo[ja][jc,pos]+=V1_hBN*exp(-im*ψ_hBN)*(spinor_set[ja,jc]'*op_2*spinor_set[ja,pos])
            single_MoirePo[ja][jc,pos]+=V2_scalar*exp(-im*ϕ)*(spinor_set[ja,jc]'*op_5*spinor_set[ja,pos])
        end


        pos=findfirst(item->item==wave[jc]+(b1_ps_T+b2_ps_T),wave)
        if pos≠nothing

            single_MoirePo[ja][jc,pos]+=V1_hBN*exp(-im*ψ_hBN)*(spinor_set[ja,jc]'*op_3*spinor_set[ja,pos])
            single_MoirePo[ja][jc,pos]+=V2_scalar*exp(-im*ϕ)*(spinor_set[ja,jc]'*op_5*spinor_set[ja,pos])
        end
    
        single_MoirePo[ja][jc,jc]+=V0_hBN/2*(spinor_set[ja,jc]'*op_4*spinor_set[ja,jc])
     end

     # adding symmetry breaking potential

     for jc in eachindex(wave)
    
        pos=findfirst(item->item==wave[jc]-b1T,wave)
        if pos≠nothing
       
      
          single_MoirePo[ja][jc,pos]+=V3_scalar*exp(-im*ϕ3)*exp(-im*dot([T1 T2]*b1T,rcenter))*(spinor_set[ja,jc]'*op_5*spinor_set[ja,pos])
        end
        
        pos=findfirst(item->item==wave[jc]-R120_b1T,wave)
        if pos≠nothing
            single_MoirePo[ja][jc,pos]+=V3_scalar*exp(-im*ϕ3)*exp(-im*dot([T1 T2]*R120_b1T,rcenter))*(spinor_set[ja,jc]'*op_5*spinor_set[ja,pos])
        end


        pos=findfirst(item->item==wave[jc]+(b1T+R120_b1T),wave)
        if pos≠nothing
            single_MoirePo[ja][jc,pos]+=V3_scalar*exp(-im*ϕ3)*exp(-im*dot([T1 T2]*(-b1T-R120_b1T),rcenter))*(spinor_set[ja,jc]'*op_5*spinor_set[ja,pos])
        end
    
     end





    
     single_MoirePo[ja]=single_MoirePo[ja]+single_MoirePo[ja]'
     FFF=eigen(single_MoirePo[ja]+single_Ham[ja])
    
      single_eigenvalue[ja]=real(FFF.values)
      single_eigenvector[ja]=FFF.vectors
    
    end


    input_DensityMatrix=[zeros(ComplexF64,dimension,dimension) for _ in eachindex(allowedq)]
 
    
     for ja in eachindex(allowedq)
      A=randn(ComplexF64,length(wave),length(wave))
      input_DensityMatrix[ja]=(A+A')*1.0
     end

    

    
    return overlapmatrix, wave, input_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,single_eigenvector,allowedq, T1, T2, a1m, a2m, b1,b2,spinor_set,Area,b1T,b2T,a1m_ps,a2m_ps,b1_ps,b2_ps,b1_ps_T,b2_ps_T
      

       
end



mutable struct ConstructDM_Workspace
    HartreeMatrix::Vector{Matrix{ComplexF64}}
    FockMatrix::Vector{Matrix{ComplexF64}}
    H_phys::Vector{Matrix{ComplexF64}}
    output_DensityMatrix::Vector{Matrix{ComplexF64}}
    DeltaMatrix::Vector{Matrix{ComplexF64}}
    NewDensityMatrix::Vector{Matrix{ComplexF64}}
    HF_eigenvalue::Vector{Vector{Float64}}
    HF_eigenvector::Vector{Matrix{ComplexF64}}

    Hartree_Density::Matrix{ComplexF64}   # dimension×dimension
    all_eigs::Vector{Float64}             # length = length(allowedq)*dimension


    DIIS_input_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}}
    DIIS_input_DeltaMatrix::Vector{Vector{Matrix{ComplexF64}}}
    diis_head::Int
    diis_len::Int
end


function init_ConstructDM_Workspace(allowedq, wave,DIIS_size)
    dimension = length(wave)
    Nk = length(allowedq)

    HartreeMatrix = [zeros(ComplexF64, dimension, dimension) for _ in 1:Nk]
    FockMatrix    = [zeros(ComplexF64, dimension, dimension) for _ in 1:Nk]
   H_phys    = [zeros(ComplexF64, dimension, dimension) for _ in 1:Nk]
    NewDensityMatrix = [zeros(ComplexF64, dimension, dimension) for _ in 1:Nk]

    output_DensityMatrix = [zeros(ComplexF64, dimension, dimension) for _ in 1:Nk]
    DeltaMatrix          = [zeros(ComplexF64, dimension, dimension) for _ in 1:Nk]

    HF_eigenvalue  = [Vector{Float64}(undef, dimension) for _ in 1:Nk]
    HF_eigenvector = [Matrix{ComplexF64}(undef, dimension, dimension) for _ in 1:Nk]

    Hartree_Density = zeros(ComplexF64, dimension, dimension)
    all_eigs = Vector{Float64}(undef, Nk * dimension)
    DIIS_input_DensityMatrix=[[zeros(ComplexF64, dimension, dimension) for _ in 1:Nk] for _ in 1:DIIS_size]
     DIIS_input_DeltaMatrix=[[zeros(ComplexF64, dimension, dimension) for _ in 1:Nk] for _ in 1:DIIS_size]


    return ConstructDM_Workspace(
        HartreeMatrix, FockMatrix,H_phys, output_DensityMatrix, DeltaMatrix,
        NewDensityMatrix, HF_eigenvalue, HF_eigenvector, Hartree_Density, all_eigs,
        DIIS_input_DensityMatrix,  DIIS_input_DeltaMatrix,1, 0)
end

struct LoopDic
    # Fock quadruples (flat)
    g1F::Vector{Int}
    g2F::Vector{Int}
    g3F::Vector{Int}
    g4F::Vector{Int}

    # Hartree quadruples (flat)
    g1H::Vector{Int}
    g2H::Vector{Int}
    g3H::Vector{Int}
    g4H::Vector{Int}

    # Hartree Coulomb per Hartree-entry
    hartree_val::Vector{Float64}

    # For fast Fock Coulomb lookup in Construct_DensityMatrix:
    dqidF::Vector{Int}           # length = length(g1F)
    dkid::Matrix{Int}            # Nk×Nk
    coulomb_dkdq::Matrix{Float64}# Ndk×Ndq
end





function construct_loop_dic(wave, allowedq, T1, T2, ϵr, constq,gateD)::LoopDic
    Ng = length(wave)
    Nk = length(allowedq)
    nthreads = Threads.maxthreadid()

    # -------------------------
    # Fock quadruples (threaded build)
    # -------------------------
    g1F_buf = [Int[] for _ in 1:nthreads]
    g2F_buf = [Int[] for _ in 1:nthreads]
    g3F_buf = [Int[] for _ in 1:nthreads]
    g4F_buf = [Int[] for _ in 1:nthreads]

    Threads.@threads for ja in 1:Ng
        tid = Threads.threadid()
        b1 = g1F_buf[tid]; b2 = g2F_buf[tid]; b3 = g3F_buf[tid]; b4 = g4F_buf[tid]

        for jb in 1:Ng, jc in 1:Ng, jd in 1:ja
            if wave[ja] + wave[jb] == wave[jc] + wave[jd]
                push!(b1, ja); push!(b2, jb); push!(b3, jc); push!(b4, jd)
            end
        end
    end

    nF = 0
    @inbounds for t in 1:nthreads
        nF += length(g1F_buf[t])
    end
    g1F = Vector{Int}(undef, nF)
    g2F = Vector{Int}(undef, nF)
    g3F = Vector{Int}(undef, nF)
    g4F = Vector{Int}(undef, nF)

    off = 0
    @inbounds for t in 1:nthreads
        len = length(g1F_buf[t])
        if len > 0
            copyto!(g1F, off+1, g1F_buf[t], 1, len)
            copyto!(g2F, off+1, g2F_buf[t], 1, len)
            copyto!(g3F, off+1, g3F_buf[t], 1, len)
            copyto!(g4F, off+1, g4F_buf[t], 1, len)
            off += len
        end
    end
    g1F_buf=nothing
    g2F_buf=nothing
    g3F_buf=nothing
    g4F_buf=nothing

    # -------------------------
    # Hartree quadruples (threaded build)
    # -------------------------
    g1H_buf = [Int[] for _ in 1:nthreads]
    g2H_buf = [Int[] for _ in 1:nthreads]
    g3H_buf = [Int[] for _ in 1:nthreads]
    g4H_buf = [Int[] for _ in 1:nthreads]

    Threads.@threads for ja in 1:Ng
        tid = Threads.threadid()
        b1 = g1H_buf[tid]; b2 = g2H_buf[tid]; b3 = g3H_buf[tid]; b4 = g4H_buf[tid]

        for jb in 1:Ng, jc in 1:ja, jd in 1:Ng
            if wave[ja] + wave[jb] == wave[jc] + wave[jd]
                push!(b1, ja); push!(b2, jb); push!(b3, jc); push!(b4, jd)
            end
        end
    end

    nH = 0
    @inbounds for t in 1:nthreads
        nH += length(g1H_buf[t])
    end
    g1H = Vector{Int}(undef, nH)
    g2H = Vector{Int}(undef, nH)
    g3H = Vector{Int}(undef, nH)
    g4H = Vector{Int}(undef, nH)

    off = 0
    @inbounds for t in 1:nthreads
        len = length(g1H_buf[t])
        if len > 0
            copyto!(g1H, off+1, g1H_buf[t], 1, len)
            copyto!(g2H, off+1, g2H_buf[t], 1, len)
            copyto!(g3H, off+1, g3H_buf[t], 1, len)
            copyto!(g4H, off+1, g4H_buf[t], 1, len)
            off += len
        end
    end

    g1H_buf=nothing
    g2H_buf=nothing
    g3H_buf=nothing
    g4H_buf=nothing

    # -------------------------
    # Hartree Coulomb per Hartree-entry (threaded)
    # -------------------------
    hartree_val = Vector{Float64}(undef, nH)
    Threads.@threads for w in 1:nH
        q = wave[g3H[w]] - wave[g1H[w]]
        hartree_val[w] = Coulomb(q, T1, T2, gateD)/ϵr + constq
    end

    # -------------------------
    # dq indexing for Fock (setup; keep simple, single-thread)
    # -------------------------
    #pack(q1,q2) = (UInt64(reinterpret(UInt32, Int32(q1))) << 32) | UInt64(reinterpret(UInt32, Int32(q2)))

    dq_map = Dict{Tuple{Int,Int},Int}()
    dq1_list = Int[]; dq2_list = Int[]
    dqidF = Vector{Int}(undef, nF)

    @inbounds for w in 1:nF
        dq = wave[g3F[w]] - wave[g1F[w]]
        key = (dq[1], dq[2])
        id = get(dq_map, key, 0)
        if id == 0
            push!(dq1_list, dq[1]); push!(dq2_list, dq[2])
            id = length(dq1_list)
            dq_map[key] = id
        end
        dqidF[w] = id
    end
    Ndq = length(dq1_list)

    # -------------------------
    # dk indexing for (jk,jk1) (setup; single-thread)
    # -------------------------
    dk_map = Dict{Tuple{Int,Int},Int}()
    dk1_list = Int[]; dk2_list = Int[]
    dkid = Matrix{Int}(undef, Nk, Nk)

    @inbounds for jk in 1:Nk, jk1 in 1:Nk
        dk1 = allowedq[jk1][1] - allowedq[jk][1]
        dk2 = allowedq[jk1][2] - allowedq[jk][2]
         key = (dk1, dk2)
        id = get(dk_map, key, 0)
        if id == 0
            push!(dk1_list, dk1); push!(dk2_list, dk2)
            id = length(dk1_list)
            dk_map[key] = id
        end
        dkid[jk,jk1] = id
    end
    Ndk = length(dk1_list)

    # -------------------------
    # Coulomb table over (dk,dq) (threaded)
    # -------------------------
    coulomb_dkdq = Matrix{Float64}(undef, Ndk, Ndq)
    Threads.@threads for idk in 1:Ndk
        dk1 = dk1_list[idk]
        dk2 = dk2_list[idk]
        @inbounds for idq in 1:Ndq
            q1 = dk1 + dq1_list[idq]
            q2 = dk2 + dq2_list[idq]
            coulomb_dkdq[idk,idq] = Coulomb([q1,q2], T1, T2,gateD)/ϵr + constq
        end
    end

    return LoopDic(g1F,g2F,g3F,g4F, g1H,g2H,g3H,g4H, hartree_val, dqidF, dkid, coulomb_dkdq)
end












function Coulomb(k::Vector{Int},T1::Vector{Float64},T2::Vector{Float64},gateD::Float64)::Float64
   D=gateD
   return k==[0,0] ? D*9047.5636 : tanh(norm([T1 T2]*k*D))/norm(k[1]*T1+k[2]*T2)*9047.5636
end


function Construct_DensityMatrix(work::ConstructDM_Workspace,loop::LoopDic,
                              allowedq::Vector{Vector{Int}},
                              wave::Vector{Vector{Int64}},
                               input_DensityMatrix::Vector{Matrix{ComplexF64}},single_Ham::Vector{Matrix{ComplexF64}},
                               single_MoirePo::Vector{Matrix{ComplexF64}},overlapmatrix::Array{ComplexF64,4},
                               energy_input::Float64,filling::Int,Area::Float64)
  
  
    dimension=length(wave)
    Nk = length(allowedq)
 
   HartreeMatrix = work.HartreeMatrix
    FockMatrix    = work.FockMatrix
    output_DensityMatrix = work.output_DensityMatrix
    DeltaMatrix          = work.DeltaMatrix
    NewDensityMatrix     = work.NewDensityMatrix
    HF_eigenvalue        = work.HF_eigenvalue
    HF_eigenvector       = work.HF_eigenvector
    Hartree_Density      = work.Hartree_Density
    all_eigs             = work.all_eigs
    H_phys               = work.H_phys

      @inbounds for ja in 1:Nk
          fill!(HartreeMatrix[ja], 0)
          fill!(FockMatrix[ja],    0)
          fill!(NewDensityMatrix[ja], 0)
        
      end
      fill!(Hartree_Density, 0)
      


   nF = length(loop.g1F)
    Threads.@threads for jk in 1:Nk
        Fk = FockMatrix[jk]
        for jk1 in 1:Nk
            dmk = input_DensityMatrix[jk1]
            opf = @view overlapmatrix[jk1, :, jk,  :]
            opm = @view overlapmatrix[jk,  :, jk1, :]

            idk = loop.dkid[jk, jk1]

            @inbounds for w in 1:nF
                g1 = loop.g1F[w]
                g2 = loop.g2F[w]
                g3 = loop.g3F[w]
                g4 = loop.g4F[w]

                cc = loop.coulomb_dkdq[idk, loop.dqidF[w]]

                Fk[g1, g4] += dmk[g3, g2] * opm[g1, g3] * opf[g2, g4] * cc
            end
        end
    end







   #Hartree_Density=sum([input_DensityMatrix[ja] .* transpose(overlapmatrix[ja,:,ja,:]) for ja in eachindex(allowedq)])
   @inbounds for ja in 1:Nk
        dmk = input_DensityMatrix[ja]
        oph = @view overlapmatrix[ja, :, ja, :]
        for j in 1:dimension, i in 1:dimension
            # transpose(oph)[j,i] == oph[i,j]
            Hartree_Density[j,i] += dmk[j,i] * oph[i,j]
        end
    end

 nH = length(loop.g1H)
    Threads.@threads for jk in 1:Nk
        Hk  = HartreeMatrix[jk]
        oph = @view overlapmatrix[jk, :, jk, :]

        @inbounds for w in 1:nH
            g1 = loop.g1H[w]
            g2 = loop.g2H[w]
            g3 = loop.g3H[w]
            g4 = loop.g4H[w]

            Hk[g1, g3] += Hartree_Density[g4, g2] * oph[g1, g3] * loop.hartree_val[w]
        end
    end



  Threads.@threads for jk in eachindex(allowedq)
    symmetrize_from_lower!(HartreeMatrix[jk])
    symmetrize_from_lower!(FockMatrix[jk])
    HartreeMatrix[jk] ./= Area
    FockMatrix[jk]    ./= Area
  end
  Threads.@threads for jk in eachindex(allowedq)
    copy!(H_phys[jk] , single_MoirePo[jk])
        H_phys[jk]  .+= single_Ham[jk]
       H_phys[jk]  .+= work.HartreeMatrix[jk]
       H_phys[jk] .-= work.FockMatrix[jk]
  end


   idx = 0
  @inbounds for ja in eachindex(allowedq)
        FFF = eigen(Hermitian(work.H_phys[ja]))
        vals = real(FFF.values)
        copy!(HF_eigenvalue[ja],vals)
        copy!(HF_eigenvector[ja], FFF.vectors)

        for j in 1:dimension
            idx += 1
            all_eigs[idx] = vals[j]
        end
  end

    nocc = filling * length(allowedq)
    partialsort!(all_eigs, 1:nocc+1)
    bound = 0.5 * (all_eigs[nocc] + all_eigs[nocc+1])

     mix_ratio=0.5

  for ja in eachindex(allowedq)
    occ=searchsortedlast(HF_eigenvalue[ja],bound)
    @views Vocc= HF_eigenvector[ja][:,1:occ]
     mul!(work.NewDensityMatrix[ja], Vocc, adjoint(Vocc))
     
       symmetrize_from_lower!(NewDensityMatrix[ja])
       @. DeltaMatrix[ja]=NewDensityMatrix[ja]-input_DensityMatrix[ja]
     
       @. output_DensityMatrix[ja]=mix_ratio*input_DensityMatrix[ja]+(1-mix_ratio)*NewDensityMatrix[ja]
  end


  
  e1=0.0
  for ja in eachindex(allowedq)
    e1+=sum(abs2,DeltaMatrix[ja])
  end
  eout=e1/length(allowedq)
  
  
  energy=0
  for ja in eachindex(allowedq)
  energy += real(dot(input_DensityMatrix[ja], single_MoirePo[ja])) +
          real(dot(input_DensityMatrix[ja], single_Ham[ja])) +
          0.5*real(dot(input_DensityMatrix[ja], HartreeMatrix[ja])) -
          0.5*real(dot(input_DensityMatrix[ja], FockMatrix[ja]))
  end

   energy_change=real(energy-energy_input)



 return  eout, energy_change, energy
end

@inline function symmetrize_from_lower!(A::Matrix{ComplexF64})
    n = size(A,1)
    @inbounds for i in 1:n
        A[i,i] = complex(real(A[i,i]), 0.0)
        for j in (i+1):n
            A[i,j] = conj(A[j,i])   # fill UPPER from LOWER
        end
    end
    return A
end

function iteration_loop(initial_DensityMatrix::Vector{Matrix{ComplexF64}},
                       allowedq::Vector{Vector{Int64}},T1::Vector{Float64},T2::Vector{Float64},
                       wave::Vector{Vector{Int64}},single_Ham::Vector{Matrix{ComplexF64}},
                       single_MoirePo::Vector{Matrix{ComplexF64}},constq::Float64,ϵr::Float64,overlapmatrix::Array{ComplexF64,4},
                       filling::Int,Area::Float64,gateD::Float64)
    eout=1.0
    itcount=0
    dimension=length(wave)
    DIIS_size=3

     #loop_dic_Fock,loop_dic_Fock_val,loop_dic_Hartree,loop_dic_Hartree_val=construct_loop_dic(wave,allowedq,T1,T2,ϵr,constq)
    loop = construct_loop_dic(wave, allowedq, T1, T2, ϵr, constq,gateD)
     work=init_ConstructDM_Workspace(allowedq, wave,DIIS_size)
    
   
    
      input_DensityMatrix=initial_DensityMatrix
    bad_count=0
    energy=0.0
    energy_change=0.0

    eout_hist = Float64[]
     PLATEAU_N = 10
     PLATEAU_FRAC = 0.10
     E_EPS = 1e-30

    diis_fire_once = false
    diis_cooldown = 0              # prevent immediate re-trigger after DIIS
     DIIS_COOLDOWN = 10 


    while (eout>1*10^(-18)) || (bad_count<DIIS_size) || (abs(energy_change)>1*10^(-9))

      if eout<1*10^(-18)
       bad_count+=1
      else
        bad_count=0
      end
       
         if (itcount > 500 && abs(eout) > 10)
                itcount = 0

                # forget DIIS/plateau state completely
                bad_count = 0
                energy = 0.0
                energy_change = 0.0
                empty!(eout_hist)
                diis_fire_once = false
                diis_cooldown = 0
                work.diis_head = 1
                work.diis_len  = 0

                # random restart DM (your style)
                for ja in eachindex(allowedq)
                  A = randn(dimension, dimension) + im*randn(dimension, dimension)
                  input_DensityMatrix[ja] = (A + A') * 0.01
                end
                dmk_used=input_DensityMatrix
                println("random start again")
                
        end

      
      tic=time()
      dmk_used = input_DensityMatrix
      
      if (itcount>100 && abs(eout)>10^(-2)) || (itcount>30 && abs(eout)<10^(-10)) || diis_fire_once
      
        dmk=implement_DIIS(work.DIIS_input_DensityMatrix,work.DIIS_input_DeltaMatrix,allowedq,DIIS_size)
        dmk_used=dmk

        eout, energy_change, energy=Construct_DensityMatrix(work, loop,allowedq,wave,dmk,single_Ham,single_MoirePo,overlapmatrix,energy,filling,Area)
       
       
        println("using DIIS")
        

         if diis_fire_once
            diis_fire_once = false
            empty!(eout_hist)          # <-- yes: clear history after firing
            diis_cooldown = DIIS_COOLDOWN
          end
      else
  
      

        eout, energy_change, energy=Construct_DensityMatrix(work,loop,allowedq,wave,input_DensityMatrix,single_Ham,single_MoirePo,overlapmatrix,energy,filling,Area)
        
        
       
         
     

      end

    
      work.diis_head, work.diis_len = diis_push!(
        work.DIIS_input_DensityMatrix,
             work.DIIS_input_DeltaMatrix,
              dmk_used,              # <--- this is the fix: store the DM that was used
              work.DeltaMatrix,
              work.diis_head,
              work.diis_len
                )

          
      copy!(input_DensityMatrix, work.output_DensityMatrix)


      itcount+=1
     
      toc=time()
      println(toc-tic,"eout=$eout","energy_change=$energy_change","itcount=$itcount")
      flush(stdout)

      if diis_cooldown > 0
            diis_cooldown -= 1
      end


       push!(eout_hist, eout)  # keep sign; we'll use abs where needed
        if length(eout_hist) > PLATEAU_N
            popfirst!(eout_hist)
        end

        # plateau detection (only if not cooling down and not already scheduled)
        if !diis_fire_once && diis_cooldown == 0 && length(eout_hist) == PLATEAU_N
            e0 = eout_hist[1]
            e1 = eout_hist[end]
            rel_change = abs(e1 - e0) / max(abs(e0), E_EPS)

            if rel_change < PLATEAU_FRAC
                diis_fire_once = true
                println("Plateau detected: |Δe|/|e| ≈ $(rel_change). Will fire DIIS once.")
            end
        end
     
    
  end







  
    return work.DIIS_input_DensityMatrix,
       work.DIIS_input_DeltaMatrix,
       work.HF_eigenvalue,
       work.HF_eigenvector,
       energy, eout,
       work.HartreeMatrix,
       work.FockMatrix
  
end



@inline function diis_push!(
    DIIS_input_DensityMatrix::Vector{Vector{Matrix{ComplexF64}}},
    DIIS_input_DeltaMatrix::Vector{Vector{Matrix{ComplexF64}}},
    input_DensityMatrix::Vector{Matrix{ComplexF64}},
    DeltaMatrix::Vector{Matrix{ComplexF64}},
    diis_head::Int,
    diis_len::Int
)
    copy!(DIIS_input_DensityMatrix[diis_head], input_DensityMatrix)
    copy!(DIIS_input_DeltaMatrix[diis_head],  DeltaMatrix)

    diis_head = (diis_head == length(DIIS_input_DensityMatrix)) ? 1 : (diis_head + 1)
    diis_len  = min(diis_len + 1, length(DIIS_input_DensityMatrix))
    return diis_head, diis_len
end








function shift_vector(vector_toshift::Array{ComplexF64,3},shiftamount::Vector{Int},wave::Vector{Vector{Int}},NL::Int,filling::Int)
    
  
  shifted_vector=zeros(ComplexF64,2*NL,length(wave),filling)
    for ja in eachindex(wave)
        pos=findfirst(item->item==wave[ja]+shiftamount,wave)
        if pos≠nothing
        shifted_vector[:,ja,:]=vector_toshift[:,pos,:]
        end
    
    end

  shifted_vector=reshape(shifted_vector,2*NL*length(wave),filling)
  for ja in 1:filling
   shifted_vector[:,ja]=shifted_vector[:,ja]/norm(shifted_vector[:,ja])
  end

    return shifted_vector

end



function triangle_chern(wave::Vector{Vector{Int}},
       allowedq::Vector{Vector{Int}},
       HF_eigenvector::Vector{Matrix{ComplexF64}},
       single_eigenvector::Vector{Matrix{ComplexF64}},spinor_set::Matrix{Vector{ComplexF64}},
       geonum::Int,b1T::Vector{Int},b2T::Vector{Int},NL::Int,filling::Int)

   Nx,Ny,_,_,_,_=Geometry(geonum)
    

    
    
    chern_allowedq=Vector{Int64}[]
    for ja in 0:Nx,jb in 0:Ny
        push!(chern_allowedq,[ja,jb])
    end

    PW_HF_eigenvector=[zeros(ComplexF64,2*NL,length(wave),filling) for _ in eachindex(allowedq)]
    PW_single_eigenvector=[zeros(ComplexF64,2*NL,length(wave),filling) for _ in eachindex(allowedq)]

    for ja in eachindex(allowedq), jb in eachindex(wave), jc in 1:filling
        PW_HF_eigenvector[ja][:,jb,jc]=HF_eigenvector[ja][jb,jc]*spinor_set[ja,jb]
        PW_single_eigenvector[ja][:,jb,jc]=single_eigenvector[ja][jb,jc]*spinor_set[ja,jb]
    end



  

    eigenvector_intermediate_bc=Vector{Array{ComplexF64,2}}(undef,(Nx+1)*(Ny+1))
    eigenvector_intermediate_single=Vector{Array{ComplexF64,2}}(undef,(Nx+1)*(Ny+1))

    Threads.@threads for ja in eachindex(chern_allowedq)
      
      
      kbraket=sendtomesh([b1T';b2T'],chern_allowedq[ja])
      kbraket_pos=findfirst(item->item==kbraket,allowedq)
      if kbraket==chern_allowedq[ja]
        eigenvector_intermediate_bc[ja]=reshape(PW_HF_eigenvector[kbraket_pos],2*NL*length(wave),filling)
        eigenvector_intermediate_single[ja]=reshape(PW_single_eigenvector[kbraket_pos],2*NL*length(wave),filling)
      else
        eigenvector_intermediate_bc[ja]=shift_vector(PW_HF_eigenvector[kbraket_pos],chern_allowedq[ja]-kbraket,wave,NL,filling)
        eigenvector_intermediate_single[ja]=shift_vector(PW_single_eigenvector[kbraket_pos],chern_allowedq[ja]-kbraket,wave,NL,filling)

      end  
    end

    eigenvector_bc=zeros(ComplexF64,length(wave)*2*NL,filling,Nx+1,Ny+1)
    eigenvector_bc_single=zeros(ComplexF64,length(wave)*2*NL,filling,Nx+1,Ny+1)
    for ja in eachindex(chern_allowedq)
       eigenvector_bc[:,:,chern_allowedq[ja][1]+1,chern_allowedq[ja][2]+1]=eigenvector_intermediate_bc[ja]
       eigenvector_bc_single[:,:,chern_allowedq[ja][1]+1,chern_allowedq[ja][2]+1]=eigenvector_intermediate_single[ja]
    end
    
    
    Uonelink=zeros(ComplexF64,Nx,Ny+1)
    Utwolink=zeros(ComplexF64,Nx+1,Ny)
    
    
    for ja in 1:Nx, jb in 1:Ny+1
       Uonelink[ja,jb]=det(eigenvector_bc[:,:,ja,jb]'*eigenvector_bc[:,:,ja+1,jb])/abs(det(eigenvector_bc[:,:,ja,jb]'*eigenvector_bc[:,:,ja+1,jb]))
    end

    
    
    for ja in 1:Nx+1, jb in 1:Ny
      Utwolink[ja,jb]=det(eigenvector_bc[:,:,ja,jb]'*eigenvector_bc[:,:,ja,jb+1])/abs(det(eigenvector_bc[:,:,ja,jb]'*eigenvector_bc[:,:,ja,jb+1]))
    end
    
    
   


    Flink=zeros(ComplexF64,Nx,Ny)
    for ja in 1:Nx, jb in 1:Ny
     Flink[ja,jb]=log(Uonelink[ja,jb]*Utwolink[ja+1,jb]/(Uonelink[ja,jb+1]*Utwolink[ja,jb]))
    end
    chern=sum(Flink)/(2*π*im)
    aveF=sum(Flink)/(Nx*Ny)
    uniform=0
    for ja in 1:Nx, jb in 1:Ny
        uniform+=(imag(Flink[ja,jb])-imag(aveF))^2*(Nx*Ny)/(2π)^2
    end


    Uonelink=zeros(ComplexF64,Nx,Ny+1)
    Utwolink=zeros(ComplexF64,Nx+1,Ny)
  
    
  for ja in 1:Nx, jb in 1:Ny+1
       Uonelink[ja,jb]=det(eigenvector_bc_single[:,:,ja,jb]'*eigenvector_bc_single[:,:,ja+1,jb])/abs(det(eigenvector_bc_single[:,:,ja,jb]'*eigenvector_bc_single[:,:,ja+1,jb]))
    end

    
    
    for ja in 1:Nx+1, jb in 1:Ny
      Utwolink[ja,jb]=det(eigenvector_bc_single[:,:,ja,jb]'*eigenvector_bc_single[:,:,ja,jb+1])/abs(det(eigenvector_bc_single[:,:,ja,jb]'*eigenvector_bc_single[:,:,ja,jb+1]))
    end
    

   

    Flink_single=zeros(ComplexF64,Nx,Ny)
    for ja in 1:Nx, jb in 1:Ny
     Flink_single[ja,jb]=log(Uonelink[ja,jb]*Utwolink[ja+1,jb]/(Uonelink[ja,jb+1]*Utwolink[ja,jb]))
    end
    chern_single=sum(Flink_single)/(2*π*im)
    aveF=sum(Flink_single)/(Nx*Ny)
    uniform_single=0
    for ja in 1:Nx, jb in 1:Ny
        uniform_single+=(imag(Flink_single[ja,jb])-imag(aveF))^2*(Nx*Ny)/(2π)^2
    end
    
    

    
 
      


    
      
   return chern,Flink,chern_single,Flink_single,uniform,uniform_single

end



function implement_DIIS(DIIS_input_projector::Vector{Vector{Matrix{ComplexF64}}},DIIS_input_DeltaMatrix::Vector{Vector{Matrix{ComplexF64}}},allowedq::Vector{Vector{Int}},DIIS_size::Int)



      Bmatrix=zeros(ComplexF64,DIIS_size+1,DIIS_size+1)
      for ja in 1:DIIS_size
       Bmatrix[ja,DIIS_size+1]=1
       Bmatrix[DIIS_size+1,ja]=1
      end
  
      for ja in 1:DIIS_size,jb in 1:DIIS_size
          for jc in eachindex(allowedq)
              Bmatrix[ja,jb]+=real(dot(DIIS_input_DeltaMatrix[ja][jc],DIIS_input_DeltaMatrix[jb][jc]))
          end
      end

      inB=safe_inverse(Bmatrix)
      if inB≠0
         onh=zeros(Float64,DIIS_size+1)
         onh[end]=1
         coeff=inB* onh
         
        AA =coeff[1]*(DIIS_input_projector[1]+DIIS_input_DeltaMatrix[1])           # alloc once
        @inbounds for ja in 2:DIIS_size
          AA+=coeff[ja]*(DIIS_input_projector[ja]+DIIS_input_DeltaMatrix[ja])  # dmk += coeff[ja] * projector[ja]
        end
         return AA
      else
        return 0
      end
end


function safe_inverse(A)
  try
      return inv(A)  # Attempt to compute inverse
  catch e
      if isa(e, SingularException)
          return pinv(A, 10^(-8))  # Use pseudoinverse as an alternative
      else
          rethrow(e)  # If another error occurs, propagate it
      end
  end
end