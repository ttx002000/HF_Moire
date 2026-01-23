using LinearAlgebra,Statistics 
using LinearAlgebra,Random

function get_Ham(k::Vector{Float64},uD::Float64,valley::Int64,NL::Int)
 
   Ham=zeros(ComplexF64,2,NL,2,NL)

  t0=3100
  t1=380
  t2=-15
  t3=-290
  t4=-141
  ff=√3/2*0.246*(valley*k[1]-im*k[2])
  Delta2=2.0
  δ=10.5
  for layer in 1:NL-1
     Ham[:,layer,:,layer+1]=[t4*ff t3*conj(ff);t1 t4*ff]
  end

  if NL>2
   for layer in 1:NL-2
       Ham[:,layer,:,layer+2]=[0.0 t2/2;0.0 0.0]
   end
 end



  for layer in 1:NL
      Ham[:,layer,:,layer]+=1/2*[uD*(layer-(NL+1)/2) t0*ff;t0*conj(ff) uD*(layer-(NL+1)/2)]
  end

   for layer in 1:NL
      if layer==1 || layer==NL
       Ham[:,layer,:,layer]+=1/2*[Delta2 0; 0 Delta2]
      else
         Ham[:,layer,:,layer]-=1/2*[Delta2 0; 0 Delta2]
      end
   end

   for layer in 1:NL
      if layer==1
         Ham[:,layer,:,layer]+=1/2*[0 0; 0 δ]
      elseif layer==NL
         Ham[:,layer,:,layer]+=1/2*[δ 0; 0 0]
      else
         Ham[:,layer,:,layer]+=1/2*[δ 0; 0 δ]
      end
   end

   Ham=reshape(Ham,2*NL,2*NL)
   Ham=Ham+Ham'
  
  
 return Ham
end

struct myparams
    kx_list::Vector{Float64}
    ky_list::Vector{Float64}
    uD::Float64
    num_layers::Int64
    num_spin::Int64
    num_valley::Int64
    num_sub::Int64
    target_density::Float64
    pristine_eigenvector::Array{ComplexF64,6}
    Hamiltonian::Array{ComplexF64,6}  
    U_stoner::Float64
    ϵr::Float64
end


function objective_function(Δ_input::Vector{Float64},p::myparams)
 
        kx_list=p.kx_list###
        ky_list=p.ky_list###
        Area=4π^2/(kx_list[2]-kx_list[1])/(ky_list[2]-ky_list[1])

        num_spin=p.num_spin
        num_valley=p.num_valley
        num_sub=p.num_sub
        uD=p.uD###
        num_layers=p.num_layers###
        U_stoner=p.U_stoner###
        ϵr=p.ϵr###
        Coulomb_constant=4523.7818/ϵr
        target_density=p.target_density###
        target_index=round(Int,Area*target_density)+length(kx_list)*length(ky_list)*num_layers*num_spin*num_valley

        pristine_eigenvector=p.pristine_eigenvector###
        Hamiltonian=p.Hamiltonian###
        
        background_density=length(kx_list)*length(ky_list)*num_layers*num_spin*num_valley/Area
        norb = num_sub * num_layers
        nt = Threads.maxthreadid()
  
        n0=ones(Float64,num_spin*num_valley*num_sub*num_layers)*background_density/(num_spin*num_valley*num_sub*num_layers)
       

        Δ=copy(reshape(Δ_input,num_spin,num_valley,num_sub,num_layers))
        eigenvector=zeros(ComplexF64,length(kx_list),length(ky_list),num_spin,num_valley,num_layers*num_sub,num_layers*num_sub)
     
        eigenvalue=zeros(Float64,length(kx_list),length(ky_list),num_spin,num_valley,num_layers*num_sub)
   
        for ja in eachindex(kx_list)
            for jb in eachindex(ky_list)
                for ispin in 1:num_spin
                    for ivalley in 1:num_valley
                        Ham=Hamiltonian[ja,jb,ispin,ivalley,:,:]
                        varH=Ham+Diagonal(vec(Δ[ispin,ivalley,:,:]))
                        FFF=eigen(varH)
                        eigenvector[ja,jb,ispin,ivalley,:,:]=FFF.vectors
                        eigenvalue[ja,jb,ispin,ivalley,:]=real(FFF.values)
                    end
                end
            end
        end
    

        fermi_level=sort(vec(eigenvalue))[target_index]
        
     
        density_tls = [zeros(Float64, num_spin, num_valley, num_sub, num_layers) for _ in 1:nt]
        E1_tls = zeros(Float64, nt)
        tmp_tls = [zeros(ComplexF64, norb) for _ in 1:nt]

        Threads.@threads for ja in eachindex(kx_list)
            tid = Threads.threadid()
            dens = density_tls[tid]
            tmp  = tmp_tls[tid]

            @inbounds for jb in eachindex(ky_list)
                for ispin in 1:num_spin, ivalley in 1:num_valley
                    H     = @view Hamiltonian[ja,jb,ispin,ivalley,:,:]
                    evals = @view eigenvalue[ja,jb,ispin,ivalley,:]
                    evecs = @view eigenvector[ja,jb,ispin,ivalley,:,:]
                    evecs_pri=@view pristine_eigenvector[ja,jb,ispin,ivalley,:,:]
                    nocc  = searchsortedlast(evals, fermi_level)
                    for bandi in 1:nocc
                         evals[bandi] < fermi_level 
                            v = @view evecs[:, bandi] 
                            v2  = reshape(v, num_sub, num_layers)
                            mul!(tmp, H, v)
                           E1_tls[tid] += real(dot(v, tmp))
                           dens[ispin,ivalley,:,:] .+= abs2.(v2)
                    end
                    for bandi in 1:num_layers
                            v_pri= @view evecs_pri[:, bandi]
                            mul!(tmp,H,v_pri)
                            E1_tls[tid]-= real(dot(v_pri, tmp))
                       
                    end
                end
            end
        end
        density_resolved = reduce(+, density_tls)/Area
        E1 = sum(E1_tls)/Area

        
        
        layer_position=[ja*0.335 for ja in 1:num_layers]
        HartreeMatrix=zeros(Float64,num_spin,num_valley,num_sub,num_layers,num_spin,num_valley,num_sub,num_layers)
        for i in 1:num_layers, j in 1:num_layers
        HartreeMatrix[:,:,:,i,:,:,:,j].=-Coulomb_constant*abs(layer_position[i]-layer_position[j])
        end
        nresolved=reshape(density_resolved,num_spin*num_valley*num_sub*num_layers)
        HartreeMatrix=reshape(HartreeMatrix,num_spin*num_valley*num_sub*num_layers,num_spin*num_valley*num_sub*num_layers)
        E2=nresolved'*HartreeMatrix*nresolved-2*n0'*HartreeMatrix*nresolved+n0'*HartreeMatrix*n0
        

    
        Stoner_matrix=zeros(Float64,num_spin*num_valley,num_sub,num_layers,num_spin*num_valley,num_sub,num_layers)
        for i in 1:num_spin*num_valley, j in 1:num_spin*num_valley
            if i≠j
               Stoner_matrix[i,:,:,j,:,:].=U_stoner
            end
        end
        Stoner_matrix=reshape(Stoner_matrix,num_spin*num_valley*num_sub*num_layers,num_spin*num_valley*num_sub*num_layers)

        E3=nresolved'*Stoner_matrix*nresolved-2*n0'*Stoner_matrix*nresolved+n0'*Stoner_matrix*n0

        

  #start to calculate gradient
        ggr_tls = [zeros(Float64, num_spin, num_valley, norb) for _ in 1:nt]
        dE1dDelta=zeros(Float64, num_spin, num_valley, norb)
        # scratch per thread (avoid allocations inside loops)
        tmp_tls   = [Matrix{ComplexF64}(undef, norb, norb) for _ in 1:nt]  # HRMG*evecs
        A_tls     = [Matrix{ComplexF64}(undef, norb, norb) for _ in 1:nt]  # evecs'*HRMG*evecs
        w_tls     = [Vector{ComplexF64}(undef, norb) for _ in 1:nt]
        conjw_tls = [Vector{ComplexF64}(undef, norb) for _ in 1:nt]
        t_tls     = [Vector{ComplexF64}(undef, norb) for _ in 1:nt]


        Threads.@threads for jx in eachindex(kx_list)
            tid   = Threads.threadid()
            g_loc = ggr_tls[tid]

            tmp   = tmp_tls[tid]
            A     = A_tls[tid]
            w     = w_tls[tid]
            conjw = conjw_tls[tid]
            t     = t_tls[tid]

            @inbounds for jy in eachindex(ky_list)
                for js in 1:num_spin, jv in 1:num_valley
                    evals = @view eigenvalue[jx, jy, js, jv, :]          # (norb)
                    HRMG  = @view Hamiltonian[jx, jy, js, jv, :, :]      # (norb,norb)
                    evecs = @view eigenvector[jx, jy, js, jv, :, :]      # (norb,norb), columns = eigenvectors

                    # We update ggr at fixed (js,jv) through a view + reshape (debug-friendly)
                    gsv  = @view g_loc[js, jv, :,]                     # (num_sub,num_layers) view
                                             

                    # A = evecs' * HRMG * evecs
                    # tmp = HRMG*evecs
                    mul!(tmp, HRMG, evecs)
                    # A = evecs'*tmp
                    mul!(A, adjoint(evecs), tmp)

                    # Occupied-unoccupied sums
                    nocc  = searchsortedlast(evals, fermi_level)
                    if nocc≠0
                        for bn in 1:nocc
                    

                            εn = evals[bn]

                            # w[bm] = A[bn,bm] / (εn - εm), with w[bn]=0
                            for bm in 1:norb
                                w[bm] = (bm == bn) ? (0.0 + 0.0im) : A[bn, bm] / (εn - evals[bm])
                            end

                            # We need s[p] = Σ_m conj(evecs[p,m]) * w[m]
                            # Do: t = evecs * conj(w)  => conj(t[p]) = Σ_m conj(evecs[p,m]) * w[m]
                        
                            mul!(t, evecs, conj.(w))

                            v_n = @view evecs[:, bn]  # eigenvector n (view)

                            # g[p] += 2 * real( v_n[p] * conj(t[p]) )
                            for p in 1:norb
                                gsv[p] += 2 * real(v_n[p] * conj(t[p]))
                            end
                
                        end
                    end
                end
            end
        end

        # reduce TLS into final ggr
        
        for tid in 1:nt
            dE1dDelta .+= ggr_tls[tid]./Area
        end
        dE1dDelta=vec(dE1dDelta)





        # D_tls as before
        D_tls = [zeros(Float64, num_spin, num_valley, norb, norb) for _ in 1:nt]

        # scratch per thread
        Cocc_tls = [Matrix{Float64}(undef, norb, norb) for _ in 1:nt]         # use first nocc rows
        A_tls    = [Matrix{ComplexF64}(undef, norb, norb) for _ in 1:nt]
        Aconj_tls= [Matrix{ComplexF64}(undef, norb, norb) for _ in 1:nt]
        Zocc_tls = [Matrix{ComplexF64}(undef, norb, norb) for _ in 1:nt]      # use first nocc rows

        Threads.@threads for jx in eachindex(kx_list)
            tid   = Threads.threadid()
            D     = D_tls[tid]
            Cocc  = Cocc_tls[tid]
            A     = A_tls[tid]
            Aconj = Aconj_tls[tid]
            Zocc  = Zocc_tls[tid]

            @inbounds for jy in eachindex(ky_list)
                for js in 1:num_spin, jv in 1:num_valley
                    evals = @view eigenvalue[jx,jy,js,jv,:]
                    V     = @view eigenvector[jx,jy,js,jv,:,:]
                    nocc  = searchsortedlast(evals, fermi_level)
                if nocc≠0
                        # Cocc[n,m] for n=1:nocc only (diagonal set 0)
                        for n in 1:nocc
                            En = evals[n]
                            for m in 1:norb
                                Cocc[n,m] = (m==n) ? 0.0 : (1.0/(En - evals[m]))
                            end
                        end

                        Dsv = @view D[js,jv,:,:]

                        for α in 1:norb
                            # A[β,m] = V[β,m]*conj(V[α,m])  and Aconj = conj(A)
                            for m in 1:norb
                                cc = conj(V[α,m])
                                @simd for β in 1:norb
                                    
                                    A[β,m] = V[β,m]*cc
                                
                                end
                            end

                            mul!(@view(Zocc[1:nocc,:]), @view(Cocc[1:nocc,:]), transpose(A))

                            for β in 1:norb
                            
                            
                                ss=dot(A[β,1:nocc],Zocc[1:nocc,β])
                                
                                Dsv[β,α] += 2.0 * real(ss)
                            end
                        end
                    end
                end
            end
        end
        dnddelta=sum(D_tls)/Area

        dnddelta_matrix=zeros(num_spin,num_valley,norb,num_spin,num_valley,norb)
        for ja in 1:num_spin, jb in 1:num_valley
        dnddelta_matrix[ja,jb,:,ja,jb,:]=dnddelta[ja,jb,:,:]
        end
        dnddelta_matrix=reshape(dnddelta_matrix,num_spin*num_valley*norb,num_spin*num_valley*norb)
        dE3dDelta=dnddelta_matrix'*Stoner_matrix*nresolved+vec(nresolved'*Stoner_matrix*dnddelta_matrix)-vec(2*n0'*Stoner_matrix*dnddelta_matrix)
        dE2dDelta=dnddelta_matrix'*HartreeMatrix*nresolved+vec(nresolved'*HartreeMatrix*dnddelta_matrix)-vec(2*n0'*HartreeMatrix*dnddelta_matrix)

        
                
    return E1+E2+E3, dE1dDelta+dE2dDelta+dE3dDelta
end

function construct_parameters(radius::Float64,num_kpoints::Int64,
    uD::Float64,num_layers::Int64,target_density::Float64,Us::Float64,ϵr::Float64)

    kx_list=LinRange(-radius,radius,num_kpoints)
    ky_list=LinRange(-radius,radius,num_kpoints)
   

    num_spin=2
    num_valley=2
    num_sub=2
  



    U_stoner=Us*√3/2*0.246^2*0.5
    

    pristine_eigenvector=zeros(ComplexF64,length(kx_list),length(ky_list),num_spin,num_valley,num_layers*num_sub,num_layers*num_sub)
    Hamiltonian=zeros(ComplexF64,length(kx_list),length(ky_list),num_spin,num_valley,num_layers*num_sub,num_layers*num_sub)
    for ja in eachindex(kx_list)
        for jb in eachindex(ky_list)
            for ispin in 1:num_spin
                for ivalley in 1:num_valley
                    Ham=get_Ham([kx_list[ja],ky_list[jb]],uD,ivalley==1 ? 1 : -1,num_layers)
            
                    Hamiltonian[ja,jb,ispin,ivalley,:,:]=Ham
            

                    pristine_eigenvector[ja,jb,ispin,ivalley,:,:]=eigen(Ham).vectors
                end
            end
        end
    end

    pp=myparams(kx_list,ky_list,uD,num_layers,num_spin,num_valley,num_sub,target_density,pristine_eigenvector,Hamiltonian,U_stoner,ϵr)
   return pp
end


function GR_descent(tole::Float64,p::myparams)
        num_spin=p.num_spin
        num_valley=p.num_valley
        num_sub=p.num_sub
        num_layers=p.num_layers
        initial_Δ=10*randn(Float64,num_spin*num_valley*num_sub*num_layers)
        Gr_record=[]
        val_record=[]
        Δ_current=copy(initial_Δ)
        step_sz=0.5
        Gr_current=ones(Float64,length(Δ_current))
        itcount=0
        
        while norm(Gr_current)>tole
            itcount+=1
            println("Iteration number: ",itcount)
            val_current,Gr_current=objective_function(Δ_current,p)
            push!(Gr_record,copy(Gr_current))
            println("Current gradient norm: ",norm(Gr_current))
            push!(val_record,val_current)
            println("Current objective value: ",val_current)

            Δ_current.=Δ_current.-step_sz*Gr_current

            if itcount>20
                if val_record[end]>val_record[end-1]-0.8*step_sz*dot(Gr_record[end],Gr_record[end])
                    step_sz=step_sz*0.8
                    println("Reducing step size to ",step_sz)
                else
                    step_sz=step_sz*1.2
                    println("Increasing step size to ",step_sz)
                end
            end
            flush(stdout)
        end


        return val_record,Gr_record,Δ_current

end