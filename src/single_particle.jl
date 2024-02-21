 
 using LinearAlgebra


function energy_state(flux::Float64)
    
    g1=4*π/(√3)*[1,0]
    g2=4*π/(√3)*[1/2,√3/2]
    g3=4*π/(√3)*[-1/2,√3/2]
    g4=-g1
    g5=-g2
    g6=-g3
    #a1m=[√3/2,1/2]
    #a2m=[0,1]
    a1m=[√3/2,1/2]
    a2m=[-√3/2,1/2]
    Lb=(√3/(4π))^(1/2)
    
    
    
    
    l1=[5,0]
    l2=[0,6]
    Nx=5;
    Ny=6;
    Nparticle=10;
    L1=l1[1]*a1m+l1[2]*a2m;
    L2=l2[1]*a1m+l2[2]*a2m;
    area=abs(L1[1]*L2[2]-L2[1]*L1[2]);
    Rotminus90=[0 1;-1 0]
    T1=2*π/area*Rotminus90*L2
    T2=-2*π/area*Rotminus90*L1
    g1T=Int.(round.(inv([T1 T2])*g1))
    g3T=Int.(round.(inv([T1 T2])*g3))
    

    wave=Vector{Int64}[]
    cutoff=18
    cutoffstandard=8.01*4*π/(√3)
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        gtest=ja*g1+jb*g3;
        if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
            push!(wave,ja*g1T+jb*g3T)
        end
    end
    dimension=length(wave)
    LL_num=1


    allowedq=Vector{Int64}[]
    for jb in 1:Ny, ja in 1:Nx
    push!(allowedq,[ja-1,jb-1])
    end
    
    
    eigenvalue_single=Vector{Float64}(undef,Nx*Ny)
    eigenvector_single=Matrix{ComplexF64}(undef,LL_num,Nx*Ny)
    for ja in 1:Nx*Ny
       
        eigenvalue_single[ja]=0.5
        eigenvector_single[:,ja]=[1.0] 
    end
    
    return eigenvalue_single, eigenvector_single, wave, T1, T2, g1T, g3T,allowedq,Nx,Ny,Nparticle,dimension,Lb
    
     
    



end