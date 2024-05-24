function get_Ham(k::Vector{Float64},δ::Float64)
    Ham=[norm(k)^2+δ/2 k[1]-im*k[2];k[1]+im*k[2] -(norm(k)^2+δ/2)]
    return Ham
end

function get_MoireHam(k::Vector{Float64},δ::Float64,wave::Vector{Vector{Int64}},T1::Vector{Float64},T2::Vector{Float64},V0::Float64,ϕ::Float64)
    dimension=2*length(wave)
    diagHam=zeros(ComplexF64,dimension,dimension)
    for ja in eachindex(wave)
       kvec=k+[T1 T2]*wave[ja]
       diagHam[2*(ja-1)+1:2*ja,2*(ja-1)+1:2*ja]=get_Ham(kvec,δ)
    end

    MoirePo=zeros(ComplexF64,dimension,dimension)
    for ja in eachindex(wave)
    
       pos=findfirst(item->item==wave[ja]-g2T,wave)
       if pos≠nothing
           MoirePo[2*(ja-1)+1:2*(ja-1)+2,2*(pos-1)+1:2*(pos-1)+2]=[V0*exp(im*ϕ) 0;0 V0*exp(im*ϕ)]
       end
       
     
       pos=findfirst(item->item==wave[ja]-g1T,wave)
       if pos≠nothing
           MoirePo[2*(ja-1)+1:2*(ja-1)+2,2*(pos-1)+1:2*(pos-1)+2]=[V0*exp(im*ϕ) 0;0 V0*exp(im*ϕ)]
       end
   
       pos=findfirst(item->item==wave[ja]+g1T+g2T,wave)
       if pos≠nothing
           MoirePo[2*(ja-1)+1:2*(ja-1)+2,2*(pos-1)+1:2*(pos-1)+2]=[V0*exp(im*ϕ) 0;0 V0*exp(im*ϕ)]
       end
    end

    return diagHam+MoirePo+MoirePo'
 
end




function kpprojected_Ham(L1::Float64)
    
    #=
    #CdAS

    #L1=80;
    
    M1=0.0205;
    A1=0;
    A2=0.889;
    B1=18.77;
    B2=13.5;
    C1=-0.0145;
    D1=10.59;
    D2=11.5;
    =#



  #SbTe

    #L1=25;
    M1=-0.22;
    A1=0.84;
    A2=3.4;
    B1=-19.64;
    B2=-48.51;
    C1=0.001;
    D1=-12.39;
    D2=-10.78;

    
    sx=[0 1;1 0]
    sy=[0 -im;im 0]
    sz=[1 0;0 -1]
    Ncut=8;
     

  tsz=diagm([1.0,-1.0,-1.0,1.0])
  OP=zeros(Float64,4*Ncut,4*Ncut);
    for jc in 1:Ncut, jd in 1:4
    
         OP[4*(jc-1)+jd,4*(jc-1)+jd]=(-1)^(jc)*(tsz[jd,jd]);

    end
     FFF=eigen(OP);
     Pmatrix=FFF.vectors[:,1:2*Ncut];


     #Let's construct the basis states at k=0




     Ham0=zeros(ComplexF64,4*Ncut,4*Ncut);
     kx=0;
     ky=0;
     k=0;
 
 
 
     for jb in 1:Ncut, jc in 1:Ncut
             if jc≠jb
                 fac=jb*jc/(jc^2-jb^2)*((-1)^(jc+jb)-1)*2/L1;
                 Ham0[4*(jb-1)+1:4*(jb-1)+4,4*(jc-1)+1:4*(jc-1)+4]+=
                 [-im*A1*fac*sx zeros(ComplexF64,2,2);zeros(ComplexF64,2,2) +im*A1*fac*sx]; 
         end
     end
 
 
         for jc in 1:Ncut
             fac=-jc^2*pi^2/(2*L1)*2/L1;
               Ham0[4*(jc-1)+1:4*(jc-1)+4,4*(jc-1)+1:4*(jc-1)+4]+=(C1-D1*fac+D2*k^2)*Matrix{Float64}(I,4,4)+[(M1-B2*k^2+B1*fac)*sz A2*(kx-im*ky)*sx; A2*(kx+im*ky)*sx (M1-B2*k^2+B1*fac)*sz]
         end



         Ham0_per=Pmatrix'*Ham0*Pmatrix;
         FFF=eigen(Ham0_per[1:2*Ncut,1:2*Ncut]);
         basis=FFF.vectors[:,Ncut:Ncut+1];

         Ham0_proj=basis'*Ham0_per[1:2*Ncut,1:2*Ncut]*basis;




         #Constructing Hx now


        Hamx=zeros(ComplexF64,4*Ncut,4*Ncut);
        kx=0;
        ky=0;
        k=0;
         
        for jc in 1:Ncut
            fac=-jc^2*pi^2/(2*L1)*2/L1
            Hamx[4*(jc-1)+1:4*(jc-1)+4,4*(jc-1)+1:4*(jc-1)+4]+=[zeros(ComplexF64,2,2) A2*sx;A2*sx zeros(ComplexF64,2,2)];
        end
         Hamx_per=Pmatrix'*Hamx*Pmatrix
         Hamx_proj=basis'*Hamx_per[1:2*Ncut,1:2*Ncut]*basis

         #Constructing Hy now
         Hamy=zeros(ComplexF64,4*Ncut,4*Ncut);
            kx=0;
            ky=0;
            k=0;
         
         
         
                 for jc in 1:Ncut
                     fac=-jc^2*pi^2/(2*L1)*2/L1;
                     Hamy[4*(jc-1)+1:4*(jc-1)+4,4*(jc-1)+1:4*(jc-1)+4]+=[zeros(ComplexF64,2,2) -im*A2*sx; im*A2*sx zeros(ComplexF64,2,2)]; 
                 end

         Hamy_per=Pmatrix'*Hamy*Pmatrix;
         Hamy_proj=basis'*Hamy_per[1:2*Ncut,1:2*Ncut]*basis;


         #Constructing Hxx now



   Hamxx=zeros(ComplexF64,4*Ncut,4*Ncut);
      kx=0;
      ky=0;
      k=0;
         for jc in 1:Ncut
            fac=-jc^2*π^2/(2*L1)*2/L1;
            Hamxx[4*(jc-1)+1:4*(jc-1)+4,4*(jc-1)+1:4*(jc-1)+4]+=(2*D2)*Matrix{Float64}(I,4,4)+[(-2*B2)*sz zeros(Float64,2,2);zeros(Float64,2,2) (-2*B2)*sz];
        end

Hamxx_per=Pmatrix'*Hamxx*Pmatrix;
Hamxx_proj=basis'*Hamxx_per[1:2*Ncut,1:2*Ncut]*basis;
Hamyy_proj=Hamxx_proj;
HamVV=zeros(ComplexF64,4*Ncut,4*Ncut);
for ja in 1:Ncut, jb in 1:Ncut
    
     fac1=g^2*L1^2+(ja-jb)^2*pi^2;
     fac2=g^2*L1^2+(ja+jb)^2*pi^2;
     HamVV[4*(ja-1)+1:4*(ja-1)+4,4*(jb-1)+1:4*(jb-1)+4]+=2*Matrix{Float64}(I,4,4)*exp(-g*L1)*(exp(g*L1)-(-1)^(ja+jb))*g*L1^2*ja*jb*pi^2/(fac1*fac2)*2/L1;
 end

HamVV_per=Pmatrix'*HamVV*Pmatrix;
HamVV_proj=basis'*HamVV_per(1:2*Ncut,1:2*Ncut)*basis;
              
return Ham0_proj,Hamx_proj,Hamy_proj,Hamxx_proj,Hamyy_proj,HamVV_proj


end





function get_MoireHam_thinfilm(k::Vector{Float64},wave::Vector{Vector{Int64}},T1::Vector{Float64},T2::Vector{Float64},g1T::Vector{Int},g2T::Vector{Int},V0::Float64,ϕ::Float64,Ham0::Matrix{ComplexF64},Hamx::Matrix{ComplexF64},Hamy::Matrix{ComplexF64},Hamxx::Matrix{ComplexF64},Hamyy::Matrix{ComplexF64},HamV::Matrix{ComplexF64})
    dimension=2*length(wave)
    diagHam=zeros(ComplexF64,dimension,dimension)


    for ja in eachindex(wave)
       kvec=k+[T1 T2]*wave[ja]
       diagHam[2*(ja-1)+1:2*ja,2*(ja-1)+1:2*ja]=Ham0+kvec[1]*Hamx+kvec[2]*Hamy+kvec[1]^2/2*Hamxx+kvec[2]^2/2*Hamyy
    end

    MoirePo=zeros(ComplexF64,dimension,dimension)
    for ja in eachindex(wave)
    
       pos=findfirst(item->item==wave[ja]-g2T,wave)
       if pos≠nothing
           MoirePo[2*(ja-1)+1:2*(ja-1)+2,2*(pos-1)+1:2*(pos-1)+2]=V0*exp(im*ϕ)*HamV
       end
       
     
       pos=findfirst(item->item==wave[ja]-g1T,wave)
       if pos≠nothing
           MoirePo[2*(ja-1)+1:2*(ja-1)+2,2*(pos-1)+1:2*(pos-1)+2]=V0*exp(im*ϕ)*HamV
       end
   
       pos=findfirst(item->item==wave[ja]+g1T+g2T,wave)
       if pos≠nothing
           MoirePo[2*(ja-1)+1:2*(ja-1)+2,2*(pos-1)+1:2*(pos-1)+2]=V0*exp(im*ϕ)*HamV
       end
    end

    return diagHam+MoirePo+MoirePo'
 
end




function get_chern(Nq::Int64,T1::Vector{Float64},T2::Vector{Float64},chern_eigenvector::Array{ComplexF64})
    

    Uonelink=zeros(ComplexF64,Nq,Nq+1)
    Utwolink=zeros(ComplexF64,Nq+1,Nq)
    tra=zeros(Float64,Nq,Nq)
        
    for ja in 1:Nq, jb in 1:Nq+1
        Uonelink[ja,jb]=dot(chern_eigenvector[:,ja,jb],chern_eigenvector[:,ja+1,jb])/abs(dot(chern_eigenvector[:,ja,jb],chern_eigenvector[:,ja+1,jb]))
    end
    
        
        
    for ja in 1:Nq+1, jb in 1:Nq
        Utwolink[ja,jb]=dot(chern_eigenvector[:,ja,jb],chern_eigenvector[:,ja,jb+1])/abs(dot(chern_eigenvector[:,ja,jb],chern_eigenvector[:,ja,jb+1]))
    end
        
    Flink=zeros(ComplexF64,Nq,Nq)
    for ja in 1:Nq, jb in 1:Nq
     Flink[ja,jb]=log(Uonelink[ja,jb]*Utwolink[ja+1,jb]/(Uonelink[ja,jb+1]*Utwolink[ja,jb]))
    end
    chern=sum(Flink)/(2*π*im)
    
    aveF=sum(Flink)/Nq^2
    uniform=0.0
    for ja in 1:Nq, jb in 1:Nq
        uniform+=(imag(Flink[ja,jb])-imag(aveF))^2*Nq^2/(2π)^2
    end
    
    dG=norm(T2)
    
    for ja in 1:Nq, jb in 1:Nq
       
        A1=dot(chern_eigenvector[:,ja,jb],chern_eigenvector[:,ja,jb+1])
        B1=dot(chern_eigenvector[:,ja,jb],chern_eigenvector[:,ja+1,jb])
        C1=dot(chern_eigenvector[:,ja,jb],chern_eigenvector[:,ja+1,jb+1])
       gyy=(1-abs(A1)^2)/dG^2
       gxx=(2-1/2*gyy*dG^2-abs(B1)^2-abs(C1)^2)/(1.5*dG^2)
       tra[ja,jb]+=(gxx+gyy)*3^(1/2)/2*dG^2
    
      
    end
     
    trace_condition=sum(tra)-sum(abs.(Flink))

    return trace_condition,tra,Flink,chern,uniform

end

function main()
    L1=25.0
    Ham0,Hamx,Hamy,Hamxx,Hamyy,HamV=kpprojected_Ham(L1)
    HamV=Matrix{ComplexF64}(I,2,2)
    #V0space=collect(0.005:0.0005:0.015)
    V0space=collect(0.04:0.0005:0.06)
    am=210
    scale=4π/(√3*am)
    Nq=30
    g1=scale*[1,0]
    g2=scale*[-1/2,√3/2]
    T1=g1/Nq
    T2=g2/Nq
    
    chern_allowedq=Vector{Int64}[]
    for ja in 0:Nq, jb in 0:Nq
       push!(chern_allowedq,[ja,jb])
    end
    
    trace_condition=zeros(ComplexF64,length(V0space))
    chern=zeros(ComplexF64,length(V0space))
    uniform=zeros(ComplexF64,length(V0space))
    gapup=zeros(ComplexF64,length(V0space))
    gapdown=zeros(ComplexF64,length(V0space))
    direct_gapup=zeros(ComplexF64,length(V0space))
    direct_gapdown=zeros(ComplexF64,length(V0space))
    bandwidth=zeros(ComplexF64,length(V0space))

    
    ϕ=π/3
    g1T=Int.(round.(inv([T1 T2])*g1))
    g2T=Int.(round.(inv([T1 T2])*g2))
    
    
    
    wave=Vector{Int64}[]
    cutoff=18
    cutoffstandard=6.01*scale
    for ja in -cutoff:cutoff, jb in -cutoff:cutoff
        gtest=ja*g1+jb*g2;
        if (gtest[1]^2+gtest[2]^2)<cutoffstandard^2
            push!(wave,ja*g1T+jb*g2T)
        end
    end
    dimension=length(wave)*2
    bandindex=length(wave)+1

    Threads.@threads for jv in eachindex(V0space)
        V0=V0space[jv]
        chern_eigenvector=zeros(ComplexF64,dimension,Nq+1,Nq+1)
        eigenvalues=zeros(Float64,dimension,(Nq+1)^2)
     for ja in eachindex(chern_allowedq)
         k=[T1 T2]*chern_allowedq[ja]
         MoireHam=get_MoireHam_thinfilm(k,wave,T1,T2,g1T,g2T,V0,ϕ,Ham0,Hamx,Hamy,Hamxx,Hamyy,HamV)
         FFF=eigen(MoireHam)
         chern_eigenvector[:,chern_allowedq[ja][1]+1,chern_allowedq[ja][2]+1]=FFF.vectors[:,bandindex]
         eigenvalues[:,ja]=real.(FFF.values)
     end
     trace_condition[jv],_,_,chern[jv],uniform[jv]=get_chern(Nq,T1,T2,chern_eigenvector)
     gapup[jv]=sort(eigenvalues[bandindex+1,:])[1]-sort(eigenvalues[bandindex,:])[(Nq+1)^2]
     gapdown[jv]=sort(eigenvalues[bandindex,:])[1]-sort(eigenvalues[bandindex-1,:])[(Nq+1)^2]
     bandwidth[jv]=-sort(eigenvalues[bandindex,:])[1]+sort(eigenvalues[bandindex,:])[(Nq+1)^2]
     direct_gapup[jv]=sort(eigenvalues[bandindex+1,:]-eigenvalues[bandindex,:])[1]
     direct_gapdown[jv]=sort(eigenvalues[bandindex,:]-eigenvalues[bandindex-1,:])[1]
   end


   return V0space,am,L1,trace_condition,chern,uniform,gapup,gapdown,direct_gapup,direct_gapdown,bandwidth

end