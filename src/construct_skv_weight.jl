using EllipticFunctions
using LinearAlgebra
using Plots
using JLD2,StaticArrays
using Combinatorics
using Random


function c_dot(z1::ComplexF64,z2::ComplexF64)
    return 1/2*(z1*conj(z2)+z2*conj(z1))
end
function c_cross(z1::ComplexF64,z2::ComplexF64)
   return 1/2*im*(z1*conj(z2)-z2*conj(z1))
end

function wrapparallel(k::Complex, b1::Complex, b2::Complex)
    A = [real(b1) real(b2);
         imag(b1) imag(b2)]
    rhs = [real(k), imag(k)]
    uv = A \ rhs
    u, v = uv[1], uv[2]

    frac(t) = t - round(t)   # matches Mathematica t - Round[t]
    return frac(u)*b1 + frac(v)*b2
end


function wrap_n1n2(k::Complex, b1::Complex, b2::Complex)
    A = [real(b1) real(b2);
         imag(b1) imag(b2)]
    rhs = [real(k), imag(k)]
    uv = A \ rhs
    u, v = uv[1], uv[2]

    frac(t) = t - round(t)   # matches Mathematica t - Round[t]
    return frac(u)*b1 + frac(v)*b2, Int(round(u)), Int(round(v))
end



function get_spinor_norm(NL::Int,q::Vector{Float64})
    qh=q[1]+im*q[2]
    aGr=0.246
    t0=3100
    t1=380
    vf=√3/2*aGr*t0
    spinor=zeros(ComplexF64,NL)
    for ja in 1:NL
      spinor[ja]=qh^(ja-1)*(vf/t1)^(ja-1)
    end
  return 1/norm(spinor)
end


function get_spinor(NL::Int,q::Vector{Float64})
    qh=q[1]+im*q[2]
    aGr=0.246
    t0=3100
    t1=380
    vf=√3/2*aGr*t0
    spinor=zeros(ComplexF64,NL)
    for ja in 1:NL
      spinor[ja]=qh^(ja-1)*(vf/t1)^(ja-1)
    end
  return spinor/norm(spinor)
end


struct sigma_related
    pv::ComplexF64
    ev::ComplexF64
end


function torusSigma(z,L1,L2,Lb,bigA)
    w_z,n1,n2=wrap_n1n2(z,L1/Lb,L2/Lb)

    pv=(-1)^(n1+n2)*wsigma(w_z, omega=(L1/(2*Lb), L2/(2*Lb)))*exp(-bigA*w_z^2)
    ev=(L1*conj(L1)*n1^2+L2*conj(L2)*n2^2)/(4*Lb^2)+(conj(L2)*n2+conj(L1)*n1)/(2*Lb)*w_z+n1*n2*conj(L2)*L1/(2*Lb^2)
    
    return sigma_related(pv,ev)

end


function modsigma(z,a1,a2,lb,α)
    w_z,n1,n2=wrap_n1n2(z,a1/lb,a2/lb)
    pv=(-1)^(n1+n2)*wsigma(w_z, omega=(a1/(2*lb), a2/(2*lb)))*exp(-α*w_z^2)
    ev=(a1*conj(a1)*n1^2+a2*conj(a2)*n2^2)/(4*lb^2)+(conj(a2)*n2+conj(a1)*n1)/(2*lb)*w_z+n1*n2*conj(a2)*a1/(2*lb^2)

    return sigma_related(pv,ev)
end


function modsigma_f(z,a1f,a2f,lbf,αf)
    w_z,n1,n2=wrap_n1n2(z,a1f/lbf,a2f/lbf)
    pv=(-1)^(n1+n2)*wsigma(w_z, omega=(a1f/(2*lbf), a2f/(2*lbf)))*exp(-αf*w_z^2)
    ev=(a1f*conj(a1f)*n1^2+a2f*conj(a2f)*n2^2)/(4*lbf^2)+(conj(a2f)*n2+conj(a1f)*n1)/(2*lbf)*w_z+n1*n2*conj(a2f)*a1f/(2*lbf^2)

    return sigma_related(pv,ev)
end


function landaulevel(z, k, a1, a2, lb, α, b1, b2)
    ss = modsigma(
        z/lb + im*lb*(k - b1/2 - b2/2),
        a1, a2, lb, α
    )

    ev_add =
        im/2 * conj(k - b1/2 - b2/2) * z -
        z*conj(z)/(4*lb^2) -
        lb^2/4 * conj(k)*k +
        lb^2*k*conj(b1+b2)/4

    # Normalize so that
    #
    # (1 / area_cell) ∫cell d²r |ϕ_k(r)|² = 1.
    #
    # Equivalently, ∫cell d²r |ϕ_k(r)|² = area_cell.
    τ = a2 / a1
    ηD = etaDedekind(τ)

    @assert imag(τ) > 0

    log_NLLL =
        0.5 * (
            log(2π) +
            0.5 * log(2.0) +
            1.5 * log(imag(τ)) +
            6 * log(abs(ηD)) -
            lb^2 * abs2(b1 + b2) / 8
        )

    return sigma_related(
        ss.pv,
        ss.ev + ev_add + log_NLLL
    )
end


function landaulevel_f(z, k, a1f, a2f, lbf, αf, b1f, b2f)
    ss = modsigma_f(
        z/lbf + im*lbf*(k - b1f/2 - b2f/2),
        a1f, a2f, lbf, αf
    )

    ev_add =
        im/2 * conj(k - b1f/2 - b2f/2) * z -
        z*conj(z)/(4*lbf^2) -
        lbf^2/4 * conj(k)*k +
        lbf^2*k*conj(b1f+b2f)/4

    # Normalize so that
    #
    # (1 / area_cell_f) ∫cell_f d²r |ϕ_k(r)|² = 1.
    τf = a2f / a1f
    ηDf = etaDedekind(τf)

    @assert imag(τf) > 0

    log_NLLLf =
        0.5 * (
            log(2π) +
            0.5 * log(2.0) +
            1.5 * log(imag(τf)) +
            6 * log(abs(ηDf)) -
            lbf^2 * abs2(b1f + b2f) / 8
        )

    return sigma_related(
        ss.pv,
        ss.ev + ev_add + log_NLLLf
    )
end

function test_Landau_level_sigma(L1,L2,Lb,bigA
                                 ,a1,a2,lb,α,b1,b2,
                                 a1f,a2f,lbf,αf,b1f,b2f)

  
   
    for attempt in 1:10
        z_test=rand(ComplexF64)*abs(L1)
        t1=torusSigma(z_test/Lb,L1,L2,Lb,bigA)
        t2=torusSigma(z_test/Lb+2*L1/Lb+L2/Lb,L1,L2,Lb,bigA)
        if abs(t1.pv * exp(t1.ev)) ≥ 1e-4
          @assert abs(t2.pv/t1.pv*exp(t2.ev-t1.ev)/((-1)^(2+1)*exp(L1*conj(L1)*2^2/(4*Lb^2)+L2*conj(L2)*1^2/(4*Lb^2)+(conj(L2)*1+conj(L1)*2)/(2*Lb)*z_test/Lb+2*1*conj(L2)*L1/(2*Lb^2)))-1)<10^(-9)
            break
        end

        if attempt == 10
            error("Failed to find a valid sample after 10 attempts")
        end
    end





    for attempt in 1:10
        z_test = rand(ComplexF64)*abs(L1)
        k_test = rand(ComplexF64)*abs(L1)

        t1 = landaulevel(z_test, k_test,a1,a2,lb,α,b1,b2)
        t2 = landaulevel(z_test + 5*a1 + 3*a2, k_test,a1,a2,lb,α,b1,b2)

        if abs(t1.pv * exp(t1.ev)) ≥ 1e-4
            @assert abs(t2.pv/t1.pv*exp(t2.ev-t1.ev)-exp(im*c_dot(k_test,5*a1+3*a2))*exp(-im*c_cross(z_test,5*a1+3*a2)/(2*lb^2))*(-1)^(3*5))<10^(-8)

            break
        end

        if attempt == 10
            error("Failed to find a valid sample after 10 attempts")
        end
    end




    for attempt in 1:10
        z_test = rand(ComplexF64)*abs(L1)
        k_test = rand(ComplexF64)*abs(L1)

        t1=landaulevel(z_test,k_test,a1,a2,lb,α,b1,b2)
        t2=landaulevel(z_test,k_test+3*b1+5*b2,a1,a2,lb,α,b1,b2)

        if abs(t1.pv * exp(t1.ev)) ≥ 1e-4
           @assert abs(t2.pv/t1.pv*exp(t2.ev-t1.ev)-exp(-im*lb^2*c_cross(k_test,3*b1+5*b2)/2)*(-1)^(3*5))<10^(-8)

            break
        end

        if attempt == 10
            error("Failed to find a valid sample after 10 attempts")
        end
    end


    
    for attempt in 1:10
        z_test = rand(ComplexF64)*abs(L1)
        k_test = rand(ComplexF64)*abs(L1)

        t1=landaulevel_f(z_test,k_test,a1f,a2f,lbf,αf,b1f,b2f)
        t2=landaulevel_f(z_test+5*a1f+3*a2f,k_test,a1f,a2f,lbf,αf,b1f,b2f)

        if abs(t1.pv * exp(t1.ev)) ≥ 1e-4
           @assert abs(t2.pv/t1.pv*exp(t2.ev-t1.ev)-exp(im*c_dot(k_test,5*a1f+3*a2f))*exp(-im*c_cross(z_test,5*a1f+3*a2f)/(2*lbf^2))*(-1)^(3*5))<10^(-8)

            break
        end

        if attempt == 10
            error("Failed to find a valid sample after 10 attempts")
        end
    end



     
    for attempt in 1:10
           
            z_test=rand(ComplexF64)*abs(L1)
            k_test=rand(ComplexF64)*abs(L1)
            t1=landaulevel_f(z_test,k_test,a1f,a2f,lbf,αf,b1f,b2f)
            t2=landaulevel_f(z_test,k_test+3*b1f+5*b2f,a1f,a2f,lbf,αf,b1f,b2f)

        if abs(t1.pv * exp(t1.ev)) ≥ 1e-4
           @assert abs(t2.pv/t1.pv*exp(t2.ev-t1.ev)-exp(-im*lbf^2*c_cross(k_test,3*b1f+5*b2f)/2)*(-1)^(3*5))<10^(-8)

            break
        end

        if attempt == 10
            error("Failed to find a valid sample after 10 attempts")
        end
    end
    


   



end




function F1_withSigma(Z::ComplexF64,kappa::ComplexF64,
           eta_set::Vector{ComplexF64},Nelectron::Int,xi_set::Vector{ComplexF64},
           T1::ComplexF64,T2::ComplexF64,
           L1::ComplexF64,L2::ComplexF64,
           Lb::Float64,bigA::ComplexF64,lbf::Float64,Deltatheta::ComplexF64) # I use label kappa rather than minus kappa, I put the Sigma in here
    zarg=Z/Lb+2*im*Lb*(-kappa/2+Nelectron/4*(T1+T2))
    ss=torusSigma(zarg,L1,L2,Lb,bigA)
    ev_Sigma=-sum(eta_set)*conj(sum(eta_set))/(4*Lb^2)+sum(conj(eta_set))/(2*Lb^2)*sum(xi_set)+sum(eta_set.*conj(eta_set))/(4*lbf^2)+Deltatheta*im/2*conj(sum(eta_set))


    return sigma_related(ss.pv,ss.ev+im*Z*(-conj(kappa)/2+Nelectron/4*(conj(T1)+conj(T2)))+ev_Sigma)
end



function evalute_determinant_part_product_form(z_pos::Vector{ComplexF64},kappa::ComplexF64,eta_set::Vector{ComplexF64},xi_set::Vector{ComplexF64},
                                              lbf::Float64,Deltatheta::ComplexF64,
                                              L1::ComplexF64,L2::ComplexF64,Lb::Float64,bigA::ComplexF64,
                                              Nelectron::Int,T1::ComplexF64,T2::ComplexF64)  
   
   
   
    area_torus =0.5 * abs(L1 * conj(L2) - L2 * conj(L1))
    ev_ini=-sum(z_pos.*conj(z_pos))/(4*lbf^2)+Nelectron/2 * log(area_torus)
   pv_ini=1.0+0.0*im
   for p1 in 1:length(z_pos)
     for p2 in p1+1:length(z_pos)
         ss=torusSigma((z_pos[p1]-z_pos[p2])/Lb,L1,L2,Lb,bigA)
         ev_ini+=ss.ev
         pv_ini*=ss.pv
     end
    end


    sf=F1_withSigma(sum(z_pos),kappa,eta_set,Nelectron,xi_set,T1,T2,L1,L2,Lb,bigA,lbf,Deltatheta)
     ev_ini+=sf.ev
     pv_ini*=sf.pv


    return sigma_related(pv_ini,ev_ini)

end



function evalute_determinant_part_determinant_form(z_pos::Vector{ComplexF64},kappa::ComplexF64,
                                                  kf_set::Vector{ComplexF64},
                                                  a1f::ComplexF64,a2f::ComplexF64,lbf::Float64,αf::ComplexF64,
                                                  b1f::ComplexF64,b2f::ComplexF64)  
   MM=zeros(ComplexF64,length(kf_set),length(z_pos))

   for kk in eachindex(kf_set), pp in eachindex(z_pos)
      ss=landaulevel_f(z_pos[pp],kf_set[kk]-kappa,a1f,a2f,lbf,αf,b1f,b2f)
      MM[kk,pp]=ss.pv*exp(ss.ev)
   end

   logabs, phase = logabsdet(MM)
   log_sqrt_factorial =
        0.5 * sum(log, 1:Nelectron)

    return sigma_related(1.0+0.0*im,logabs+im*angle(phase) - log_sqrt_factorial)

end









function evalute_phi_vac(z_pos::ComplexF64,xi_set::Vector{ComplexF64},
                    eta_set::Vector{ComplexF64},xi_set_partial::Vector{ComplexF64},
                    shift_set::Vector{Vector{Int64}},
                    L1::ComplexF64,L2::ComplexF64,Lb::Float64,bigA::ComplexF64,
                    lbf::Float64)
    ev_ini=0.0+0.0*im
    pv_ini=1.0+0.0*im


    ev_ini+=-z_pos*conj(z_pos)/(4*lbf^2)


    ev_ini+=(sum(xi_set)-sum(eta_set))/(2*Lb^2)*conj(z_pos)



    for xx in eachindex(xi_set_partial)
            ss=torusSigma((z_pos-xi_set_partial[xx])/Lb,L1,L2,Lb,bigA)
            ev_ini+=conj(ss.ev)
            pv_ini*=conj(ss.pv)
    end

    

    ev_intermidiate=0.0
    for xx in eachindex(shift_set)
            n11=shift_set[xx][1]
            n22=shift_set[xx][2]
            zz=z_pos-eta_set[xx]
             
            pv_ini*=(-1)^(n11+n22)
            ev_intermidiate+=(-L1*conj(L1)*n11^2/(4*Lb^2)-L2*conj(L2)*n22^2/(4*Lb^2)-n11*n22*conj(L2)*L1/(2*Lb^2)-conj(L2)*n22/(2*Lb^2)*zz-conj(L1)*n11/(2*Lb^2)*zz)
    end
    ev_ini+=-conj(ev_intermidiate)
    
    
    
    
    
    return sigma_related(pv_ini,ev_ini)



end

#=
function testBlochF(z::ComplexF64,k1::ComplexF64,k2::ComplexF64,grid_cutoff::Float64,
                b1f::ComplexF64,b2f::ComplexF64,lbf::Float64,Deltatheta::ComplexF64,kappa::ComplexF64,
                b1::ComplexF64,b2::ComplexF64)
    
    Nmax=Int(round.(max(20,grid_cutoff*abs(b1)/abs(b1f)*2,grid_cutoff*abs(b1)/abs(b2f)*2)))
   
    
    dd=0.0
    R0=-im*lbf^2*(k2+Deltatheta-kappa)
    for ja in -Nmax:Nmax, jb in -Nmax:Nmax
      g=ja*b1f+jb*b2f
      dd+=exp(im*π*ja*jb)*exp(-lbf^2/4*(g*conj(g)+2*(k1-k2-Deltatheta)*conj(g)))*exp(im*c_dot(z-R0,k1-k2+g-Deltatheta))
    end
    return dd
end
=#


function testBlochF(
    z::ComplexF64,
    k1::ComplexF64,
    k2::ComplexF64,
    k1_n::Vector{Int},
    k2_n::Vector{Int},
    Tgrid::Vector{Vector{Int}},
    BfT_matrix,
    T1::ComplexF64,
    T2::ComplexF64,
    lbf::Float64,
    Deltatheta::ComplexF64,
    kappa::ComplexF64,
)

    dd = 0.0 + 0.0im

    q = k1 - k2 - Deltatheta
    R0 = -im*lbf^2*(k2 + Deltatheta - kappa)

    delta_k_n = k1_n - k2_n

    for t_n in Tgrid
        # t_n = k1_n - k2_n + g_n
        g_n = t_n - delta_k_n

        g_coeff = BfT_matrix \ g_n

        n1g = round(Int, g_coeff[1])
        n2g = round(Int, g_coeff[2])

        if maximum(abs.(g_coeff .- [n1g, n2g])) > 1e-10
            continue
        end

        gg = g_n[1]*T1 + g_n[2]*T2
        Q = q + gg

        # The momentum-pair envelope exp(+lbf²|q|²/4)
        # has been removed.
        ct_2 =
            exp(im*π*n1g*n2g) *
            exp(
                -lbf^2/4 *
                (
                    abs2(gg) +
                    2*q*conj(gg) +
                    abs2(q)
                )
            )

        ct_3 = exp(im*c_dot(z - R0, Q))

        dd += ct_2*ct_3
    end

    return dd
end





function evalute_vacancy_part_withSigma(
                    eta_set::Vector{ComplexF64},xi_set::Vector{ComplexF64},
                    topo_sec::Int,
                    L1::ComplexF64,L2::ComplexF64,Lb::Float64,bigA::ComplexF64,
                    Deltatheta::ComplexF64,lbf::Float64,
                    Nelectron::Int,T1::ComplexF64,T2::ComplexF64,type::Int)

    @assert in(topo_sec,[0,1,2])
                    

    ev_ini=0.0+0.0*im
    pv_ini=1.0+0.0*im


    if type==1
            ev_ini+=-sum(eta_set)*conj(sum(eta_set))/(2*Lb^2)

            for p1 in 1:length(eta_set)
            for p2 in p1+1:length(eta_set)
                ss=torusSigma((eta_set[p1]-eta_set[p2])/Lb,L1,L2,Lb,bigA)
                ev_ini+=conj(ss.ev)*2
                pv_ini*=conj(ss.pv)^2
        
            end
            end
            
            alphas=Lb^2/3*(Nelectron*(T1+T2)/2+topo_sec*T1)+Lb^2/3*Deltatheta
            betas=Deltatheta/2-im*sum(xi_set)/(2*Lb^2)+3*alphas/(2*Lb^2)-Deltatheta/2

            ev_ini+=im*betas*conj(sum(eta_set))
            ev_ini-=L1/(2*Lb^2)*conj(sum(eta_set))
            ev_ini-=im*sum(eta_set)*conj(Deltatheta)/2

            for jj in 0:2
            ll=torusSigma((sum(eta_set)-im*alphas+jj*L1/3)/Lb,L1,L2,Lb,bigA)
            
            ev_ini+=conj(ll.ev)
            pv_ini*=conj(ll.pv)
            
            
            end

            ev_ini-=-sum(eta_set)*conj(sum(eta_set))/(4*Lb^2)+sum(conj(eta_set))/(2*Lb^2)*sum(xi_set)+sum(eta_set.*conj(eta_set))/(4*lbf^2)+Deltatheta*im/2*conj(sum(eta_set))
            
            
    
             return sigma_related(pv_ini,ev_ini)
    elseif type==2
            ev_ini+=-sum(eta_set.*conj.(eta_set))/(4*lbf^2)

            for p1 in 1:length(eta_set)
            for p2 in p1+1:length(eta_set)
                ss=torusSigma((eta_set[p1]-eta_set[p2])/Lb,L1,L2,Lb,bigA)
                ev_ini+=conj(ss.ev)*3+ss.ev
                pv_ini*=conj(ss.pv)^3*ss.pv
        
            end
            end
            
            alphas=Lb^2/3*(Nelectron*(T1+T2)/2+topo_sec*T1)+Lb^2/3*Deltatheta
            betas=Deltatheta/2-im*sum(xi_set)/(2*Lb^2)+3*alphas/(2*Lb^2)-Deltatheta/2

            ev_ini+=im*betas*conj(sum(eta_set))
            ev_ini-=L1/(2*Lb^2)*conj(sum(eta_set))
            ev_ini-=im*sum(eta_set)*conj(Deltatheta)/2

            for jj in 0:2
            ll=torusSigma((sum(eta_set)-im*alphas+jj*L1/3)/Lb,L1,L2,Lb,bigA)
            
            ev_ini+=conj(ll.ev)
            pv_ini*=conj(ll.pv)
            
            
            end

            ev_ini-=-sum(eta_set)*conj(sum(eta_set))/(4*Lb^2)+sum(conj(eta_set))/(2*Lb^2)*sum(xi_set)+sum(eta_set.*conj(eta_set))/(4*lbf^2)+Deltatheta*im/2*conj(sum(eta_set))
            
            
    
             return sigma_related(pv_ini,ev_ini)


    else

        error("type not implemented")
    end

end







function evalute_full_formula_withSigma(z_pos::Vector{ComplexF64},kappa::ComplexF64,xi_set::Vector{ComplexF64},
                    eta_set::Vector{ComplexF64},xi_set_partial::Vector{ComplexF64},
                    shift_set::Vector{Vector{Int64}},
                    L1::ComplexF64,L2::ComplexF64,Lb::Float64,bigA::ComplexF64,
                    Deltatheta::ComplexF64,lbf::Float64,
                    Nelectron::Int,T1::ComplexF64,T2::ComplexF64)

    ev_ini=0.0+0.0*im
    pv_ini=1.0+0.0*im

    ev_ini+=-sum(z_pos.*conj(z_pos))/(2*lbf^2)

    for p1 in 1:length(z_pos)
     for p2 in p1+1:length(z_pos)
         ss=torusSigma((z_pos[p1]-z_pos[p2])/Lb,L1,L2,Lb,bigA)
         ev_ini+=ss.ev
         pv_ini*=ss.pv
     end
    end


    sf=F1_withSigma(sum(z_pos),kappa,eta_set,Nelectron,xi_set,T1,T2,L1,L2,Lb,bigA,lbf,Deltatheta)
    ev_ini+=sf.ev
    pv_ini*=sf.pv



    ev_ini+=(sum(xi_set)-sum(eta_set))/(2*Lb^2)*conj(sum(z_pos))


    
    for  p1 in 1:length(z_pos)
        for xx in eachindex(xi_set_partial)
            ss=torusSigma((z_pos[p1]-xi_set_partial[xx])/Lb,L1,L2,Lb,bigA)
            ev_ini+=conj(ss.ev)
            pv_ini*=conj(ss.pv)
        end
    end
    
    for p1 in 1:length(z_pos)
        ev_intermidiate=0.0
        for xx in eachindex(shift_set)
            n11=shift_set[xx][1]
            n22=shift_set[xx][2]
            zz=z_pos[p1]-eta_set[xx]
             
            pv_ini*=(-1)^(n11+n22)
            ev_intermidiate+=(-L1*conj(L1)*n11^2/(4*Lb^2)-L2*conj(L2)*n22^2/(4*Lb^2)-n11*n22*conj(L2)*L1/(2*Lb^2)-conj(L2)*n22/(2*Lb^2)*zz-conj(L1)*n11/(2*Lb^2)*zz)
        end
        ev_ini+=-conj(ev_intermidiate)
    end
    
    
    
    
    #ev_ini+=-sum(eta_set)*conj(sum(eta_set))/(4*Lb^2)+sum(conj(eta_set))/(2*Lb^2)*sum(xi_set)+sum(eta_set.*conj(eta_set))/(4*lbf^2)+Deltatheta*im/2*conj(sum(eta_set))


    return sigma_related(pv_ini,ev_ini)
end





function fixed_configuration_result(eta_set::Vector{ComplexF64},shift_set::Vector{Vector{Int}},xi_set_partial::Vector{ComplexF64},params::NamedTuple)

   
   
                a1=params.a1
                a2=params.a2
                lb=params.lb
                b1=params.b1
                b2=params.b2
                α=params.α

                L1=params.L1
                L2=params.L2
                Lb=params.Lb
                T1=params.T1
                T2=params.T2
                bigA=params.bigA

                a1f=params.a1f
                a2f=params.a2f
                lbf=params.lbf
                b1f=params.b1f
                b2f=params.b2f
                αf=params.αf

                N1=params.N1
                N2=params.N2
                N1f=params.N1f
                N2f=params.N2f
                Nvac=params.Nvac
                Nphi=params.Nphi
                Nelectron=params.Nelectron

                Deltatheta=params.Deltatheta
               
                T1_vec=params.T1_vec
                T2_vec=params.T2_vec
                b1_vec=params.b1_vec
                b2_vec=params.b2_vec
                b1f_vec=params.b1f_vec
                b2f_vec=params.b2f_vec
                b1T=params.b1T
                b2T=params.b2T
                b1fT=params.b1fT
                b2fT=params.b2fT
                BfT_matrix = hcat(b1fT, b2fT)

                a1_vec=params.a1_vec
                a2_vec=params.a2_vec
                L1_vec=params.L1_vec
                L2_vec=params.L2_vec

                Tgrid=params.Tgrid
                Tgrid_dict=params.Tgrid_dict
                gfgrid=params.gfgrid
                gfgrid_dict=params.gfgrid_dict

                kf_set=params.kf_set
                kf_set_n=params.kf_set_n

                xi_set=params.xi_set

                topo_sec=params.topo_sec
                NL=params.NL
                grid_cutoff=  params.grid_cutoff
                type=params.type
   
   
   
   
   
   
   
    kappa=2*(im*sum(eta_set)-im*sum(xi_set))/(2*Lb^2)+Nelectron*(T1+T2)/2+Deltatheta


 #getting M_eta
    M_eta=0.0  
    for attempt in 1:10
      
        z_pos=wrapparallel.(rand(ComplexF64,N1f*N2f)*3*abs(L1),L1,L2)

        s1=evalute_determinant_part_product_form(z_pos,kappa,eta_set,xi_set,
                                                    lbf,Deltatheta,
                                                    L1,L2,Lb,bigA,
                                                    Nelectron,T1,T2)

        s2=evalute_determinant_part_determinant_form(z_pos,kappa,kf_set,
                                                        a1f,a2f,lbf,αf,
                                                        b1f,b2f)




    println(s1.pv * exp(s1.ev),"s1")
    println(s2.pv * exp(s2.ev),"s2")
        
        if abs(s1.pv * exp(s1.ev)) ≥ 10^(-8) && abs(s2.pv * exp(s2.ev)) ≥ 10^(-8) && abs(s1.pv * exp(s1.ev)) < 10^(18) && abs(s2.pv * exp(s2.ev)) < 10^(18)
        
            M_eta=sigma_related(s1.pv/s2.pv,s1.ev-s2.ev)
            break

           
        end

        if attempt == 10
            error("Failed to find a valid sample for Meta after 10 attempts")
        end
    end




        Nmatrix=zeros(ComplexF64,length(xi_set_partial),length(kf_set))
        for ja in eachindex(xi_set_partial), jb in eachindex(kf_set)
        ss=landaulevel_f(xi_set_partial[ja],kf_set[jb]+Deltatheta-kappa,
                            a1f,a2f,lbf,αf,b1f,b2f)
        Nmatrix[ja,jb]=conj(ss.pv*exp(ss.ev))
        end


        F = svd(Nmatrix)          # Mmatrix = U * Diagonal(S) * V'
        σmin = F.S[end]           # smallest singular value

        nullvec = F.V[:, end]        # right singular vector for σmin (unit-norm)
        # (svd already returns V with orthonormal columns; normalize again if you want)
        nullvec ./= norm(nullvec)

        # quick sanity check: ‖M vmin‖ ≈ σmin
        res = norm(Nmatrix * nullvec)
        if abs(res)>10^(-9)
            error("This singular value is not zero")
        end


        N_coeff=copy(nullvec)

  

        try_p1=xi_set_partial[1]+1/2*(a1+a2)
        phi_vac_p1=evalute_phi_vac(try_p1,xi_set,
                            eta_set,xi_set_partial,shift_set,
                            L1,L2,Lb,bigA,
                            lbf)
        sum_p1=0.0
        for ja in eachindex(kf_set)
            ss=landaulevel_f(try_p1,kf_set[ja]+Deltatheta-kappa,
                            a1f,a2f,lbf,αf,b1f,b2f)
            sum_p1+=N_coeff[ja]*conj(ss.pv*exp(ss.ev))
        end
        if abs(phi_vac_p1.pv*exp(phi_vac_p1.ev)/sum_p1)>10^(10) || abs(phi_vac_p1.pv*exp(phi_vac_p1.ev)/sum_p1)<10^(-10)
            error("too large for N_coeff, Consider using alternative form")
        else

           N_coeff*=(phi_vac_p1.pv*exp(phi_vac_p1.ev)/sum_p1)
        end


        #Mkkmatrix=zeros(ComplexF64,length(kf_set),length(kf_set))
        Mkkmatrix=Matrix{sigma_related}(undef,length(kf_set),length(kf_set))
        for j1 in eachindex(kf_set), j2 in eachindex(kf_set)
            
            k11=kf_set[j1]
            k22=kf_set[j2]
            k11_n = kf_set_n[j1]
            k22_n = kf_set_n[j2]
            z_test=rand(ComplexF64)*abs(L1)
            itcount=0
            while abs(wrapparallel(z_test+im*lbf^2*(k11-kappa-b1f/2-b2f/2),a1f,a2f))<abs(1/8*a1f) || abs(wrapparallel(z_test+im*lbf^2*(k22-kappa+Deltatheta-b1f/2-b2f/2),a1f,a2f))<abs(1/8*a1f) 
                z_test=rand(ComplexF64)*abs(L1)
                itcount+=1
                if itcount==10
                    error("sampling error in calculating MKK")
                end
            end
            
                ll1=landaulevel_f(z_test,k11-kappa,a1f,a2f,lbf,αf,b1f,b2f)
                ll2=landaulevel_f(z_test,k22-kappa+Deltatheta,a1f,a2f,lbf,αf,b1f,b2f)
     
                BlochF = testBlochF(
                z_test,
                k11,
                k22,
                k11_n,
                k22_n,
                Tgrid,
                BfT_matrix,
                T1,
                T2,
                lbf,
                Deltatheta,
                kappa,
                 )
                println(abs(BlochF))
                @assert 10^(12)>abs(BlochF)>10^(-10)

                Mkkmatrix[j1,j2]=sigma_related(ll1.pv*conj(ll2.pv),ll1.ev+conj(ll2.ev)-log(abs(BlochF))-im*angle(BlochF))
          


        end


        orbital_basis = zeros(
        ComplexF64,
        length(Tgrid),
        Nelectron,
        )

    for k1i in eachindex(kf_set)

        k1_n = kf_set_n[k1i]
        k1 = kf_set[k1i]

        for k2i in eachindex(kf_set)

            k2_n = kf_set_n[k2i]
            k2 = kf_set[k2i]

            Mkk = Mkkmatrix[k1i, k2i]

            ct_1 =
                Mkk.pv *
                exp(Mkk.ev) *
                N_coeff[k2i]

            delta_k_n = k1_n - k2_n
            q = k1 - k2 - Deltatheta

            Rk2 =
                -im * lbf^2 *
                (
                    k2 +
                    Deltatheta -
                    kappa
                )

            for pp in eachindex(Tgrid)

                t_n = Tgrid[pp]

                # t_n = k1_n - k2_n + g_n
                g_n = t_n - delta_k_n

                # Check whether g_n belongs to the reciprocal lattice
                # generated by b1fT and b2fT.
                g_coeff = BfT_matrix \ g_n

                n1g = round(Int, g_coeff[1])
                n2g = round(Int, g_coeff[2])

                if maximum(
                    abs.(g_coeff .- [n1g, n2g])
                ) > 1e-10
                    continue
                end

                gg =
                    g_n[1] * T1 +
                    g_n[2] * T2

                # Physical Fourier momentum retained in Tgrid.
                Q =
                    t_n[1] * T1 +
                    t_n[2] * T2 -
                    Deltatheta

                @assert abs(Q - (q + gg)) < 1e-9

                # The envelope exp(+lbf^2*abs2(q)/4) has been removed.
                # Consequently, abs(ct_2) decays as
                # exp(-lbf^2*abs2(Q)/4).
                ct_2 =
                    exp(im * π * n1g * n2g) *
                    exp(
                        -lbf^2 / 4 *
                        (
                            abs2(gg) +
                            2 * q * conj(gg) +
                            abs2(q)
                        )
                    )

                ct_3 =
                    exp(
                        -im * c_dot(Q, Rk2)
                    )

                ct_4 =
                    1 / get_spinor_norm(
                        NL,
                        [real(Q), imag(Q)],
                    )

                orbital_basis[pp, k1i] +=
                    ct_1 * ct_2 * ct_3 * ct_4
            end
        end
    end


 
        orbital_basis_norm=ones(Float64,length(kf_set))
        for ja in 1:Nelectron
            nn=norm(orbital_basis[:,ja])
            orbital_basis[:,ja]=orbital_basis[:,ja]/nn
        
            orbital_basis_norm[ja]=nn
        end
   


        F = qr(orbital_basis)                 # QR factorization
        orbital_basis_orthogonal= Matrix(F.Q)[:, 1:size(orbital_basis,2)]
        orbital_Rmatrix=Matrix(F.R)

       
        vacancy_factor=evalute_vacancy_part_withSigma(eta_set,xi_set,
                    topo_sec,
                    L1,L2,Lb,bigA,
                    Deltatheta,lbf,
                    Nelectron,T1,T2,type)




  return vacancy_factor,orbital_Rmatrix,orbital_basis_norm,orbital_basis_orthogonal,M_eta,Mkkmatrix,N_coeff

end




function big_func(flux1::Float64,flux2::Float64,NL::Int,moiream::Float64,
                  N1::Int,N2::Int,N1f::Int,N2f::Int,xi_00::ComplexF64,
                  grid_cutoff::Float64,topo_sec::Int,type::Int,sample_num::Int)

        θ = π/3;
        
        l1 = moiream;
        l2 = moiream;

        a1 = ComplexF64(l1)
        a2 = l2*(cos(θ) + im*sin(θ))
        lb=(1/(2π)*1/2*abs(a1*conj(a2)-a2*conj(a1)))^(1/2)

        b1=-im*a2/lb^2
        b2=im*a1/lb^2

        
        Nphi=N1*N2
      
        Nvac=Nphi-N1f*N2f
        Nelectron=N1f*N2f

        @assert Int(N1*N2*2/3)==N1f*N2f


        Lb=(N1*N2)^(1/2)*lb
        L1=N1*a1
        L2=N2*a2
        lbf=Lb/(N1f*N2f)^(1/2)

        a1f=L1/N1f
        a2f=L2/N2f
        b1f=-im*a2f/lbf^2
        b2f=im*a1f/lbf^2
        T1=b1/N1
        T2=b2/N2
    
        Deltatheta=(T1*flux1+T2*flux2)


        T1_vec=[real(T1),imag(T1)]
        T2_vec=[real(T2),imag(T2)]
        b1_vec=[real(b1),imag(b1)]
        b2_vec=[real(b2),imag(b2)]
        b1f_vec=[real(b1f),imag(b1f)]
        b2f_vec=[real(b2f),imag(b2f)]
        b1T=Int.(round.(inv([T1_vec T2_vec])*b1_vec))
        b2T=Int.(round.(inv([T1_vec T2_vec])*b2_vec))
        b1fT=Int.(round.(inv([T1_vec T2_vec])*b1f_vec))
        b2fT=Int.(round.(inv([T1_vec T2_vec])*b2f_vec))

        a1_vec=[real(a1),imag(a1)]
        a2_vec=[real(a2),imag(a2)]
        L1_vec=[real(L1),imag(L1)]
        L2_vec=[real(L2),imag(L2)]



        eta1 = wzeta(a1 / (2 * lb), omega=(a1 / (2 * lb), a2 / (2 * lb)))
        Eta1 = wzeta(L1 / (2 * Lb), omega=(L1 / (2 * Lb), L2 / (2 * Lb)))
        eta1f = wzeta(a1f / (2 * lbf), omega=(a1f / (2 * lbf), a2f / (2 * lbf)))
        
        α = eta1 * lb / a1 - conj(a1) / (4 * a1)
        αf = eta1f * lbf / a1f - conj(a1f) / (4 * a1f)
        bigA = Eta1 * Lb / L1 - conj(L1) / (4 * L1)

        
        Tgrid=Vector{Int}[]
      
        cutoffstandard=norm(b1_vec)*grid_cutoff
        cutoff=5*Int(round.(max(cutoffstandard/norm(T1_vec),cutoffstandard/norm(T2_vec))))
        for ja in -cutoff:cutoff, jb in -cutoff:cutoff
            gg=ja*[1,0]+jb*[0,1]
            gg_vec=[T1_vec T2_vec]*gg

            if norm(gg_vec)<cutoffstandard
                push!(Tgrid,gg) 
            end

        end

        Tgrid_dict=Dict{Vector{Int64},Int64}()
        for ja in eachindex(Tgrid)
            Tgrid_dict[Tgrid[ja]]=ja
        end


        gfgrid=Vector{Int}[]
        cutoffstandard=norm(b1_vec)*(grid_cutoff+2)
        cutoff=5*Int(round.(max(cutoffstandard/norm(b1f_vec),cutoffstandard/norm(b2f_vec))))
        for ja in -cutoff:cutoff, jb in -cutoff:cutoff
            gg=ja*b1fT+jb*b2fT
            gg_vec=[T1_vec T2_vec]*gg

            if norm(gg_vec)<cutoffstandard
                push!(gfgrid,gg) 
            end

        end

        gfgrid_dict=Dict{Vector{Int64},Int64}()
        for ja in eachindex(gfgrid)
            gfgrid_dict[gfgrid[ja]]=ja
        end



        kf_set=vec(wrapparallel.([n1*T1+n2*T2 for n1 in 0:N1f-1, n2 in 0:N2f-1],b1f,b2f))
        kf_set_n=Vector{Int}[]
        for ja in eachindex(kf_set)
          nn=Int.(round.(inv([T1_vec T2_vec])*[real(kf_set[ja]), imag(kf_set[ja])]))
          push!(kf_set_n,nn)
        end


        test_Landau_level_sigma(L1,L2,Lb,bigA
                                 ,a1,a2,lb,α,b1,b2,
                                 a1f,a2f,lbf,αf,b1f,b2f)



        
        xi_set=vec([wrapparallel(xi_00+ja*a1+jb*a2,L1,L2) for ja in 0:N1-1, jb in 0:N2-1])






           k_set=vec([wrapparallel(n1*T1+n2*T2,b1,b2) for n1 in 1:N1, n2 in 1:N2])
                
            k_set_n=Vector{Int}[]
            for ja in eachindex(k_set)
            nn=round.(inv([T1_vec T2_vec])*[real(k_set[ja]), imag(k_set[ja])])
            push!(k_set_n,nn)
            end


            ggrid=Vector{Int}[]
            cutoffstandard=norm(b1_vec)*(grid_cutoff+2)
            cutoff=5*Int(round.(max(cutoffstandard/norm(b1_vec),cutoffstandard/norm(b2_vec))))
            for ja in -cutoff:cutoff, jb in -cutoff:cutoff
                gg=ja*b1T+jb*b2T
                gg_vec=[T1_vec T2_vec]*gg

                if norm(gg_vec)<cutoffstandard
                    push!(ggrid,gg) 
                end

            end

            ggrid_dict=Dict{Vector{Int64},Int64}()
            for ja in eachindex(ggrid)
                ggrid_dict[ggrid[ja]]=ja
            end


            qq=-((b1+b2)/2+im*xi_00/lb^2)

            Bloch_states=zeros(ComplexF64,length(Tgrid),N1*N2)

            for k1i in eachindex(k_set)
                k1_n=k_set_n[k1i]
                k1=k_set[k1i]


                    for gi in eachindex(ggrid)
                        g_n=ggrid[gi]
                        gg=g_n[1]*T1+g_n[2]*T2
                        g_n_in_b=inv([b1_vec b2_vec])*[real(gg),imag(gg)]

                        Rk2=-im*lb^2*(-qq)
                        pp=get(Tgrid_dict,k1_n+g_n,0)
                    
                        ct_2=exp(im*π* g_n_in_b[1]*g_n_in_b[2])*exp(-lb^2/4*(conj(gg)*gg+2*(k1-Deltatheta)*conj(gg)))
                        ct_3=exp(-im*c_dot(k1+gg-Deltatheta,Rk2))
                        ct_4=1/get_spinor_norm(NL,[real(k1-Deltatheta+gg),imag(k1-Deltatheta+gg)])

                        if pp≠0

                        Bloch_states[pp,k1i]+=ct_2*ct_3*ct_4
                        
                        end
                    
                    
                    end

            end

            for ja in 1:N1*N2
            Bloch_states[:,ja]=Bloch_states[:,ja]/norm(Bloch_states[:,ja])
            end








        
        params=(a1=a1,
                a2=a2,
                lb=lb,
                b1=b1,
                b2=b2,
                α=α,

                L1=L1,
                L2=L2,
                Lb=Lb,
                T1=T1,
                T2=T2,
                bigA=bigA,

                a1f=a1f,
                a2f=a2f,
                lbf=lbf,
                b1f=b1f,
                b2f=b2f,
                αf=αf,

                N1=N1,
                N2=N2,
                N1f=N1f,
                N2f=N2f,
                Nvac=Nvac,
                Nphi=Nphi,
                Nelectron=Nelectron,

                Deltatheta=Deltatheta,
               
                T1_vec=T1_vec,
                T2_vec=T2_vec,
                b1_vec=b1_vec,
                b2_vec=b2_vec,
                b1f_vec=b1f_vec,
                b2f_vec=b2f_vec,
                b1T=b1T,
                b2T=b2T,
                b1fT=b1fT,
                b2fT=b2fT,

                a1_vec=a1_vec,
                a2_vec=a2_vec,
                L1_vec=L1_vec,
                L2_vec=L2_vec,

                Tgrid=Tgrid,
                Tgrid_dict=Tgrid_dict,
                gfgrid=gfgrid,
                gfgrid_dict=gfgrid_dict,
                grid_cutoff=grid_cutoff,
                ggrid= ggrid,
                ggrid_dict=ggrid_dict,
                Bloch_states= Bloch_states,
                qq=qq,


                kf_set=kf_set,
                kf_set_n=kf_set_n,

                k_set=k_set,
                k_set_n=k_set_n,


                xi_set=xi_set,
                topo_sec=topo_sec,
                NL=NL,
                type=type
               
        )




        possible_config=[randperm(N1*N2)[1:(N1*N2-N1f*N2f)] for _ in 1:sample_num]
        vac_fac_list=Vector{sigma_related}(undef,length(possible_config))
        orbital_Rmatrix_list=Vector{Matrix{ComplexF64}}(undef,length(possible_config))
        orbital_basis_norm_list=Vector{Vector{ComplexF64}}(undef,length(possible_config))
        orbital_basis_orthogonal_list=Vector{Matrix{ComplexF64}}(undef,length(possible_config))
        M_eta_list=Vector{sigma_related}(undef,length(possible_config))
        Mkkmatrix_list=Vector{Matrix{sigma_related}}(undef,length(possible_config))
        N_coeff_list=Vector{Vector{ComplexF64}}(undef,length(possible_config))

        for ja in eachindex(possible_config)
        
            eta_set=xi_set[possible_config[ja]]
            xi_set_partial=xi_set[setdiff(collect(1:N1*N2),possible_config[ja])]

            # We now shift a bit
            shift_set=[[0,0] for _ in 1:Nvac]
         



            for jb in 1:Nvac
                eta_set[jb]+=shift_set[jb][1]*L1+shift_set[jb][2]*L2
            end
            @assert length(eta_set)==Nvac
            @assert length(xi_set_partial)==Nelectron





              vac_fac_list[ja],  orbital_Rmatrix_list[ja], orbital_basis_norm_list[ja],  orbital_basis_orthogonal_list[ja], M_eta_list[ja], Mkkmatrix_list[ja], N_coeff_list[ja]=fixed_configuration_result(eta_set,shift_set,xi_set_partial,params)

        end


       
        spinor_set=Vector{ComplexF64}[]
        for ja in eachindex(Tgrid)
            tvec=Tgrid[ja][1]*T1+Tgrid[ja][2]*T2-Deltatheta
            vv=get_spinor(NL,[real(tvec),imag(tvec)])
            push!(spinor_set,vv)

        end










     


       

        weight_list=zeros(ComplexF64,length(possible_config))
        PN_list=zeros(ComplexF64,length(possible_config))
        
        for ja in eachindex(possible_config)
         weight_list[ja]=tr(orbital_basis_orthogonal_list[ja]'*(Bloch_states*Bloch_states')*orbital_basis_orthogonal_list[ja])/tr(orbital_basis_orthogonal_list[ja]'*orbital_basis_orthogonal_list[ja])
         PN_list[ja]=tr(orbital_basis_orthogonal_list[ja]'*orbital_basis_orthogonal_list[ja])
        end







       



        return params,overall_mag_list,spinor_set,possible_config, vac_fac_list,orbital_Rmatrix_list,
               orbital_basis_norm_list, orbital_basis_orthogonal_list,M_eta_list,Mkkmatrix_list,N_coeff_list,weight_list,PN_list
     


end





