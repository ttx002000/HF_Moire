using LinearAlgebra,Statistics 
using Random

function get_Ham(k::Vector{Float64},potential_profile::Vector{Float64},valley::Int64,NL::Int,Ham_ver::Int)
    if Ham_ver==1
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
            Ham[:,layer,:,layer]+=1/2*[potential_profile[layer] t0*ff;t0*conj(ff) potential_profile[layer]]
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
    elseif Ham_ver==2
            Ham=zeros(ComplexF64,2*NL,2*NL)
            
            t0=3100
            t1=380
            t2=-21
            t3=290
            t4=141
                ff=√3/2*0.246*(valley*k[1]-im*k[2])
            for layer in 1:NL-1
                Ham[2*layer-1:2*layer,2*layer+1:2*layer+2]=[t4*ff t3*conj(ff);t1 t4*ff]
            end

            if NL>2
            for layer in 1:NL-2
                Ham[2*layer-1:2*layer,2*layer+3:2*layer+4]=[0.0 t2/2;0.0 0.0]
            end
            end

            Ham=Ham+Ham'

            for layer in 1:NL
                Ham[2*layer-1:2*layer,2*layer-1:2*layer]=[potential_profile[layer] -t0*ff;-t0*conj(ff) potential_profile[layer]]
            end

    

    elseif Ham_ver==3
            Ham=zeros(ComplexF64,2*NL,2*NL)
            
            t0=2600
            t1=356.1
            t2=-8.3
            t3=293
            t4=144
                ff=√3/2*0.246*(valley*k[1]+im*k[2])
            for layer in 1:NL-1
                Ham[2*layer-1:2*layer,2*layer+1:2*layer+2]=[-t4*conj(ff) -t3*ff;t1 -t4*conj(ff)]
            end

            if NL>2
            for layer in 1:NL-2
                Ham[2*layer-1:2*layer,2*layer+3:2*layer+4]=[0.0 t2;0.0 0.0]
            end
            end

            Ham=Ham+Ham'

            for layer in 1:NL
                Ham[2*layer-1:2*layer,2*layer-1:2*layer]=[potential_profile[layer] t0*conj(ff); t0*(ff) potential_profile[layer]]
            end

    end
        
  
 return Ham
end







function get_potential_profile(density_profile::Vector{Float64},top_gate::Float64,bottom_gate::Float64,tg_dis::Float64,bg_dis::Float64,ϵr::Float64)
   tot_density_profile=zeros(Float64,length(density_profile)+2)
   tot_density_profile[1]=bottom_gate
   tot_density_profile[end]=top_gate
   tot_density_profile[2:end-1]=density_profile

   
   tot_pos=zeros(Float64,length(tot_density_profile))
   tot_pos[2:end-1]=[0.335*ja for ja in 1:length(density_profile)]
   tot_pos[1]=tot_pos[2]-bg_dis
   tot_pos[end]=tot_pos[end-1]+tg_dis
   

   field_profile=zeros(Float64,length(tot_density_profile)-1)
   for ja in eachindex(field_profile)
     field_profile[ja]=(sum(tot_density_profile[ja+1:end])-sum(tot_density_profile[1:ja]))*9047.5636/ϵr # In units of mV/nm, pointing up
   end
   
   potential_profile=zeros(Float64,length(tot_density_profile))
   for ja in 2:length(potential_profile)
     potential_profile[ja]=potential_profile[ja-1]+field_profile[ja-1]*(tot_pos[ja]-tot_pos[ja-1])
   end
      #println( potential_profile)

    potential_energy=1/2*sum(potential_profile.*tot_density_profile)

   return potential_profile,potential_energy

end


function find_FL(quasi_particle_energy::Vector{Float64},target_density::Float64,val_s::Float64,val_e::Float64,temp::Float64,Area::Float64,bg_particle_density::Float64)
  try_FL=(val_s+val_e)/2

 
  if target_density==0.0
    stan=10^(-9)
  else
    stan=abs(10^(-8)*target_density)
  end
  fermifactor=[1/(exp((quasi_particle_energy[ja]-try_FL)/temp)+1) for ja in eachindex(quasi_particle_energy)]
 
  fl=sum(fermifactor)/Area-bg_particle_density



 if abs(fl-target_density)<stan
  
    return try_FL,fl
  elseif fl-target_density>=stan
 
    return find_FL(quasi_particle_energy,target_density,val_s, try_FL,temp,Area,bg_particle_density)
  elseif fl-target_density<=-stan
 
    return find_FL(quasi_particle_energy,target_density,try_FL,val_e,temp,Area,bg_particle_density)
   end
 

end


function safe_inverse(A)
  try
      return inv(A)  # Attempt to compute inverse
  catch e
      if isa(e, SingularException)
          println("Matrix is singular, doing randomstart again.")
          return pinv(A, 0.1)  # Use pseudoinverse as an alternative
      else
          rethrow(e)  # If another error occurs, propagate it
      end
  end
end



function implement_DIIS(DIIS_delta::Vector{Vector{Float64}},DIIS_potential::Vector{Vector{Float64}},DIIS_size::Int)


      Bmatrix=zeros(Float64,DIIS_size+1,DIIS_size+1)
      for ja in 1:DIIS_size
       Bmatrix[ja,DIIS_size+1]=1
       Bmatrix[DIIS_size+1,ja]=1
      end
  
      for ja in 1:DIIS_size,jb in 1:DIIS_size
          
             Bmatrix[ja,jb]+=real(DIIS_delta[ja]'*DIIS_delta[jb])
       
      end

      inB=safe_inverse(Bmatrix)
      if inB≠0
         inv_vec=zeros(Float64,DIIS_size+1)
         inv_vec[end]=1
         coeff=inB*inv_vec
         dmk=coeff[1]*(DIIS_potential[1]+DIIS_delta[1])
          for ja in 2:DIIS_size
            dmk+=coeff[ja]*(DIIS_potential[ja]+DIIS_delta[ja])
          end
         return dmk
      else
        return 0
      end
end


function get_density_profile(kx_list::Vector{Float64},ky_list::Vector{Float64},num_layers::Int,num_sub::Int,num_spin::Int,num_valley::Int,
                             Area::Float64,target_density_list::Matrix{Float64},temp::Float64,current_potential_profile::Vector{Float64},Ham_ver::Int)
     
        valley_set=[1,-1]
        eigenvector=zeros(ComplexF64,length(kx_list),length(ky_list),num_valley,num_layers*num_sub,num_layers*num_sub)
        eigenvalue=zeros(Float64,length(kx_list),length(ky_list),num_valley,num_layers*num_sub)

        Hamiltonian=zeros(ComplexF64,length(kx_list),length(ky_list),num_valley,num_layers*num_sub,num_layers*num_sub)


     for vi in 1:num_valley
        for ja in eachindex(kx_list)
            for jb in eachindex(ky_list)
            
                        Ham=get_Ham([kx_list[ja],ky_list[jb]],current_potential_profile,valley_set[vi],num_layers,Ham_ver)
                    
                        Hamiltonian[ja,jb,vi,:,:]=Ham
                        FFF=eigen(Ham)
                        eigenvector[ja,jb,vi,:,:]=FFF.vectors
                        eigenvalue[ja,jb,vi,:]=real(FFF.values)
        
            end
        end
     end
    
       fermi_energy_list=zeros(Float64,num_spin,num_valley)
        for vi in 1:num_valley, si in 1:num_spin
      
         sorted_e=sort(vec(eigenvalue[:,:,vi,:]))
         fermi_energy_list[si,vi],_=find_FL(sorted_e,target_density_list[si,vi],sorted_e[1],sorted_e[end],temp,Area,num_layers*length(kx_list)*length(ky_list)/Area)
        end

        density_layer_resolved=zeros(Float64,num_spin,num_valley,num_layers)

        for ja in eachindex(kx_list)
            for jb in eachindex(ky_list)
                      for vi in 1:num_valley
                        evecs=@view eigenvector[ja,jb,vi,:,:]
                        evecs=reshape(evecs,num_sub,num_layers,num_sub*num_layers)
                        evals=@view eigenvalue[ja,jb,vi,:]
                        for si in 1:num_spin
                         for bi in 1:2*num_layers
                        
                            density_layer_resolved[si,vi,:]+=vec(sum(abs2,evecs[:,:,bi]; dims=1))*1/(1+exp((evals[bi]-fermi_energy_list[si,vi])/temp))
                            
                         end
                        end
                    end

                
        
            end
        end
        density_bg=length(kx_list)*length(ky_list)*num_layers*num_spin*num_valley/(Area*num_layers)*ones(Float64,num_layers)
        density_profile=dropdims(sum(density_layer_resolved, dims=(1,2)), dims=(1,2))/Area.-density_bg

        return density_profile,fermi_energy_list,density_layer_resolved
end



function get_kinetic_energy(kx_list::Vector{Float64},ky_list::Vector{Float64},num_layers::Int,num_sub::Int,num_spin::Int,num_valley::Int,
                            Area::Float64,target_density_list::Matrix{Float64},temp::Float64,current_potential_profile::Vector{Float64},Ham_ver::Int64)
# printing the kinetic energy

        valley_set=[1,-1]
        eigenvector=zeros(ComplexF64,length(kx_list),length(ky_list),num_valley,num_layers*num_sub,num_layers*num_sub)
        eigenvalue=zeros(Float64,length(kx_list),length(ky_list),num_valley,num_layers*num_sub)

        no_diag_Hamiltonian=zeros(ComplexF64,length(kx_list),length(ky_list),num_valley,num_layers*num_sub,num_layers*num_sub)


     for vi in 1:num_valley
        for ja in eachindex(kx_list)
            for jb in eachindex(ky_list)
            
                        Ham=get_Ham([kx_list[ja],ky_list[jb]],current_potential_profile,valley_set[vi],num_layers,Ham_ver)
                        nodiag_Ham=get_Ham([kx_list[ja],ky_list[jb]],zeros(Float64,length(current_potential_profile)),valley_set[vi],num_layers,Ham_ver)
                        no_diag_Hamiltonian[ja,jb,vi,:,:]=nodiag_Ham
                        FFF=eigen(Ham)
                        eigenvector[ja,jb,vi,:,:]=FFF.vectors
                        eigenvalue[ja,jb,vi,:]=real(FFF.values)
        
            end
        end
    end
    
        
        fermi_energy_list=zeros(Float64,num_spin,num_valley)
        for vi in 1:num_valley, si in 1:num_spin
      
         sorted_e=sort(vec(eigenvalue[:,:,vi,:]))
         fermi_energy_list[si,vi],_=find_FL(sorted_e,target_density_list[si,vi],sorted_e[1],sorted_e[end],temp,Area,num_layers*length(kx_list)*length(ky_list)/Area)
        end

        kE_resolved=zeros(Float64,num_spin,num_valley)
    for vi in 1:num_valley
        for ja in eachindex(kx_list)
            for jb in eachindex(ky_list)
            
                        evecs=@view eigenvector[ja,jb,vi,:,:]
                      
                        evals=@view eigenvalue[ja,jb,vi,:]
                        hh=@view no_diag_Hamiltonian[ja,jb,vi,:,:]

                    
                        for si in 1:num_spin, bi in 1:2*num_layers
                            
                            kE_resolved[si,vi]+=real((evecs[:,bi]'*hh*evecs[:,bi]))*1/(1+exp((evals[bi]-fermi_energy_list[si,vi])/temp))
                            
                        end
            
        
            end
        end
    end
    

        return sum(kE_resolved)/Area,kE_resolved/Area
end

function get_Ham_dx(k::Vector{Float64},potential_profile::Vector{Float64},valley::Int64,NL::Int,Ham_ver::Int)
  step=1e-5
 return (get_Ham(k.+[step,0.0],potential_profile,valley,NL,Ham_ver)-get_Ham(k.-[step,0.0],potential_profile,valley,NL,Ham_ver))/(2*step)
end

function get_Ham_dy(k::Vector{Float64},potential_profile::Vector{Float64},valley::Int64,NL::Int,Ham_ver::Int)
  step=1e-5
 return (get_Ham(k.+[0.0,step],potential_profile,valley,NL,Ham_ver)-get_Ham(k.-[0.0,step],potential_profile,valley,NL,Ham_ver))/(2*step)
end


function get_DOS(num_kpoints_forOBM::Int,radius::Float64,num_layers::Int,potential_profile::Vector{Float64},
               target_density_list::Matrix{Float64},temp::Float64,Ham_ver::Int)

    num_spin=2
    num_valley=2
    num_sub=2
    valley_set=[1,-1]
    kxrange=collect(LinRange(-radius,radius,num_kpoints_forOBM))
    kyrange=collect(LinRange(-radius,radius,num_kpoints_forOBM))
    Area=4*π^2/(kxrange[2]-kxrange[1])/(kyrange[2]-kyrange[1])


    eigenvalue=zeros(Float64,length(kxrange),length(kyrange),num_valley,num_sub*num_layers)
    for vi in 1:num_valley
    for ja in eachindex(kxrange)
          for jb in eachindex(kyrange)
            HH=get_Ham([kxrange[ja],kyrange[jb]],potential_profile,valley_set[vi],num_layers,Ham_ver)
            eigenvalue[ja,jb,vi,:]=real.(eigen(HH).values)
          end
    end
   end
       
    fermi_energy_list=zeros(Float64,num_spin,num_valley)
        for vi in 1:num_valley, si in 1:num_spin
      
         sorted_e=sort(vec(eigenvalue[:,:,vi,:]))
         fermi_energy_list[si,vi],_=find_FL(sorted_e,target_density_list[si,vi],sorted_e[1],sorted_e[end],temp,Area,num_layers*length(kxrange)*length(kyrange)/Area)
        end

    nt=Threads.maxthreadid()
    println(fermi_energy_list)
    
    DOS_current=zeros(Float64,num_spin,num_valley)
    DOS_old=zeros(Float64,num_spin,num_valley)
    difference=zeros(Float64,num_spin,num_valley)
    energy_cut=0.1

    each_sc_count=10^6

    for si in 1:num_spin, vi in 1:num_valley
        bin_count=zeros(Int64,nt)
        for big_count in 1:3
                Threads.@threads for sc in 1:each_sc_count
                
                    kx=(rand()-0.5)/0.5*radius
                    ky=(rand()-0.5)/0.5*radius
                    HH=get_Ham([kx,ky],potential_profile,valley_set[vi],num_layers,Ham_ver)
                    vals=eigen(HH).values
                    for jb in 1:2*num_layers
                        if abs(vals[jb]-fermi_energy_list[si,vi])<energy_cut
                            
                            bin_count[Threads.threadid()]+=1
                        end
                    end


                end
                
                Area_DOS=4*π^2/radius^2*each_sc_count*big_count
                DOS_current[si,vi]=1/Area_DOS*sum(bin_count)/(2*energy_cut)
                difference[si,vi]=DOS_current[si,vi]-DOS_old[si,vi]
                DOS_old[si,vi]=DOS_current[si,vi]
                        
                println(difference[si,vi],"difference",sum(bin_count),"binacount",each_sc_count*big_count,"samplecount")
                flush(stdout)
        end
    end

   






  return DOS_current,difference
end














function get_OBM(num_kpoints_forOBM::Int,radius::Float64,num_layers::Int,potential_profile::Vector{Float64},
               target_density_list::Matrix{Float64},temp::Float64,Ham_ver::Int)

    num_spin=2
    num_valley=2
    num_sub=2
    valley_set=[1,-1]
    kxrange=collect(LinRange(-radius,radius,num_kpoints_forOBM))
    kyrange=collect(LinRange(-radius,radius,num_kpoints_forOBM))
    Area=4*π^2/(kxrange[2]-kxrange[1])/(kyrange[2]-kyrange[1])


    eigenvalue=zeros(Float64,length(kxrange),length(kyrange),num_valley,num_sub*num_layers)
    eigenvector=zeros(ComplexF64,length(kxrange),length(kyrange),num_valley,num_sub*num_layers,num_sub*num_layers)
    for vi in 1:num_valley
    for ja in eachindex(kxrange)
          for jb in eachindex(kyrange)
            HH=get_Ham([kxrange[ja],kyrange[jb]],potential_profile,valley_set[vi],num_layers,Ham_ver)
            FFF=eigen(HH)
            eigenvalue[ja,jb,vi,:]=real(FFF.values)
            eigenvector[ja,jb,vi,:,:]=FFF.vectors
          end
    end
    end

   
    fermi_energy_list=zeros(Float64,num_spin,num_valley)
        for vi in 1:num_valley, si in 1:num_spin
      
         sorted_e=sort(vec(eigenvalue[:,:,vi,:]))
         fermi_energy_list[si,vi],_=find_FL(sorted_e,target_density_list[si,vi],sorted_e[1],sorted_e[end],temp,Area,num_layers*length(kxrange)*length(kyrange)/Area)
        end


    M1=zeros(ComplexF64,length(kxrange),length(kyrange),num_valley,num_layers)
    M2=zeros(ComplexF64,length(kxrange),length(kyrange),num_valley,num_layers)
 
    for vi in 1:num_valley
      for ja in eachindex(kxrange)
            for jb in eachindex(kyrange)
                k=[kxrange[ja],kyrange[jb]]
           
              HH_dx=get_Ham_dx(k,potential_profile,valley_set[vi],num_layers,Ham_ver)
              HH_dy=get_Ham_dy(k,potential_profile,valley_set[vi],num_layers,Ham_ver)
              evals=@view eigenvalue[ja,jb,vi,:]
              evecs=@view eigenvector[ja,jb,vi,:,:]
            for ll in num_layers+1:2*num_layers
                for n in 1:2*num_layers
                  if n≠ll
                    vn=evecs[:,n]
                    vtarget=evecs[:,ll]
                    velx=vtarget'*HH_dx*vn
                    vely=vn'*HH_dy*vtarget
                    M1[ja,jb,vi,ll-num_layers]+=-imag(velx*vely)/(evals[ll]-evals[n])
                    M2[ja,jb,vi,ll-num_layers]+=-imag(velx*vely)/(evals[ll]-evals[n])^2
                  end
                end
            end
          

            end
      end
    end

     total_M_resolved=zeros(Float64,num_spin,num_valley)
     for si in 1:num_spin, vi in 1:num_valley
      for ja in eachindex(kxrange)
        for jb in eachindex(kyrange)
            for ll in num_layers+1:2*num_layers
              total_M_resolved[si,vi]+=(M1[ja,jb,vi,ll-num_layers]+2*(fermi_energy_list[si,vi]-eigenvalue[ja,jb,vi,ll])*M2[ja,jb,vi,ll-num_layers])*(1/(1+exp((eigenvalue[ja,jb,vi,ll]-fermi_energy_list[si,vi])/temp)))
           end
       end
      end
    end

    conversion_factor=2*9.109383*10^(-31)*(10^(-9))^2*(10^(-3)* 1.60217663 * 10^(-19))/(1.0545718*10^(-34))^2        
    
     orbital_magnetization=zeros(Float64,num_spin,num_valley)
     for si in 1:num_spin, vi in 1:num_valley
        if target_density_list[si,vi]≠0
        orbital_magnetization[si,vi]=total_M_resolved[si,vi]/(Area*target_density_list[si,vi])*conversion_factor
        end
     end



  return orbital_magnetization,total_M_resolved
end




function iteration_loop(num_kpoints::Int,radius::Float64,uD::Float64,temp::Float64,
                  num_layers::Int,ϵr::Float64,target_density_list::Matrix{Float64},tg_dis::Float64,
                  bg_dis::Float64,Ham_ver::Int64)
        
   
        kx_list=collect(LinRange(-radius,radius,num_kpoints))
        ky_list=collect(LinRange(-radius,radius,num_kpoints))
        Area=4π^2/(kx_list[2]-kx_list[1])/(ky_list[2]-ky_list[1])

        num_spin=2
        num_valley=2
        num_sub=2

        eout=1.0

        DIIS_size=8


        top_gate=-sum(target_density_list)/2+uD*ϵr*5.52635/0.335*10^(-5)
        bottom_gate=-sum(target_density_list)/2-uD*ϵr*5.52635/0.335*10^(-5)


    
        
       
       
        
        itcount=0

        DIIS_density=Vector{Vector{Float64}}(undef,DIIS_size)
        DIIS_delta=Vector{Vector{Float64}}(undef,DIIS_size)
        current_density=randn(num_layers)
        current_density.-=mean(current_density)
        current_density .+= sum(target_density_list)/num_layers
        potential_profile=zeros(Float64,num_layers)
        bad_count=0
        potential_energy=0.0
        fermi_energy_list=zeros(Float64,num_spin,num_valley,num_layers)
        density_layer_resolved=zeros(Float64,num_spin,num_valley,num_layers)

    while eout>10^(-9) || bad_count<DIIS_size
 
        if eout>10^(-8)
            bad_count=0
        end

        if eout<10^(-8)
            bad_count+=1
        end

       if itcount>DIIS_size*2
            println("using DIIS")
            current_density=implement_DIIS(DIIS_delta,DIIS_density,DIIS_size)      
    

        end

        
        DIIS_density[mod(itcount,DIIS_size)+1]=copy(current_density)
        ppp,potential_energy=get_potential_profile(current_density,top_gate,bottom_gate,
                                            tg_dis,bg_dis,ϵr)
        potential_profile=ppp[2:end-1]
        updated_density,fermi_energy_list,density_layer_resolved=get_density_profile(kx_list,ky_list,num_layers,num_sub,num_spin,num_valley,
                                Area,target_density_list,temp,potential_profile,Ham_ver)

        residual_vec=updated_density-current_density
        DIIS_delta[mod(itcount,DIIS_size)+1]=residual_vec

        
        

        current_density=0.5*updated_density+0.5* current_density

        

        





        
        eout=norm(residual_vec)
        println("diff", eout)
        flush(stdout)
        itcount+=1
    end
    
    kinetic_energy,kinetic_energy_resolved=get_kinetic_energy(kx_list,ky_list,num_layers,num_sub,num_spin,num_valley,
                               Area,target_density_list,temp,potential_profile,Ham_ver)                           
   


  return DIIS_density,density_layer_resolved,eout,potential_profile,kinetic_energy,kinetic_energy_resolved,potential_energy,fermi_energy_list

end



