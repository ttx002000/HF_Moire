

function get_f(k::Vector{Float64})
   delta1=1/√3*0.246*[0,1]
   delta2=1/√3*0.246*[√3/2,-1/2]
   delta3=1/√3*0.246*[-√3/2,-1/2]
 
   return exp(im*dot(k,delta1))+exp(im*dot(k,delta2))+exp(im*dot(k,delta3))
end
 
 
function Hamiltonian(k::Vector{Float64},valley::Int64,uD::Float64,)
   NL=5
   Ham=zeros(ComplexF64,2*NL,2*NL)
   Kac=4π/(3*0.246)*[1,0]*valley
   t0=3100
   t1=380
   t2=-21
   t3=290
   t4=141
   for layer in 1:NL-1
      Ham[2*layer-1:2*layer,2*layer+1:2*layer+2]=[t4*get_f(k+Kac) t3*conj(get_f(k+Kac));t1 t4*get_f(k+Kac)]
   end
 
   for layer in 1:NL-2
       Ham[2*layer-1:2*layer,2*layer+3:2*layer+4]=[0.0 t2/2;0.0 0.0]
   end
 
   Ham=Ham+Ham'
 
   for layer in 1:NL
       Ham[2*layer-1:2*layer,2*layer-1:2*layer]=[uD*(layer-(NL+1)/2) -t0*get_f(k+Kac);-t0*conj(get_f(k+Kac)) uD*(layer-(NL+1)/2)]
   end
  return Ham
 end



function Coulomb(k1::Vector{Float64},k2::Vector{Float64},ϵr::Float64)::Float64
    d=1500
    if k1==k2
      return 9047.5636/(ϵr)*d
    else
      return 9047.5636/(ϵr*norm(k1-k2))*(1-exp(-norm(k1-k2)*d))
    end
  
end

function find_chemical_potential(carrier_density::Vector{Float64},uD::Float64,bandindex::Int)

    radius=2.0
    num_points=400
    kx_grid=collect(range(-radius/2, stop=+radius/2, length=num_points))
    ky_grid=collect(range(-radius/2, stop=+radius/2, length=num_points))
    Area=4*π^2/((kx_grid[2]-kx_grid[1])*(ky_grid[2]-ky_grid[1]))
    eig_set=zeros(Float64,length(kx_grid),length(ky_grid))
    
    for ja in eachindex(kx_grid),jb in eachindex(ky_grid)
       Ham=Hamiltonian([kx_grid[ja],ky_grid[jb]],1,uD)
       FFF=eigen(Ham)
       eig_set[ja,jb]=-real(FFF.values[bandindex])
       
    end

    density_list=[ja/Area for ja in 1:num_points^2]

    
    chemical_potential=zeros(Float64,length(carrier_density))
    for ja in eachindex(carrier_density)
     p1=searchsortedfirst(density_list,carrier_density[ja])
     chemical_potential[ja]=sort(vec(eig_set))[p1]
    end
    
    
  

   return chemical_potential

end


function find_Qlength(uD::Float64,bandindex::Int)

    radius=2.0
    num_points=400
    kx_grid=collect(range(-radius/2, stop=+radius/2, length=num_points))
    ky_grid=collect(range(-radius/2, stop=+radius/2, length=num_points))
    
    eig_set=zeros(Float64,length(kx_grid),length(ky_grid))
    k_set=Matrix{Vector{Float64}}(undef,length(kx_grid),length(ky_grid))
    
    for ja in eachindex(kx_grid),jb in eachindex(ky_grid)
       Ham=Hamiltonian([kx_grid[ja],ky_grid[jb]],1,uD)
       FFF=eigen(Ham)
       eig_set[ja,jb]=-real(FFF.values[bandindex])
       k_set[ja,jb]=[kx_grid[ja],ky_grid[jb]]
    end



    
    s1=sortperm(vec(eig_set))

    Qlength=norm(vec(k_set)[s1[1]])

   return   Qlength 

end


function sample_states(uD::Float64,num_u_grid::Int,num_theta_grid::Int,density_start::Float64,density_stop::Float64,cutoff::Float64)

  num_density_grid=15
  density_point=collect(range(density_start, stop=density_stop, length=num_density_grid))*0.01 #The 0.01 is for unit conversion
  bandindex=5
  chemical_potential=find_chemical_potential(density_point,uD,bandindex)
  


 Qlength=find_Qlength(uD,bandindex)
 umin=-(Qlength)^(1/3)
 umax=(1.5*Qlength)^(1/3)


 ugrid=collect(range(umin, stop=umax, length=num_u_grid))
 θgrid=collect(range(0, stop=2π, length=num_theta_grid+1))
 θgrid=θgrid[1:end-1]


 eig_set=Float64[]
 eig_vec_set=Vector{ComplexF64}[]
 k_set=Vector{Float64}[]
 theta_set=Float64[]
 u_set=Float64[]

 for ja in eachindex(ugrid),jb in eachindex(θgrid)
   kvec=(Qlength+ugrid[ja]^3)*[cos(θgrid[jb]),sin(θgrid[jb])]

   Ham=Hamiltonian(kvec,1,uD)
   FFF=eigen(Ham)
   push!(eig_set,-real(FFF.values[bandindex]))
   push!(eig_vec_set,FFF.vectors[:,bandindex])
   push!(k_set,kvec)
   push!(theta_set,θgrid[jb])
   push!(u_set,ugrid[ja]) 
 end

 eig_set_final=Float64[]
 eig_vec_set_final=Vector{ComplexF64}[]
 k_set_final=Vector{Float64}[]
 theta_set_final=Float64[]
 u_set_final=Float64[]
 
 


    for ja in eachindex(eig_set)
 
        if eig_set[ja]<chemical_potential[end]+cutoff
          push!(eig_set_final,eig_set[ja])
          push!(eig_vec_set_final,eig_vec_set[ja])
          push!(k_set_final,k_set[ja])
          push!(theta_set_final,theta_set[ja])
          push!(u_set_final,u_set[ja])
       end
    end

 


 
 s1=sortperm(eig_set_final)
 eig_set_final=eig_set_final[s1]
 eig_vec_set_final=eig_vec_set_final[s1]
 k_set_final=k_set_final[s1]
 theta_set_final=theta_set_final[s1]
 u_set_final=u_set_final[s1]
 
 println("I keep this many states",length(eig_set_final))



  return  eig_set_final,eig_vec_set_final,k_set_final,theta_set_final,u_set_final,umin,umax,Qlength,density_point,chemical_potential,ugrid,θgrid
end


function get_formfactors(eig_set_final::Vector{Float64},eig_vec_set_final::Vector{Vector{ComplexF64}},
                         k_set_final::Vector{Vector{Float64}},theta_set_final::Vector{Float64},
                         u_set_final::Vector{Float64},ϵr::Float64,ugrid::Vector{Float64},θgrid::Vector{Float64})


    formfactors=zeros(Float64,length(eig_set_final),length(eig_set_final))
   
    measureone=zeros(Float64,length(eig_set_final))
    
    tic=time()
    
    Threads.@threads for ja in eachindex(eig_set_final)
      for jb in eachindex(eig_set_final)
        formfactors[ja,jb]=abs(dot(eig_vec_set_final[ja],eig_vec_set_final[jb]))^2*Coulomb(k_set_final[ja],k_set_final[jb],ϵr)
      end
    end
    toc=time()
    println(toc-tic,"form factors takes this amount of time")
    flush(stdout)
    
    Threads.@threads for ja in eachindex(eig_set_final)
      
        f1=3*u_set_final[ja]^2*(u_set_final[ja]^3+Qlength)
        measureone[ja]=(1/(2*π)^2*(θgrid[2]-θgrid[1])*(ugrid[2]-ugrid[1]))*f1
    end
    
   
   
   

    return formfactors,measureone


end



function find_FL(quasi_particle_energy::Vector{Float64},measureone::Vector{Float64},target_density::Float64,val_s::Float64,val_e::Float64,temp::Float64)
   try_FL=(val_s+val_e)/2
 
   stan=10^(-5)*target_density
 
   fermifactor=[1/(exp((quasi_particle_energy[ja]-try_FL)/temp)+1) for ja in eachindex(measureone)]
  
   fl=sum(measureone.*fermifactor)
 
 
 
  if abs(fl-target_density)<stan
     return try_FL,fermifactor,fl
   elseif fl-target_density>=stan
     return find_FL(quasi_particle_energy,measureone,target_density,val_s, try_FL,temp)
   elseif fl-target_density<=-stan
     return find_FL(quasi_particle_energy,measureone,target_density,try_FL,val_e,temp)
    end
  
 
end



function calculate_energy(measureone::Vector{Float64},target_density::Float64,
                          formfactors::Matrix{Float64},temp::Float64,m1::Matrix{Float64},eig_set_final::Vector{Float64})
 
 

 quasi_particle_energy_old=randn(length(eig_set_final))
 quasi_particle_energy_new=randn(length(eig_set_final))
 itnum=0
 energy_new=10.0
 energy_old=0.0
 fermilevel=0.0
 fermifactor=[]
 renormalized_density=0.0
 

 while abs(energy_new-energy_old)>10^(-5)
  itnum+=1
  energy_old=energy_new
  fermilevel,fermifactor,renormalized_density=find_FL(quasi_particle_energy_old,
                                                      measureone,target_density,
                                                      sort(quasi_particle_energy_old)[1],sort(quasi_particle_energy_old)[end],temp)
  
 
  quasi_particle_energy_new=[eig_set_final[ja]-sum(formfactors[ja,:].*measureone.*fermifactor) for ja in eachindex(eig_set_final)]
 
  quasi_particle_energy_old=quasi_particle_energy_new
  energy_new=sum(eig_set_final.*measureone.*fermifactor)/renormalized_density
  energy_new-=transpose(fermifactor)*m1*fermifactor/(2*renormalized_density)

 
  println("iterations",itnum)
  println("energy change",energy_new-energy_old)
  println("new energy",energy_new)
  println("renormalized_density",renormalized_density)
  println("fermilevel",fermilevel)
  end

  return fermifactor,fermilevel,renormalized_density,energy_new,quasi_particle_energy_new
end

function do_iterations(measureone::Vector{Float64},density_point::Vector{Float64},formfactors::Matrix{Float64},temp::Float64,eig_set_final::Vector{Float64})
   m1=(formfactors.*(measureone*transpose(measureone)))
   fermifactor_final=Vector{Vector{Float64}}(undef,length(density_point))
   fermilevel_final=Vector{Float64}(undef,length(density_point))
   renormalized_density_final=Vector{Float64}(undef,length(density_point))
   energy_final=Vector{Float64}(undef,length(density_point))
   quasi_particle_energy_final=Vector{Vector{Float64}}(undef,length(density_point))

   Threads.@threads for ja in eachindex(density_point)
    fermifactor_final[ja],fermilevel_final[ja],renormalized_density_final[ja],energy_final[ja],quasi_particle_energy_final[ja]=calculate_energy(measureone,density_point[ja],formfactors,temp,m1,eig_set_final)
   end


   return  fermifactor_final,fermilevel_final,renormalized_density_final,energy_final,quasi_particle_energy_final
end