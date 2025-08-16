
using LinearAlgebra

function get_f(k::Vector{Float64})
  delta1=1/√3*0.246*[1,0]
  delta2=1/√3*0.246*[-1/2,√3/2]
  delta3=1/√3*0.246*[-1/2,-√3/2]

  return exp(im*dot(k,delta1))+exp(im*dot(k,delta2))+exp(im*dot(k,delta3))
end




function find_FL(quasi_particle_energy::Vector{Float64},target_density::Float64,val_s::Float64,val_e::Float64,temp::Float64,Area::Float64,bg_particle_density::Float64)
  try_FL=(val_s+val_e)/2

  
  if target_density==0.0
    stan=10^(-11)
  else
    stan=abs(10^(-10)*target_density)
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





function get_ABCA_Ham(k::Vector{Float64},uD::Float64)
 

 
  γ0=3120
  γ1=377
  γ2=-20.6
  γ3=290
  γ4=120
  γ5=25

 
  ff=get_f(k)
  D1=[1.5*uD γ0*ff;γ0*ff' 1.5*uD]
  D2=conj.([0.5*uD γ0*ff;γ0*ff' 0.5*uD]) # This is the one from literature (no conjugate, but I think there should be a conjugate), but is this correct?

  D3=conj.([-0.5*uD γ0*ff;γ0*ff' -0.5*uD])
  D4=[-1.5*uD γ0*ff;γ0*ff' -1.5*uD]

  H12=[γ1 -γ4*conj(ff);-γ4*conj(ff) γ3*ff]
  H21=H12'
  H32=[-γ4*ff γ1;γ3*ff' -γ4*ff] 
  H23=H32'

  H31=[0 0;0 γ2/2]
  H13=H31'


  H14=[0 0;0 0]
  H41=H14'

  H24=[γ2/2 0;0 0]
  H34=[γ3*ff -γ4*ff'; -γ4*ff' γ1]
  H42=H24'
  H43=H34'

  


  Ham=[D1 H12 H13 H14;H21 D2 H23 H24; H31 H32 D3 H34; H41 H42 H43 D4]


   return Ham
end




function get_ABAB_Ham(k::Vector{Float64},uD::Float64)
 
  Ham=zeros(ComplexF64,2*4,2*4)

 
  γ0=3120
  γ1=377
  γ2=-20.6
  γ3=290
  γ4=120
  γ5=25
  δ=36.6-20.6-25.0
  #δ=36.6

  ff=get_f(k)
   Ham[1,2]+=γ0*ff
   Ham[1,3]+=γ1
   Ham[1,4]+=-γ4*conj(ff)
   Ham[1,5]+=γ5/2
   
   Ham[2,3]+=-γ4*conj(ff)
   Ham[2,4]+=γ3*ff
   Ham[2,6]+=γ2/2
   
   Ham[3,4]+=γ0*conj(ff)
   Ham[3,5]+=γ1
   Ham[3,6]+=-γ4*ff  #This one is in the paper they cite
    #Ham[3,6]+=-γ4*conj(ff)  #This is in his own appendix

   Ham[4,5]+=-γ4*ff
   Ham[4,6]+=γ3*conj(ff)
   
   Ham[5,6]+=γ0*ff


   Ham[7,8]+=γ0*conj(ff)
   Ham[5,7]+=γ1
   Ham[5,8]+=-γ4*conj(ff)
   Ham[6,7]+=-γ4*conj(ff)
   Ham[6,8]+=γ3*ff
   Ham[4,8]+=γ2/2
   Ham[3,7]+=γ5/2

   Ham+=Ham'

   Ham[1,1]+=δ+γ5+1.5*uD
   Ham[2,2]+=γ2+1.5*uD
   Ham[3,3]+=δ+γ5+0.5*uD
   Ham[4,4]+=γ2+0.5*uD
   Ham[5,5]+=δ+γ5-0.5*uD
   Ham[6,6]+=γ2-0.5*uD
   Ham[7,7]+=δ+γ5-1.5*uD
   Ham[8,8]+=γ2-1.5*uD



   return Ham
end



function get_reference_CNP(
           temp::Float64,kx_grid::Vector{Int},
           ky_grid::Vector{Int},Nq::Int64,NL::Int)
            values_plot_ABC_reference=zeros(Float64,2*NL,length(kx_grid),length(kx_grid))

            values_plot_ABA_reference=zeros(Float64,2*NL,length(kx_grid),length(kx_grid))

            Threads.@threads for ja in eachindex(kx_grid)
              for jb in eachindex(ky_grid)
            
                k=kx_grid[ja]/(Nq)*G1+ky_grid[jb]/Nq*G2
            
                single_Ham=get_ABCA_Ham(k,0.0)
                
                values_plot_ABC_reference[:,ja,jb]=real(eigen(single_Ham).values)
              end
            end



            Threads.@threads for ja in eachindex(kx_grid)
              for jb in eachindex(ky_grid)
            
                k=kx_grid[ja]/(Nq)*G1+ky_grid[jb]/Nq*G2

                single_Ham=get_ABAB_Ham(k,0.0)
                values_plot_ABA_reference[:,ja,jb]=real(eigen(single_Ham).values)

            end
            end



            CNP_ABA_reference,_=find_FL(vec(values_plot_ABA_reference),0.0,sort(vec(values_plot_ABA_reference))[1],sort(vec(values_plot_ABA_reference))[end],temp,Area,Nq^2/Area*NL)
            CNP_ABC_reference,_=find_FL(vec(values_plot_ABC_reference),0.0,sort(vec(values_plot_ABC_reference))[1],sort(vec(values_plot_ABC_reference))[end],temp,Area,Nq^2/Area*NL)




            aba_record=copy(CNP_ABA_reference)
            abc_record=copy(CNP_ABC_reference)

           

    return aba_record, abc_record
end 


function big_func(uD::Float64,aba_record::Float64,abc_record::Float64,
         workfunction::Float64,temp::Float64,kx_grid::Vector{Int},
         ky_grid::Vector{Int},Nq::Int64,NL::Int64,density_list::Vector{Float64})

            
            values_plot_ABC=zeros(Float64,2*NL,length(kx_grid),length(ky_grid))

            values_plot_ABA=zeros(Float64,2*NL,length(kx_grid),length(ky_grid))

            Threads.@threads for ja in eachindex(kx_grid)
              for jb in eachindex(ky_grid)
            
            
                
                k=kx_grid[ja]/(Nq)*G1+ky_grid[jb]/Nq*G2
                single_Ham=get_ABCA_Ham(k,uD)
                
                values_plot_ABC[:,ja,jb]=real(eigen(single_Ham).values)
            
              end
            end



            Threads.@threads for ja in eachindex(kx_grid)
              for jb in eachindex(ky_grid)
            
            
                
                
                k=kx_grid[ja]/(Nq)*G1+ky_grid[jb]/Nq*G2
                single_Ham=get_ABAB_Ham(k,uD)
                values_plot_ABA[:,ja,jb]=real(eigen(single_Ham).values)
            

              end
            end


            v_ABA=sort(vec(values_plot_ABA[:,:,:]))
            v_ABC=sort(vec(values_plot_ABC[:,:,:]))



            CNP_ABA=aba_record
            CNP_ABC=abc_record







            values_plot_ABC=values_plot_ABC.-CNP_ABC
            values_plot_ABA=values_plot_ABA.-CNP_ABA

            v_ABA=sort(vec(values_plot_ABA[:,:,:]))
            v_ABC=sort(vec(values_plot_ABC[:,:,:]))


            CNP_ABA=0.0
            CNP_ABC=0.0

            

            bg_ABA=Nq^2/Area*NL
            bg_ABC=Nq^2/Area*NL



            energy_diff=zeros(Float64,length(density_list))

            Threads.@threads for ja in eachindex(density_list)
                println(ja)
            mu1,_=find_FL(v_ABA,density_list[ja],v_ABA[1],v_ABA[end],temp,Area,bg_ABA)
            mu2,_=find_FL(v_ABC,density_list[ja],v_ABC[1],v_ABC[end],temp,Area,bg_ABC)
              
              jj1=0.0
              jj2=0.0
               for jb in eachindex(v_ABA)
                jj1+=1/(1+exp((v_ABA[jb]-mu1)/temp))*v_ABA[jb]
                jj2+=1/(1+exp((v_ABC[jb]-mu2)/temp))*v_ABC[jb]
                #jj1=sum(([1/(1+exp((v_ABA[ja]-mu1)/temp)) for ja in eachindex(v_ABA)]/Area).*v_ABA)
                #jj2=sum(([1/(1+exp((v_ABC[ja]-mu2)/temp)) for ja in eachindex(v_ABC)]/Area).*v_ABC)
               end
               jj1=jj1/Area
               jj2=jj2/Area
            energy_diff[ja]=(jj2-jj1-density_list[ja]*(workfunction))*2
            end
            println(uD)
            flush(stdout)
            return energy_diff
end



