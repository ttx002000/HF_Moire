using LinearAlgebra
function get_ABCA_Ham(k::Vector{Float64},uD::Float64)
 

 
  γ0=3120
  γ1=377
  γ2=-20.6
  γ3=290
  γ4=120
  γ5=25

 
  ff=0.246*√3/2*(-k[1]+im*k[2])
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


  ff=0.246*√3/2*(-k[1]+im*k[2])
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


function get_MLG_Ham(k::Vector{Float64})
    kp=k[1]+im*k[2]
    km=k[1]-im*k[2]
    Ham=zeros(ComplexF64,2,2)

 

    γ=√3/2*0.246*3160
 

    Ham[1,2]+=γ*km


    Ham=Ham+Ham'

    return Ham

end




function main(radius::Float64,Γ::Float64,uD::Float64,whichstack::Int,ϵspace::Vector{Float64},num_points::Int64)
   tic=time()
   kx_grid=collect(range(-radius/2, stop=+radius/2, length=num_points))
   ky_grid=collect(range(-radius/2, stop=+radius/2, length=num_points))
   Area=4*π^2/((ky_grid[2]-ky_grid[1])*(kx_grid[2]-kx_grid[1]))

    Hxmatrix=Matrix{Matrix{ComplexF64}}(undef,length(kx_grid),length(ky_grid))
    Hymatrix=Matrix{Matrix{ComplexF64}}(undef,length(kx_grid),length(ky_grid))
    eigenvalue_record=Matrix{Vector{Float64}}(undef,length(kx_grid),length(ky_grid))
   
    if whichstack==1
      get_Ham = k -> get_ABAB_Ham(k, uD)
    elseif whichstack==2
      get_Ham = k -> get_ABCA_Ham(k, uD)
    elseif whichstack==3
      get_Ham = k -> get_MLG_Ham(k)
    else
      error("Invalid stack type specified")
    end

    dHx=(get_Ham([0.5,0.0])-get_Ham([0.0,0.0]))/0.5
    dHy=(get_Ham([0.0,0.5])-get_Ham([0.0,0.0]))/0.5


    Threads.@threads for ja in eachindex(kx_grid)
        for jb in eachindex(ky_grid)
    
   
  
        k = [kx_grid[ja], ky_grid[jb]]
        H = get_Ham(k)
   
        FFF=eigen(H)

        Hxmatrix[ja,jb]=(FFF.vectors')*dHx*FFF.vectors
        Hymatrix[ja,jb]=(FFF.vectors')*dHy*FFF.vectors
        eigenvalue_record[ja,jb]=real(FFF.values)
        end
   end
   toc=time()
   println("getting Hamiltonian takes time",toc-tic)
   flush(stdout)

   #Gmatrix=[diagm([1/(ϵ-real(eigenvalue_record[ja,jb][ii])+im*Γ) for ii in eachindex(eigenvalue_record[ja,jb]) ]) for ja in eachindex(kx_grid),jb in eachindex(ky_grid)]
   
   #Fz=1/Area*imag(sum([tr(Gmatrix[ja,jb]*Hxmatrix[ja,jb]*Gmatrix[ja,jb]*Hymatrix[ja,jb]*Gmatrix[ja,jb]*Hxmatrix[ja,jb]*Gmatrix[ja,jb]*Hymatrix[ja,jb]) for ja in eachindex(kx_grid),jb in eachindex(ky_grid)]))
  


Fz = zeros(Float64, length(ϵspace))


# Pre-compute constant terms
inv_area = 1.0 / Area
gamma_complex = im * Γ

# Optimized version
  Threads.@threads for jc in eachindex(ϵspace)
    tic=time()
    
      
      ε = ϵspace[jc]
      Gmatrix = Array{Vector{ComplexF64}}(undef, length(kx_grid), length(ky_grid))
      # Compute Green's function matrix more efficiently
      @inbounds for ja in eachindex(kx_grid), jb in eachindex(ky_grid)
          # Vectorized computation of Green's function
          Gmatrix[ja, jb] = @. 1 / (ε - real(eigenvalue_record[ja, jb]) + gamma_complex)
      end
      
      # Compute the trace sum more efficiently
      trace_sum = zero(ComplexF64)
      @inbounds for ja in eachindex(kx_grid), jb in eachindex(ky_grid)
          # Get local references to avoid repeated indexing
          #G = Gmatrix[ja, jb]
          #Hx = Hxmatrix[ja, jb]
          #Hy = Hymatrix[ja, jb]
          
          # Compute (Hx.*G) and (Hy.*G) once
          HxG = Hxmatrix[ja, jb] .* Gmatrix[ja, jb]
          HyG = Hymatrix[ja, jb] .* Gmatrix[ja, jb]
          
          # Compute the product matrix and its trace
          # This is (Hx.*G)*(Hy.*G)*(Hx.*G)*(Hy.*G)
          product_matrix = (HxG * HyG) * (HxG * HyG)
          trace_sum += tr(product_matrix)
      end
      
      Fz[jc] = inv_area * imag(trace_sum)
        toc=time()
          println("Processing energy point: $jc/$(length(ϵspace)),takestime$(toc-tic)")
  end

    
  return Fz
end



function main_bigmemory(radius::Float64,Γ::Float64,uD::Float64,whichstack::Int,ϵspace::Vector{Float64},num_points::Int64)
 
   kx_grid=collect(range(-radius/2, stop=+radius/2, length=num_points))
   ky_grid=collect(range(-radius/2, stop=+radius/2, length=num_points))
   Area=4*π^2/((ky_grid[2]-ky_grid[1])*(kx_grid[2]-kx_grid[1]))
       # Pre-compute constant terms
    inv_area = 1.0 / Area
    gamma_complex = im * Γ
    Fz = zeros(Float64, length(ϵspace))

    Hxmatrix=Matrix{Matrix{ComplexF64}}(undef,length(kx_grid),length(ky_grid))
    Hymatrix=Matrix{Matrix{ComplexF64}}(undef,length(kx_grid),length(ky_grid))
    eigenvalue_record=Matrix{Vector{Float64}}(undef,length(kx_grid),length(ky_grid))
   
    if whichstack==1
      get_Ham = k -> get_ABAB_Ham(k, uD)
    elseif whichstack==2
      get_Ham = k -> get_ABCA_Ham(k, uD)
    elseif whichstack==3
      get_Ham = k -> get_MLG_Ham(k)
    else
      error("Invalid stack type specified")
    end
    
    dHx=(get_Ham([0.5,0.0])-get_Ham([0.0,0.0]))/0.5
    dHy=(get_Ham([0.0,0.5])-get_Ham([0.0,0.0]))/0.5

    for ja in eachindex(kx_grid)
      tic=time()
     
     
      for jb in eachindex(ky_grid)
        
        k = [kx_grid[ja], ky_grid[jb]]
        H = get_Ham(k)
        FFF=eigen(H)

        Hxmatrix=(FFF.vectors')*dHx*FFF.vectors
        Hymatrix=(FFF.vectors')*dHy*FFF.vectors
     
        eigenvalue_record=real(FFF.values)
        trace_sum=zeros(ComplexF64,length(ϵspace))
     
        Threads.@threads for jc in eachindex(ϵspace)
             ε=ϵspace[jc]
             Gmatrix = @. 1 / (ε - real(eigenvalue_record) + gamma_complex)
             HxG = Hxmatrix .* Gmatrix
             HyG = Hymatrix .* Gmatrix
             product_matrix = (HxG * HyG) * (HxG * HyG)
             trace_sum[jc] += tr(product_matrix)

             

        end
        Fz+= inv_area * imag(trace_sum)

      end
      toc=time()
       println("Processing kx: $(ja)/$(num_points),takes time $(tco-tic)")
            flush(stdout)
    end




    
  return Fz
end