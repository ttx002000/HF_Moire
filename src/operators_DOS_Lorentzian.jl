using LinearAlgebra,Plots,StatsBase

function get_f(k::Vector{Float64})
  delta1=1/√3*0.246*[0,1]
  delta2=1/√3*0.246*[√3/2,-1/2]
  delta3=1/√3*0.246*[-√3/2,-1/2]

  return exp(im*dot(k,delta1))+exp(im*dot(k,delta2))+exp(im*dot(k,delta3))
end

function big_func(type)


    
    function get_RNG_Ham(k::Vector{Float64},uD::Float64,NL::Int)
 
            Ham=zeros(ComplexF64,2*NL,2*NL)

            t0=-3120
        t1=377
        t2=-20.6
        t3=290
        t4=-120
        ff=get_f(k)
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
            Ham[2*layer-1:2*layer,2*layer-1:2*layer]=[uD*(layer-(NL+1)/2) -t0*ff;-t0*conj(ff) uD*(layer-(NL+1)/2)]
        end

        return Ham
    end



        function get_ABA_Ham(k::Vector{Float64},uD::Float64)
        
                Ham=zeros(ComplexF64,2*3,2*3)
                
                t0=-3120
                t1=377
                t2=-20.6
                t3=290
                t4=-120
                t5=25
                δ=36.6-20.6-25.0
                #δ=36.6 
                ff=get_f(k)

                
                layer=1
                Ham[2*layer-1:2*layer,2*layer+1:2*layer+2]+=[t4*ff t3*conj(ff);t1 t4*ff]
                layer=2
                Ham[2*layer-1:2*layer,2*layer+1:2*layer+2]+=[t4*conj(ff) t1; t3*ff t4*ff]
                


                layer=1
                    Ham[2*layer-1:2*layer,2*layer+3:2*layer+4]=[t2/2 0.0;0.0 t5/2]


                Ham=Ham+Ham'

                for layer in 1:3
                    Ham[2*layer-1:2*layer,2*layer-1:2*layer]=[uD*(layer-(3+1)/2) -t0*ff;-t0*conj(ff) uD*(layer-(3+1)/2)]
                end


                Ham[1,1]+=t2
                Ham[2,2]+=t5+δ
                
                Ham[3,3]+=t5
                Ham[4,4]+=t2

                Ham[5,5]+=t2
                Ham[6,6]+=t5+δ

                return Ham
        end

        if type==1
            g = (x, y) -> get_RNG_Ham(x, y, 3)
            return g
        elseif type==2
            return get_ABA_Ham
        end


end

function record_values(Nq::Int,type::Int,Γ::Float64,uD::Float64,Area::Float64)

    
    values_record=[Vector{Float64}[] for _ in 1:Threads.nthreads()]
    index_set=vec([[ja,jb] for ja in 1:Nq, jb in 1:Nq])
    Ham_func=big_func(type)



  Threads.@threads for ja in eachindex(index_set)
    k=index_set[ja][1]/Nq*G1+index_set[ja][1]/Nq*G2
    h1=Ham_func(k,uD)
    s1=real(eigen(h1).values)
    push!(values_record[Threads.threadid()],s1)
  end

  values_record=sort(reduce(vcat,reduce(vcat,values_record)))
  Elist=collect(range(values_record[1]-20*Γ, stop=values_record[end]+20*Γ, length=Int(round((values_record[end]-values_record[1]+40*Γ)*0.5))))
  DOS=zeros(Float64,length(Elist))

  
  Threads.@threads for ja in eachindex(Elist)
    println(ja)
    DOS[ja]=1/(π*Area)*sum([Γ/(Γ^2+(values_record[jb]-Elist[ja])^2) for jb in eachindex(values_record)])
  end


  return Elist,DOS

end
