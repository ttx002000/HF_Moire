using JLD2
include("../../src/operators_onetwo_HF.jl")
using LinearAlgebra

args=parse.(Float64,ARGS)

wAA=args[1]
wAB=args[2]
vF=args[3] #unit meV*nm

lambda_MDT=args[4]
Nb_up=Int(args[5])
Nb_down=Int(args[6])
θ=args[7]/180*π
ϵr=args[8]
Nband=Nb_up+Nb_down
geonum=Int(args[9])
filling=Int(args[10])
trytime=Int(args[11])
#=
wAA=75.0
wAB=110.0
vF=579.2265 #unit meV*nm

lambda_MDT=-0.2
Nb_up=3
Nb_down=3
θ=1.23/180*π
ϵr=8.0
Nband=Nb_up+Nb_down
geonum=Int(1)
filling=3*4+8
=#


eigenvector,eigenvalue,wave,wave_diff,wave_dic,allowedq,allowedq_dic,T1,T2,constq,Minv,g1mT,g2mT,a1m,a2m=single_particle(geonum,θ,wAA,wAB,vF,ϵr,Nband,lambda_MDT,Nb_down,Nb_up)
Npa=length(allowedq)*filling
formfactors=get_formfactors(allowedq,wave,wave_diff,wave_dic,Minv,Nband,eigenvector)
initial_projector, bg_projector, single_Ham=get_initial_proj(allowedq,eigenvalue,Nband)


HF_eigenvalue,HF_eigenvector,energy, DIIS_input_projector,bound=iteration(formfactors,initial_projector,bg_projector,constq,Nband,wave_diff,allowedq,allowedq_dic,T1,T2,single_Ham,Npa)
xgrid,ygrid,HF_density=(eigenvector,DIIS_input_projector[1]-bg_projector,a1m,a2m,wave,T1,T2)

savepath=joinpath(@__DIR__, "data_output/$(args[1])wAA$(args[2])wAB$(args[3])vf$(args[4])MDT$(args[5])Nup$(args[6])Ndown$(args[7])theta$(args[8])er$(args[9])geo$(args[10])fill$(args[11])try.jld2")

jldsave(savepath,HF_eigenvalue=HF_eigenvalue,eigenvalue=eigenvalue,energy=energy,xgrid=xgrid,ygrid=ygrid,HF_density=HF_density,bound=bound)