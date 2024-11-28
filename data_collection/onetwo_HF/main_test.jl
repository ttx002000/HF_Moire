using JLD2
include("../../src/operators_onetwo_HF.jl")
using LinearAlgebra

args=parse.(Float64,ARGS)

#args=[75.0,110.0,500.0,-0.2,3,3,1.23,8,1,20,10,1,0.0]

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
trytime=Int(args[13])
shift=Int(args[12])
Dfield=args[11]

#=
wAA=75.0
wAB=110.0
vF=579.2265 #unit meV*nm

lambda_MDT=-0.2
Nb_up=1
Nb_down=1
θ=1.23/180*π
ϵr=8.0
Nband=Nb_up+Nb_down
geonum=Int(2)
filling=1
Dfield=20.0
=#

eigenvector,eigenvalue,wave,wave_diff,wave_dic,allowedq,allowedq_dic,T1,T2,constq,Minv,g1mT,g2mT,a1m,a2m,g_cutoff,q_cutoff=single_particle(geonum,θ,wAA,wAB,vF,ϵr,Nband,lambda_MDT,Nb_down,Nb_up,Dfield,shift)
Npa=length(allowedq)*filling
#formfactors=get_formfactors(allowedq,wave,wave_diff,wave_dic,Minv,Nband,eigenvector)
#scratch_dir = ENV["SCRATCH"]
#savepath=joinpath(scratch_dir, "onetwo/FF_$(args[1])wAA$(args[2])wAB$(args[3])vf$(args[4])MDT$(args[5])Nup$(args[6])Ndown$(args[7])theta$(args[9])geo$(args[11])Dfield$(args[12])shift$(g_cutoff)gcut$(q_cutoff)qcut.jld2")
#jldsave(savepath,formfactors=formfactors)


scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "onetwo/FF_$(args[1])wAA$(args[2])wAB$(args[3])vf$(args[4])MDT$(args[5])Nup$(args[6])Ndown$(args[7])theta$(args[9])geo$(args[11])Dfield$(args[12])shift$(g_cutoff)gcut$(q_cutoff)qcut.jld2")
st=load(savepath)
formfactors=st["formfactors"]

perturb_Ham=get_bias(g1mT,g2mT,T1,T2,eigenvector,allowedq,a1m,a2m,Nband,shift)
initial_projector, bg_projector, single_Ham=get_initial_proj(allowedq,eigenvalue,Nband,perturb_Ham,Npa)


HF_eigenvalue,HF_eigenvector,energy, DIIS_input_projector,bound,Hartree_matrix,Fock_matrix=iteration(formfactors,initial_projector,bg_projector,constq,Nband,wave_diff,allowedq,allowedq_dic,T1,T2,single_Ham,perturb_Ham,Npa,Minv)
xgrid,ygrid,HF_density=plot_Chargedensity(eigenvector,DIIS_input_projector[1],a1m,a2m,wave,T1,T2,allowedq)


layer_pol,sub_pol=get_polarization(HF_eigenvector,eigenvector,allowedq,wave,Nband)
sub_exp_HF,chernsub_HF,sub_eig_HF,sub_eig_single=get_chernsub(HF_eigenvector,eigenvector,allowedq,wave,Nband)

chern_number,chern_num__nonabelian=get_chernnumber(HF_eigenvector,eigenvector,allowedq,allowedq_dic,wave,wave_dic,Nband,geonum,Minv)

#savepath=joinpath(@__DIR__, "data_output/shift$(Int(args[12]))/$(args[1])wAA$(args[2])wAB$(args[3])vf$(args[4])MDT$(args[5])Nup$(args[6])Ndown$(args[7])theta$(args[8])er$(args[9])geo$(args[10])fill$(args[11])Dfield$(args[12])shift$(args[13])try.jld2")
savepath=joinpath(scratch_dir, "onetwo/data_output/shift$(Int(args[12]))/$(args[1])wAA$(args[2])wAB$(args[3])vf$(args[4])MDT$(args[5])Nup$(args[6])Ndown$(args[7])theta$(args[8])er$(args[9])geo$(args[10])fill$(args[11])Dfield$(args[12])shift$(args[13])try.jld2")

jldsave(savepath,
        HF_eigenvalue=HF_eigenvalue,eigenvalue=eigenvalue,energy=energy,
        xgrid=xgrid,ygrid=ygrid,HF_density=HF_density,
        bound=bound,layer_pol=layer_pol,sub_pol=sub_pol,
       chern_number=chern_number,g_cutoff=g_cutoff,
       q_cutoff=q_cutoff,
       sub_exp_HF=sub_exp_HF,chernsub_HF=chernsub_HF,sub_eig_HF=sub_eig_HF, sub_eig_single= sub_eig_single,
       chern_num__nonabelian=chern_num__nonabelian,
       Hartree_matrix=Hartree_matrix,Fock_matrix=Fock_matrix,single_Ham=single_Ham)

