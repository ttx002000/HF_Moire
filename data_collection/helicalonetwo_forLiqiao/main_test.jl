using JLD2
include("../../src/operators_helicalonetwo_HF_Liqiao.jl")
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
Dfield=args[11]
shift=Int(args[12])
trytime=Int(args[13])
g_cutoff=args[14]
q_cutoff=args[15]
filepos=Int(args[16])

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

eigenvector,eigenvalue,wave,wave_diff,wave_dic,allowedq,allowedq_dic,T1,T2,constq,Minv,g1mT,g2mT,a1m,a2m=single_particle(geonum,θ,wAA,wAB,vF,ϵr,Nband,lambda_MDT,Nb_down,Nb_up,Dfield,shift,g_cutoff,q_cutoff)
Npa=length(allowedq)*filling
#formfactors=get_formfactors(allowedq,wave,wave_diff,wave_dic,Minv,Nband,eigenvector)
#scratch_dir = ENV["SCRATCH"]
#ffpath=joinpath(scratch_dir, "helicalonetwo_Liqiao/data_output$(Int(args[16]))/FF/FF_$(args[1])wAA$(args[2])wAB$(args[3])vf$(args[4])MDT$(args[5])Nup$(args[6])Ndown$(args[7])theta$(args[9])geo$(args[11])Dfield$(args[12])shift$(args[14])gcut$(args[15])qcut.jld2")
#jldsave(ffpath,formfactors=formfactors)


scratch_dir = ENV["SCRATCH"]
ffpath=joinpath(scratch_dir, "helicalonetwo_Liqiao/data_output$(Int(args[16]))/FF/FF_$(args[1])wAA$(args[2])wAB$(args[3])vf$(args[4])MDT$(args[5])Nup$(args[6])Ndown$(args[7])theta$(args[9])geo$(args[11])Dfield$(args[12])shift$(args[14])gcut$(args[15])qcut.jld2")
st=load(ffpath)
formfactors=st["formfactors"]

perturb_Ham=get_bias(g1mT,g2mT,T1,T2,eigenvector,allowedq,a1m,a2m,Nband,shift)
initial_projector, bg_projector, single_Ham=get_initial_proj(allowedq,eigenvalue,Nband,perturb_Ham,Npa)


HF_eigenvalue,HF_eigenvector,energy, DIIS_input_projector,bound,Hartree_matrix,Fock_matrix=iteration(formfactors,initial_projector,bg_projector,constq,Nband,wave_diff,allowedq,allowedq_dic,T1,T2,single_Ham,perturb_Ham,Npa,Minv)

chern_number,chern_num__nonabelian=get_chernnumber(HF_eigenvector,eigenvector,allowedq,allowedq_dic,wave,wave_dic,Nband,geonum,Minv)

savepath=joinpath(scratch_dir, "helicalonetwo_Liqiao/data_output$(Int(args[16]))/$(args[1])wAA$(args[2])wAB$(args[3])vf$(args[4])MDT$(args[5])Nup$(args[6])Ndown$(args[7])theta$(args[8])er$(args[9])geo$(args[10])fill$(args[11])Dfield$(args[12])shift$(args[13])try$(args[14])gcut$(args[15])qcut.jld2")

jldsave(savepath,
        HF_eigenvalue=HF_eigenvalue,eigenvalue=eigenvalue,energy=energy,final_density_matrix=DIIS_input_projector[1],
        bound=bound,chern_number=chern_number,
       chern_num__nonabelian=chern_num__nonabelian,
       Hartree_matrix=Hartree_matrix,Fock_matrix=Fock_matrix,single_Ham=single_Ham,
       T1=T1,T2=T2,g1mT=g1mT,g2mT=g2mT,allowedq=allowedq,Minv=Minv,a1m=a1m,a2m=a2m)

