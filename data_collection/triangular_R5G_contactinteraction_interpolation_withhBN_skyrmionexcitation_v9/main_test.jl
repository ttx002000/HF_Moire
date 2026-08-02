using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames

# Basically the same as v2. But we put a1 on x axis, and rather than theta we use am. Also we add in ED functional to compare with the analytic energy
include("../../src/operators_R5G_contactinteraction_interpolation_withhBN_skyrmionexcitation_v9.jl")

args=parse.(Float64,ARGS)
#args=[5,3,0.6,30.0*10^4,5,35.0,1,3.51,1.0,1,1.0,28.9,21.0,-0.29,5.0,0.0,1.0]


NL=Int(args[1])

Nq=Int(args[2]);
moiream=args[3];
constq=args[4]
ϵr=args[5]
uD=args[6]
filling=Int(args[7])
gcutoff=args[8]
λ=args[9]
trytimes=Int(args[10])
enlarge_factor=Int(args[11])
V0_hBN=args[12]
V1_hBN=args[13]
ψ_hBN=args[14]
V2_scalar=args[15]
ϕ=args[16]/180*π
gateD=args[17]
flux1=args[18]
flux2=args[19]
Nparticle=Int(args[20])
filepos=Int(args[21])






overlapmatrix, wave, wave_diff,initial_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,single_eigenvector,allowedq, T1, T2, a1m, a2m, b1,b2,spinor_set=triangle_initial_Densitymatrix(Int(NL),moiream,Nq,gcutoff,uD,λ,enlarge_factor,V0_hBN,V1_hBN,ψ_hBN,V2_scalar,ϕ,flux1,flux2)

Area=Nq^2*√3/2*norm(a1m)^2


DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,eout,HartreeMatrix,FockMatrix=iteration_loop(initial_DensityMatrix,
                                                    allowedq,T1,T2,Nq,wave,single_Ham,single_MoirePo,constq,ϵr,overlapmatrix,filling,Area,gateD)


chern,Flink,chern_single,Flink_single,trace_condition,trace_condition_single,uniform,uniform_single=triangle_chern(Nq,wave,allowedq,HF_eigenvector,single_eigenvector,spinor_set)




 values_record,eig_vec_record,  all_MB_state_can, all_MB_state_integer,Fmatrix,G_dic_record,rho_mat_record=do_ED(allowedq,wave,NL,wave_diff,
               spinor_set,HF_eigenvector,Nq,
               gateD,ϵr,Area,constq,Nparticle)









scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "triangle_R5G_contact_interpolation_withhBN_skyrmionexcitation_v9/data_output$(Int(args[21]))/$(args[1])NL$(args[2])Nq$(args[3])theta$(args[4])constq$(args[5])ϵr$(args[6])uD$(args[7])filling$(args[8])cutoff$(args[9])lambda$(args[10])trytime$(args[11])enlarge$(args[12])V0_hBN$(args[13])V1_hBN$(args[14])ψ_hBN$(args[15])V2_scalar$(args[16])ϕ$(args[17])gD$(args[18])f1$(args[19])f2$(args[20])Npa.jld2")



jldsave(savepath,single_Ham=single_Ham,single_MoirePo=single_MoirePo,
               single_eigenvalue=single_eigenvalue,spinor_set=spinor_set,
                chern=chern,Flink=Flink,chern_single=chern_single,
                Flink_single=Flink_single,arguments=args,
                densitymatrix=DIIS_input_DensityMatrix[1],energy=energy,eout=eout,
                TC=trace_condition,TCS=trace_condition_single,
                HFeigenvalue=HF_eigenvalue,
                uniform=uniform,uniform_single=uniform_single,
                HartreeMatrix=HartreeMatrix,FockMatrix=FockMatrix,HF_eigenvector=HF_eigenvector,
                single_eigenvector=single_eigenvector,allowedq=allowedq,T1=T1,T2=T2,wave=wave,
                a1m=a1m,a2m=a2m,b1=b1,b2=b2,
                 values_record=values_record,eig_vec_record=eig_vec_record, 
                  all_MB_state_can=all_MB_state_can, all_MB_state_integer=all_MB_state_integer,
                  Fmatrix= Fmatrix,
                  G_dic_record=G_dic_record,
                  rho_mat_record=rho_mat_record)


