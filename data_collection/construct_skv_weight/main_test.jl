using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2


include("../../src/construct_skv_weight.jl")

args=parse.(Float64,ARGS)

flux1=args[1]
flux2=args[2]
NL=Int(args[3])
moiream=args[4]
N1=Int(args[5])
N2=Int(args[6])
N1f=Int(args[7])
N2f=Int(args[8])
xi00_re=args[9]
xi00_im=args[10]
grid_cutoff=args[11]
topo_sec=Int(args[12])
type=Int(args[13])
sample_num=Int(args[14])
trytime=Int(args[15])
file_pos=Int(args[16])






params,overall_mag_list,spinor_set,possible_config, orbital_Rmatrix_list,
               orbital_basis_norm_list, orbital_basis_orthogonal_list,M_eta_list,Mkkmatrix_list,N_coeff_list,weight_list,PN_list=big_func(flux1,flux2,NL,moiream,
                  N1,N2,N1f,N2f,xi00_re+im*xi00_im,
                  grid_cutoff,topo_sec,type,sample_num)





     



scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "construct_skv_weight/data_output$(Int(args[16]))/$(args[1])f1$(args[2])f2$(args[3])NL$(args[4])am$(args[5])N1$(args[6])N2$(args[7])N1f$(args[8])N2f$(args[9])xir$(args[10])xii$(args[11])grid$(args[12])topo$(args[13])type$(args[14])snum$(args[15])try.jld2")




jldsave(savepath,params=params,overall_mag_list=overall_mag_list,spinor_set=spinor_set,
            possible_config=possible_config,orbital_Rmatrix_list=orbital_Rmatrix_list,
               orbital_basis_norm_list=orbital_basis_norm_list, orbital_basis_orthogonal_list=orbital_basis_orthogonal_list,M_eta_list=M_eta_list,
               Mkkmatrix_list=Mkkmatrix_list,N_coeff_list=N_coeff_list,weight_list=weight_list,PN_list=PN_list)


