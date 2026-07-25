using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2


include("../../src/construct_twothirds_Lauhglin_ED_step2.jl")

args=parse.(Float64,ARGS)


#=

     flux1=args[1]
     flux2=args[2]
     q1=args[3]
     q2=args[4]
     N1=Int(args[5])
   N2=Int(args[6])
   Npa_LL=Int(args[7])
   Nvec=Int(args[8])
   file_pos=Int(args[9])
   grid_cutoff=args[10]
   NL=Int(args[11])
   moiream=args[12]
=#

rho_matrix_all,
        final_G_q_all,
        deno_all,
        reconstruct_matrix,
        overlap_manybody,
        coeff_orthogonal,
        orbital_basis_orthogonal,
        params=do_Gq(args)



     



scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "constrcut_twothirds_Lauhglin_ED/data_output$(Int(args[9]))/onetwobody/$(args[1])f1$(args[2])f2$(args[3])q1$(args[4])q2$(args[5])N1$(args[6])N2$(args[7])Npa$(args[8])Nvec$(args[10])gcut$(args[11])NL$(args[12])am.jld2")




jldsave(savepath,rho_matrix_all=rho_matrix_all,
        final_G_q_all=final_G_q_all,
        deno_all=deno_all,
        reconstruct_matrix=reconstruct_matrix,
        overlap_manybody=overlap_manybody,
        coeff_orthogonal=coeff_orthogonal,
        orbital_basis_orthogonal=orbital_basis_orthogonal,
        params=params)


