using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2


include("../../src/construct_twothirds_Lauhglin_ED.jl")

args=parse.(Float64,ARGS)

flux1=args[1]
flux2=args[2]
q1=args[3]
q2=args[4]
Nx=Int(args[5])
Ny=Int(args[6])
Nparticle=Int(args[7])
num_vecs=Int(args[8])
file_pos=Int(args[9])





PH_eigvector,PH_state_can,allowedq,T1,T2, 
values_record,eig_vec_record, state_can_record, state_integer_record,
Laughlin_state_vector,Laughlin_state_sector,Laughlin_state_can, Laughlin_state_integer,
   a1m,a2m,g1,g2=do_ED(Nx,Ny,Nparticle,
              flux1,flux2,q1,q2,
              num_vecs)





     



scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "constrcut_twothirds_Lauhglin_ED/data_output$(Int(args[9]))/$(args[1])f1$(args[2])f2$(args[3])q1$(args[4])q2$(args[5])N1$(args[6])N2$(args[7])Npa$(args[8])Nvec.jld2")




jldsave(savepath,PH_eigvector=PH_eigvector,PH_state_can=PH_state_can,allowedq=allowedq,
values_record=values_record,eig_vec_record=eig_vec_record, state_can_record= state_can_record, 
state_integer_record=state_integer_record,
Laughlin_state_vector=Laughlin_state_vector,
Laughlin_state_sector=Laughlin_state_sector,
Laughlin_state_can=Laughlin_state_can,
 Laughlin_state_integer= Laughlin_state_integer,
  a1m=a1m,a2m=a2m,g1=g1,g2=g2,T1=T1,T2=T2)


