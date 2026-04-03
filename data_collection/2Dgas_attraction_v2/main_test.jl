using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2


# in v2, I am using unscreened interaction
include("../../src/2Dgas_attraction_v2.jl")

args=parse.(Float64,ARGS)

scale=args[1]
gcutoff=args[2]
constq=args[3]
density=args[4]
lpo=args[5]
attstr=args[6]
temp=args[7]
trytime=Int(args[8])
filepos=Int(args[9])




wave, initial_DensityMatrix,  single_Ham, T1, T2, a1m, a2m, b1,b2,Area=triangle_initial_Densitymatrix(scale,gcutoff)


scratch_dir = ENV["SCRATCH"]

seed_path=joinpath(scratch_dir, "2Dgas_attraction_v2/data_output$(Int(args[9]))/seed")
if only(rand())>0.2
   seed_file_path=pick_random_jld2_path(seed_path)
   if !(seed_file_path==nothing)
      seed_file=load(seed_file_path)
      initial_DensityMatrix=seed_file["densitymatrix"]
      println("Using seed from $seed_file_path")
   end
end


DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,freeenergy,eout,HartreeMatrix,FockMatrix,fermi_level,fermifactor=iteration_loop(initial_DensityMatrix,
                                                                            T1,T2,wave,single_Ham,
                                                                            constq,
                                                                            density,temp,Area,lpo,attstr)





scratch_dir = ENV["SCRATCH"]

savepath=joinpath(scratch_dir, "2Dgas_attraction_v2/data_output$(Int(args[9]))/$(args[1])scale$(args[2])gcut$(args[3])constq$(args[4])density$(args[5])lpo$(args[6])attstr$(args[7])temp$(args[8])trytime.jld2")
#savepath="test.jld2"


jldsave(savepath,single_Ham=single_Ham,
                arguments=args,
                densitymatrix=DIIS_input_DensityMatrix[1],energy=energy,eout=eout,
                HFeigenvalue=HF_eigenvalue,
                HartreeMatrix=HartreeMatrix,FockMatrix=FockMatrix,HF_eigenvector=HF_eigenvector,
                 T1=T1,T2=T2,wave=wave,
                a1m=a1m,a2m=a2m,b1=b1,b2=b2,freeenergy=freeenergy,fermi_level=fermi_level,fermifactor=fermifactor,Area=Area)


