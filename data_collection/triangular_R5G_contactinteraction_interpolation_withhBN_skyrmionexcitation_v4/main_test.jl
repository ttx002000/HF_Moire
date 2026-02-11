using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2


# Basically the same as v1, but this one specifically deals with Nq=1 case, and also spinless
include("../../src/operators_R5G_contactinteraction_interpolation_withhBN_skyrmionexcitation_v4.jl")

args=parse.(Float64,ARGS)



NL=args[1]
θ=args[2]/180*pi;
constq=args[3]
ϵr=args[4]
uD=args[5]
filling=Int(args[6])
gcutoff=args[7]
λ=args[8]
trytimes=Int(args[9])
enlarge_factor=Int(args[10])
V0_hBN=args[11]
V1_hBN=args[12]
ψ_hBN=args[13]
V2_scalar=args[14]
ϕ=args[15]/180*π
filepos=Int(args[16])






overlapmatrix, wave, initial_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,single_eigenvector, T1, T2, a1m, a2m, b1,b2,spinor_set=triangle_initial_Densitymatrix(Int(NL),θ,gcutoff,uD,λ,enlarge_factor,V0_hBN,V1_hBN,ψ_hBN,V2_scalar,ϕ)

Area=√3/2*norm(a1m)^2

if constq<5*10^4
  DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,eout,HartreeMatrix,FockMatrix=iteration_loop(initial_DensityMatrix,
                                                                            T1,T2,wave,single_Ham,
                                                                            single_MoirePo,constq,ϵr,
                                                                 overlapmatrix,filling,Area)
else
    println("contact interaction too large, using strategy 2")
      DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,eout,HartreeMatrix,FockMatrix=iteration_loop(initial_DensityMatrix,
                                                                            T1,T2,wave,single_Ham,
                                                                            3*single_MoirePo,5*10^4,ϵr,
                                                                 overlapmatrix,filling,Area)
        last_input=copy(DIIS_input_DensityMatrix[1])
    DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy,eout,HartreeMatrix,FockMatrix=iteration_loop(last_input,
                                                                            T1,T2,wave,single_Ham,
                                                                            single_MoirePo,constq,ϵr,
                                                                 overlapmatrix,filling,Area)
end





scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "triangle_R5G_contact_interpolation_withhBN_skyrmionexcitation_v4/data_output$(Int(args[16]))/$(args[1])NL$(args[2])theta$(args[3])constq$(args[4])ϵr$(args[5])uD$(args[6])filling$(args[7])cutoff$(args[8])lambda$(args[9])trytime$(args[10])enlarge$(args[11])V0_hBN$(args[12])V1_hBN$(args[13])ψ_hBN$(args[14])V2_scalar$(args[15])ϕ.jld2")
#savepath="test.jld2"


jldsave(savepath,single_Ham=single_Ham,single_MoirePo=single_MoirePo,
               single_eigenvalue=single_eigenvalue,spinor_set=spinor_set,
                arguments=args,
                densitymatrix=DIIS_input_DensityMatrix[1],energy=energy,eout=eout,
                HFeigenvalue=HF_eigenvalue,
                HartreeMatrix=HartreeMatrix,FockMatrix=FockMatrix,HF_eigenvector=HF_eigenvector,
                single_eigenvector=single_eigenvector,T1=T1,T2=T2,wave=wave,
                a1m=a1m,a2m=a2m,b1=b1,b2=b2)


