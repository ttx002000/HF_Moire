using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames


include("../../src/operators_finiteS_forED.jl")

args=parse.(Float64,ARGS)
spin=args[1]
V0=args[2]
ϕ=args[3]/180*π
Geonum=Int(args[4]);
Nx,Ny,l1,l2=Geometry(Geonum)
scale=args[5];
constq=args[6]/(Nx*Ny)
trytimes=Int(args[7])
vf=8*sqrt(spin*sqrt(3)/(4*π))

overlapmatrix, wave, initial_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,allowedq, T1, T2, a1m, a2m, b1T,b2T=triangle_initial_Densitymatrix(spin,vf,V0,ϕ,scale,Geonum)
NoHFdensity=Densitymap(a1m,a2m,overlapmatrix,wave,initial_DensityMatrix)
 
for ja in eachindex(allowedq)
    A=randn(ComplexF64,length(wave),length(wave))
    initial_DensityMatrix[ja]=initial_DensityMatrix[ja]+(A+A')*0.1
end

DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,HF_eigenvector,energy=iteration_loop(initial_DensityMatrix,allowedq,T1,T2,wave,single_Ham,single_MoirePo,constq,overlapmatrix)
HFdensity=Densitymap(a1m,a2m,overlapmatrix,wave,DIIS_input_DensityMatrix[1])


jldsave(joinpath(@__DIR__, "data_output/$(args[1])spin$(args[2])V0$(args[3])phi$(args[4])Geom$(args[5])scale$(args[6])constq$(args[7])try_idealvF.jld2"),HFdensity=HFdensity,NoHFdensity=NoHFdensity,arguments=args,energy=energy,HFeigenvalue=HF_eigenvalue,single_eigenvalue=single_eigenvalue,allowedq=allowedq,T1=T1,T2=T2,b1T=b1T,b2T=b2T,a1m=a1m,a2m=a2m)


