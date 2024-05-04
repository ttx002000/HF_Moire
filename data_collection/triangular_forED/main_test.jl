using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames


include("../../src/operators_forED.jl")

args=parse.(Float64,ARGS)
flux=args[1]*π
V0=args[2]
ϕ=args[3]/180*π
Nx=Int(args[4]);
Ny=Int(args[5]);
scale=args[6];
constq=args[7]/(Nx*Ny)
trytimes=Int(args[8])


overlapmatrix, wave, initial_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,allowedq, T1, T2, a1m, a2m=triangle_initial_Densitymatrix(flux,V0,ϕ,scale,Nx,Ny)
NoHFdensity=Densitymap(a1m,a2m,overlapmatrix,wave,initial_DensityMatrix)
 
for ja in 1:Nx*Ny
    A=randn(ComplexF64,length(wave),length(wave))
    initial_DensityMatrix[ja]=initial_DensityMatrix[ja]+(A+A')*0.1
end

DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue=iteration_loop(initial_DensityMatrix,allowedq,T1,T2,Nx,Ny,wave,single_Ham,single_MoirePo,constq,overlapmatrix)
HFdensity=Densitymap(a1m,a2m,overlapmatrix,wave,DIIS_input_DensityMatrix[1])

chern,Flink,chern_single,Flink_single,trace_condition,trace_condition_single,uniform,uniform_single=triangle_chern(Nx,Ny,wave,scale,ϕ,flux,DIIS_input_DensityMatrix[1],constq)
energy,ED_vec,ED_val=calculate_energy(Nx,Ny,wave,scale,ϕ,flux,DIIS_input_DensityMatrix[1],constq,overlapmatrix)
jldsave(joinpath(@__DIR__, "data_output/$(args[4])Nx$(args[5])Ny$(args[1])flux$(args[2])V0$(args[3])phi$(args[6])scale$(args[7])constq_$(args[8])try.jld2"),chern=chern,Flink=Flink,chern_single=chern_single,Flink_single=Flink_single,HFdensity=HFdensity,NoHFdensity=NoHFdensity,arguments=args,energy=energy,TC=trace_condition,TCS=trace_condition_single,HFeigenvalue=HF_eigenvalue,single_eigenvalue=single_eigenvalue,uniform=uniform,uniform_single=uniform_single,ED_vec=ED_vec,ED_val=ED_val)

