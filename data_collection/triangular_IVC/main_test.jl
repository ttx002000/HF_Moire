using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames


include(joinpath(@__DIR__,"../../src/operators_IVC.jl"))

args=parse.(Float64,ARGS)
flux=args[1]*π
V0=args[2]
ϕ=args[3]/180*π
Nq=Int(args[4]);
scale=args[5];
constq=args[6]/Nq^2
ζ=args[7]
trytimes=Int(args[8])
shift_index=[Int(args[9]),Int(args[10])]

#=
    flux=arg[1]
    V0=arg[2]
    ϕ=arg[3]
    Nq=arg[4];
    scale=arg[5];
    constq=arg[6];
    ζ=arg[7];
    trytimes=arg[8]
    shiftindex=arg[9]
=#

if shift_index≠[1,1]
overlapmatrix, wave, initial_DensityMatrix, single_MoirePo, single_Ham,seed_Ham, single_eigenvalue,allowedq, T1, T2, a1m, a2m=triangle_initial_Densitymatrix_withseed(flux,V0,ϕ,scale,Nq,shift_index)
end


if shift_index==[1,1]
    overlapmatrix, wave, initial_DensityMatrix, single_MoirePo, single_Ham,seed_Ham, single_eigenvalue,allowedq, T1, T2, a1m, a2m=triangle_initial_Densitymatrix(flux,V0,ϕ,scale,Nq,shift_index)
end


off_overlap=get_offdiagoverlap(wave,T1,T2,allowedq,flux,scale,Nq,shift_index)
   


NoHFdensity=Densitymap(a1m,a2m,overlapmatrix,wave,initial_DensityMatrix,T1,T2)
NoHF_sxmap=sxmap(a1m,a2m,off_overlap,wave,initial_DensityMatrix,T1,T2)
NoHF_symap=symap(a1m,a2m,off_overlap,wave,initial_DensityMatrix,T1,T2)
 

for ja in 1:Nq^2
   A=randn(ComplexF64,2*length(wave),2*length(wave))
    initial_DensityMatrix[ja]=initial_DensityMatrix[ja]+(A+A')*0.1
end

#if shift_index≠[1,1]
DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,energy_input=iteration_loop_withseed(initial_DensityMatrix,allowedq,T1,T2,Nq,wave,single_Ham,single_MoirePo,seed_Ham,constq,ζ,overlapmatrix)
#end

#if shift_index==[1,1]
 #   DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue,energy_input=iteration_loop(initial_DensityMatrix,allowedq,T1,T2,Nq,wave,single_Ham,single_MoirePo,constq,ζ,overlapmatrix)
#end

HFdensity=Densitymap(a1m,a2m,overlapmatrix,wave,DIIS_input_DensityMatrix[1],T1,T2)
HF_sxmap=sxmap(a1m,a2m,off_overlap,wave,DIIS_input_DensityMatrix[1],T1,T2)
HF_symap=symap(a1m,a2m,off_overlap,wave,DIIS_input_DensityMatrix[1],T1,T2)


chern,Flink,chern_single,Flink_single,trace_condition,trace_condition_single,uniform,uniform_single=triangle_chern(Nq,wave,scale,ϕ,flux,DIIS_input_DensityMatrix[1],constq,ζ,allowedq[shift_index])

energy=calculate_energy(Nq,wave,scale,ϕ,flux,DIIS_input_DensityMatrix[1],constq,ζ,overlapmatrix,allowedq[shift_index])
jldsave(joinpath(@__DIR__, "data_output/$(args[1])flux$(args[2])V0$(args[3])phi$(args[4])Nq$(args[5])scale$(args[6])constq$(args[7])zeta$(args[8])try$(Int(args[9]))$(Int(args[10]))shift.jld2"),chern=chern,Flink=Flink,chern_single=chern_single,Flink_single=Flink_single,HFdensity=HFdensity,NoHFdensity=NoHFdensity,arguments=args,energy=energy,TC=trace_condition,TCS=trace_condition_single,HFeigenvalue=HF_eigenvalue,single_eigenvalue=single_eigenvalue,uniform=uniform,uniform_single=uniform_single,DM=DIIS_input_DensityMatrix[1],energy_input=energy_input,HF_symap=HF_symap,HF_sxmap=HF_sxmap,NoHF_symap=NoHF_symap,NoHF_sxmap=NoHF_sxmap)

