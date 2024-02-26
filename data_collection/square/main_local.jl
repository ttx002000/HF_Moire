include("../../src/operators.jl")
flux=2*π
V0=0.5
ϕ=0.0
Nq=3;
scale=1.0;
constq=0.1/Nq^2



overlapmatrix, wave, initial_DensityMatrix, single_MoirePo, single_Ham, single_eigenvalue,allowedq, T1, T2, a1m, a2m=square_initial_Densitymatrix(flux,V0,ϕ,scale,Nq)
NoHFdensity=Densitymap(a1m,a2m,overlapmatrix,wave,initial_DensityMatrix)
 
#=for ja in 1:Nq^2
    A=randn(ComplexF64,length(wave),length(wave))
    initial_DensityMatrix[ja]=initial_DensityMatrix[ja]+(A+A')*0.001
end=#

DIIS_input_DensityMatrix,DIIS_input_DeltaMatrix,HF_eigenvalue=iteration_loop(initial_DensityMatrix,allowedq,T1,T2,Nq,wave,single_Ham,single_MoirePo,constq,overlapmatrix)
HFdensity=Densitymap(a1m,a2m,overlapmatrix,wave,DIIS_input_DensityMatrix[1])

chern,Flink,chern_single,Flink_single=square_chern(Nq,wave,scale,ϕ,flux,DIIS_input_DensityMatrix[1],constq)
println(Flink)
