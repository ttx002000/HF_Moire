using JLD2
include("../../src/concentratedBC.jl")
using LinearAlgebra

args=parse.(Float64,ARGS)
#args=[60.0,3.0,8.0,0.0,1.0,1.0,0.5]
Nq=Int(args[1])
spin=args[2]
scale=args[3]
ϕ=args[4]/180*π
V0=args[5]
mass=args[6]
M=args[7]


directgap,BW,indirectgap,chern_eigenvector,wave,T1,T2=get_eigenvector(Nq,spin,scale,ϕ,V0,mass,M)
chern,uniform,trace_condition=get_chern(Nq,wave,spin,M,T1,T2)

#directgap,BW,indirectgap,chern_eigenvector,wave,T1,T2=get_eigenvector_square(Nq,spin,scale,ϕ,V0,mass,M)
#chern,uniform,trace_condition=get_chern_square(Nq,wave,spin,M,T1,T2)


jldsave(joinpath(@__DIR__, "data_output/square_$(args[1])Nq$(args[2])spin$(args[3])scale$(args[4])phi$(args[5])V0$(args[6])mass$(args[7])M.jld2"),chern=chern,directgap=directgap,BW=BW,indirectgap=indirectgap,uniform=uniform,trace_condition=trace_condition)
