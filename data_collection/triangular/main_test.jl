using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames


include("../../src/operators.jl")

args=parse.(Int64,ARGS)
flux=args[2]
epsilon=args[3]

Nx=5;
Ny=6;
Nparticle=10;
onebodyreal = CSV.read(joinpath(@__DIR__, "data_input/onebodyreal$(flux)_e$(epsilon).csv"), DataFrame, header=false)
onebodyimag = CSV.read(joinpath(@__DIR__, "data_input/onebodyimag$(flux)_e$(epsilon).csv"), DataFrame, header=false)
onebodymatrix=Matrix(onebodyreal) .+ (im*Matrix(onebodyimag))

twobodyreal = CSV.read(joinpath(@__DIR__, "data_input/twobodyreal$(flux)_e$(epsilon).csv"), DataFrame, header=false)
twobodyimag = CSV.read(joinpath(@__DIR__, "data_input/twobodyimag$(flux)_e$(epsilon).csv"), DataFrame, header=false)
twobodymatrix=Matrix(twobodyreal) .+ (im*Matrix(twobodyimag))


(reduced_Vcol,reduced_Vcoor,eigenvalue_single,allowedq)=convertdata(onebodymatrix,twobodymatrix,Nx,Ny,Nparticle)
(MB_state_can, MB_state_integer)=Construct_MBstate(Nx,Ny,Nparticle,allowedq);
state_can=MB_state_can[args[1]]
state_integer=MB_state_integer[args[1]]
MB_state_can=nothing
MB_state_integer=nothing

values=Construct_Manybodymatrix(reduced_Vcol,reduced_Vcoor,state_can,state_integer,eigenvalue_single)


jldsave(joinpath(@__DIR__, "data_output/eigenvalue$(args[1])sector_ep$(epsilon)flux$(flux).jld2"),values=values)


