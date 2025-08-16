using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames


include("../../src/operators_workfunction.jl")

args=parse.(Float64,ARGS)
#args=[50.0,600,1.5/11.6045,4,-22.0,-0.1,0.1,81]
uD=args[1]
Nq=Int(args[2])
temp=args[3]
NL=Int(args[4])
wf=args[5]
den_st=args[6]
den_end=args[7]
den_num=Int(args[8])
file_pos=Int(args[9])



ac=0.246

G1=4π/(√3*ac)*[1,0]
G2=4π/(√3*ac)*[1/2,√3/2]
Area=4π^2/(norm(G1/Nq)^2*√3/2)
num_atom=Area/(√3/2*0.246^2)*2


kx_grid=collect(1:1:Nq)
ky_grid=collect(1:1:Nq)

density_list=collect(range(den_st, stop=den_end, length=den_num))

aba_record, abc_record=get_reference_CNP(
                            temp,kx_grid,
                             ky_grid,Nq,NL)

energy_diff=big_func(uD,aba_record,abc_record,
        wf,temp,kx_grid,
         ky_grid,Nq,NL,density_list)


scratch_dir = ENV["SCRATCH"]
savepath=joinpath(scratch_dir, "graphene_wf/data_output$(Int(args[9]))/$(args[1])uD$(args[2])Nq$(args[3])temp$(args[4])NL$(args[5])wf$(args[6])denstart$(args[7])denend$(args[8])dennum.jld2")



jldsave(savepath,energy_diff=energy_diff)




