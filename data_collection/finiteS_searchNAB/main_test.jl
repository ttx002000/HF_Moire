using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using Plots
using JLD2
using CSV,DataFrames
args=parse.(Float64,ARGS)

include("../../src/single_particle_finiteS.jl")
scale=1.0
standard=10.01*scale
vf=args[1]
spin=args[2]
V0=args[3]
ϕ=args[4]/180*π

am=4*π/(√3*scale)


g1=scale*[1,0]
g2=scale*[1/2,√3/2]
g3=scale*[-1/2,√3/2]
g4=-g1
g5=-g2
g6=-g3

#a1m=am*[√3/2,1/2]
#a2m=am*[0,1]

a1m=am*[√3/2,-1/2]
a2m=am*[√3/2,1/2]
Km=(g1+g6)/3
Kp=(g1+g2)/3
γ=[0,0];
mpoint=g1/2



l1=[4,0]
l2=[0,4]
Nx=4;
Ny=4;
L1=l1[1]*a1m+l1[2]*a2m;
L2=l2[1]*a1m+l2[2]*a2m;
area=abs(L1[1]*L2[2]-L2[1]*L1[2]);
Rotminus90=[0 1;-1 0]
T1=2*π/area*Rotminus90*L2
T2=-2*π/area*Rotminus90*L1
g1T=Int.(round.(inv([T1 T2])*g1))
g3T=Int.(round.(inv([T1 T2])*g3))


allowedq=Vector{Int64}[]
for jb in 1:Ny, ja in 1:Nx
push!(allowedq,[ja-1,jb-1])
end

Nqpath=60
Trqpath=Vector{Vector{Float64}}(undef,Nqpath)
for ja in 1:Nqpath
   Trqpath[ja]=(1-ja/Nqpath)*γ+(ja/Nqpath)*Km*2
end



wave=getwave(standard,g1,g3,g1T,g3T)
Moireglist=getMoireglist(wave,g1T,g3T)
formq=Vector{Vector{Float64}}(undef,Nx*Ny)
averaged_formq=zeros(Float64,Nqpath)

for ja in 1:Nx*Ny
    formq[ja]=calculateformq(allowedq[ja][1]*T1+allowedq[ja][2]*T2,Trqpath,wave,Moireglist,T1,T2,g1,g3,spin,vf,V0,ϕ)
end


for ja in 1:Nqpath,jb in 1:Nx*Ny
    averaged_formq[ja]+=formq[jb][ja]/(Nx*Ny)
end

formq1LL=zeros(Float64,Nqpath)

Lb=(4*3^(1/2)*π/3)^(1/2)
for ja in 1:Nqpath
    formq1LL[ja]=exp(-Lb^2*norm(Trqpath[ja])^2/2)*(1-Lb^2*norm(Trqpath[ja])^2/2)^2
end

diff=sum(abs.(formq1LL-averaged_formq))
diffdivedeq=0
for ja in 1:Nqpath
    diffdivideq+=abs(formq1LL[ja]-averaged_formq[ja])/norm(Trqpath[ja])
end


chern,uniform,Flink=calculatechern(wave,Moireglist,T1,T2,g1,g3,spin,vf,V0,ϕ)




jldsave(joinpath(@__DIR__, "data_output/$(args[1])vf$(args[2])spin$(args[3])V0$(args[4])phi.jld2"),chern=chern,Flink=Flink,arguments=args,uniform=uniform,diff=diff,diffdivedeq=diffdivedeq)


