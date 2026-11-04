

using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))

using JLD2
using Random

include(joinpath(@__DIR__, "../../src/Xiaomeng_RMG_sus_v2.jl"))

BLAS.set_num_threads(
    parse(Int, get(ENV, "SLURM_CPUS_PER_TASK", "1"))
)

@assert length(ARGS) == 20
args = parse.(Float64, ARGS)

NL = Int(args[1])
band = Int(args[2])
uD = args[3]
n_h = args[4]
Q = args[5:6]
εr = args[7]
interaction_type = Int(args[8])
Nθ = Int(args[9])
Nr = Int(args[10])
nmodes = Int(args[11])
model_flag = Int(args[12])
filepos = Int(args[13])
rmax = args[14]
compression_tol = args[15]
eig_tol = args[16]
maxrank = Int(args[17])
ncheck = Int(args[18])
nscan = Int(args[19])
leafsize = Int(args[20])

@assert model_flag in (0, 1)
modeltype = model_flag == 1 ? :spinor : :no_spinor

valley = 1
stacking = 1
kmax = 1.0
Nk = 501
r_tol = 1e-10
η = 1.0
krylovdim = max(60, 2 * nmodes + 1)
maxiter = 400

Random.seed!(seed)

# Chemical potential.
mu_data = get_mu(
    uD=uD, n_h=n_h,
    NL=NL, valley=valley, stacking=stacking, band=band,
    kmax=kmax, Nk=Nk
)

# Grid, band data, and compressed exchange blocks.
problem = get_compressed_problem(
    Q, mu_data.μ, uD, valley, stacking, NL, band,
    εr, interaction_type;
    Nθ=Nθ, Nr=Nr, rmax=rmax, nscan=nscan, r_tol=r_tol,
    modeltype=modeltype, leafsize=leafsize,
    compression_tol=compression_tol,
    maxrank=maxrank, η=η, ncheck=ncheck
)

# Lowest eigenmodes.
result = get_compressed_modes(
    problem, nmodes;
    eig_tol=eig_tol, krylovdim=krylovdim, maxiter=maxiter
)

println("μ: ", mu_data.μ)
println("Number of points: ", length(problem.D))
println("Eigenvalues: ", result.eigenvalues)
println("Compressed residuals: ", result.residuals)
println("Weighted normalizations: ", result.normalizations)

# Save.
output_dir = joinpath(
    ENV["SCRATCH"],
    "Xiaomeng_RMG_v2",
    "data_output$(filepos)"
)
mkpath(output_dir)

filename =
    "NL$(NL)_b$(band)_uD$(uD)_nh$(n_h)_" *
    "Qx$(Q[1])_Qy$(Q[2])_er$(εr)_it$(interaction_type)_" *
    "Nt$(Nθ)_Nr$(Nr)_nm$(nmodes)_sp$(model_flag)_R$(rmax)_" *
    "ct$(compression_tol)_et$(eig_tol)_rk$(maxrank)_" *
    "nc$(ncheck)_ns$(nscan)_lf$(leafsize).jld2"

savepath = joinpath(output_dir, filename)
@assert !ispath(savepath) "Output already exists."

jldsave(savepath;
    arguments=args,
    valley=valley,
    stacking=stacking,
    kmax=kmax,
    Nk=Nk,
    r_tol=r_tol,
    η=η,
    krylovdim=krylovdim,
    maxiter=maxiter,

    μ=mu_data.μ,
    n_h_grid=mu_data.n_h_grid,

    Q=problem.Q,
    p_points=problem.p_points,
    weights=problem.weights,
    order=problem.order,
    D=problem.D,
    Uminus=problem.Uminus,
    Uplus=problem.Uplus,
    FQ=problem.FQ,

    eigenvalues=result.eigenvalues,
    modes=result.modes,
    solver_modes=result.solver_modes,
    normalizations=result.normalizations,
    dQ=result.dQ,
    compressed_residuals=result.residuals
)

println("Saved: ", savepath)