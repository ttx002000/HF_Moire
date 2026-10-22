

using Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using JLD2

include("../../src/Xiaomeng_RMG_sus.jl")

BLAS.set_num_threads(parse(Int, get(ENV, "SLURM_CPUS_PER_TASK", "1")))

argument_names = [
    "NL", "band", "uD", "n_h",
    "Qx", "Qy", "epsilon_r", "dgate",
    "Rpair", "dk_pair", "nmodes", "model_flag", "filepos",
]

@assert length(ARGS) == length(argument_names) (
    "Usage: julia main.jl " * join(argument_names, " ")
)

args = parse.(Float64, ARGS)
@assert all(isfinite, args) "All arguments must be finite."

NL         = Int(args[1])
band       = Int(args[2])
uD         = args[3]
n_h        = args[4]
Qx         = args[5]
Qy         = args[6]
εr         = args[7]
dgate      = args[8]
Rpair      = args[9]
dk_pair    = args[10]
nmodes     = Int(args[11])
model_flag = Int(args[12])
filepos    = Int(args[13])

# Fixed parameters
valley = 1
stacking = 1

# Default grid parameters for determining μ
kmax = 1.0
Nk = 501

Q = [Qx, Qy]
modeltype = model_flag == 1 ? :spinor : :no_spinor

@assert NL >= 1 && 1 <= band <= 2NL
@assert n_h > 0 && norm(Q) > 0
@assert εr > 0 && dgate > 0
@assert Rpair > 0 && dk_pair > 0
@assert nmodes > 0
@assert model_flag in (0, 1)
@assert filepos >= 0


filename = join([
    "NL$(NL)", "v$(valley)", "s$(stacking)", "b$(band)",
    "uD$(uD)", "nh$(n_h)", "Qx$(Qx)", "Qy$(Qy)",
    "eps$(εr)", "gate$(dgate)", "R$(Rpair)", "dk$(dk_pair)",
    "kmax$(kmax)", "Nk$(Nk)", "nm$(nmodes)", "spin$(model_flag)",
], "_") * ".jld2"


# Output directory
scratch_dir = get(ENV, "SCRATCH", @__DIR__)
output_dir = joinpath(scratch_dir, "Xiaomeng_RMG", "data_output$(filepos)")
mkpath(output_dir)

savepath = joinpath(output_dir, filename)

@assert !ispath(savepath) (
    "Output already exists. Choose another filepos for a separate run."
)


# Determine μ
println("Computing chemical potential...")

mu_data = get_mu(
    uD=uD,
    n_h=n_h,
    NL=NL,
    valley=valley,
    stacking=stacking,
    band=band,
    kmax=kmax,
    Nk=Nk,
)

μ = mu_data.μ

println("μ (meV): ", μ)
println("Sampled hole density: ", mu_data.n_h_grid)


# Construct the active pairs
println("Constructing active pairs...")

pair = get_pair_grid(
    Q, μ, uD, valley, stacking, NL, band,
    Rpair, dk_pair,
)

println("Number of active pairs: ", pair.Nactive)
println("HQ storage alone (GiB): ", 16.0 * pair.Nactive^2 / 2.0^30)


# Calculate the endpoint spinors
println("Computing endpoint spinors...")

spinors = get_pair_spinors(
    pair.k_active, pair.kQ_active,
    uD, valley, stacking, NL, band, modeltype,
)


# Construct H^Q
println("Constructing HQ...")

HQ = get_HamQ(
    pair.k_active, Q, pair.D,
    spinors.U, spinors.UQ, spinors.FQ,
    pair.w_pair, εr, dgate,
)


# Extract the lowest eigenmodes
println("Computing lowest eigenmodes...")

result = get_lowest_modes(
    HQ, spinors.FQ, pair.w_pair, nmodes,
)

println("Eigenvalues (meV): ", result.eigenvalues)
println("Residuals (meV): ", result.residuals)
println("Normalizations: ", result.normalizations)
println("Charge overlaps: ", result.dQ)


jldsave(savepath;
    # Input parameters
    arguments=args,
    argument_names=argument_names,
    NL=NL,
    valley=valley,
    stacking=stacking,
    band=band,
    uD=uD,
    n_h=n_h,
    Qx=Qx,
    Qy=Qy,
    Q=Q,
    εr=εr,
    dgate=dgate,
    Rpair=Rpair,
    dk_pair=dk_pair,
    kmax=kmax,
    Nk=Nk,
    nmodes=nmodes,
    model_flag=model_flag,
    modeltype=String(modeltype),
    filepos=filepos,

    # Chemical potential and its sampling grid
    μ=μ,
    n_h_grid=mu_data.n_h_grid,
    mu_kgrid=collect(mu_data.kgrid),
    mu_εgrid=mu_data.εgrid,
    mu_fgrid=mu_data.fgrid,
    mu_dk=mu_data.dk,
    mu_w=mu_data.w,
    mu_Nholes=mu_data.Nholes,

    # Full pair grid
    eQ=pair.eQ,
    eperp=pair.eperp,
    sgrid=pair.sgrid,
    tgrid=pair.tgrid,
    w_pair=pair.w_pair,
    kmesh=pair.kmesh,
    kQmesh=pair.kQmesh,
    εpair=pair.εpair,
    εpairQ=pair.εpairQ,
    fpair=pair.fpair,
    fpairQ=pair.fpairQ,

    # Active pairs and indexing
    active=pair.active,
    active_index=pair.active_index,
    index_grid=pair.index_grid,
    k_active=pair.k_active,
    kQ_active=pair.kQ_active,
    D=pair.D,
    Nactive=pair.Nactive,

    # Spinors and form factors
    U=spinors.U,
    UQ=spinors.UQ,
    FQ=spinors.FQ,

    # Eigenmodes and diagnostics
    eigenvalues=result.eigenvalues,
    modes=result.modes,
    dQ=result.dQ,
    normalizations=result.normalizations,
    residuals=result.residuals,
    hermiticity_error=result.hermiticity_error,
    nmodes_saved=length(result.eigenvalues),
)

println("Saved: ", savepath)