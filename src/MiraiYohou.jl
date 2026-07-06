"""
社会動態シミュレータ「未来予報」

仕様は docs/SPEC.md(Single Source of Truth)。実装判断は docs/DECISIONS.md。
"""
module MiraiYohou

using Distributions
using Random
using TOML

include("coordinates.jl")
include("parameters.jl")
include("diagnostics.jl")
include("drift.jl")
include("diffusion.jl")
include("jumps.jl")
include("integrator.jl")

# coordinates (§2/§3)
export logit, sigmoid, softplus, pluspart
export to_state, from_state, to_state_var, from_state_var
export N_STATE, STATE_NAMES
export IX_P, IX_W, IX_H, IX_K, IX_G, IX_T, IX_PHI, IX_V,
       IX_TAU, IX_TAUA, IX_SIG, IX_PP, IX_LAME

# parameters (§8)
export L1Params, L2Params, L3Params, ConstantExogenous, ModelParameters
export build_params, prior_lognormal, l3_priors

# diagnostics (§4/§7)
export DimensionlessNumbers, dimensionless_numbers
export branching_ratio, deborah_number, hardening_ratio
export tfp, output_y, tech_growth, dep_ratio

# drift / diffusion / jumps / integrator (§5/§6/§10)
export drift!, drift_with_diagnostics!, diffusion!
export JumpMode, EndogenousHawkes, ExogenousEvents
export lam_b, intensity, draw_mark, apply_jump!, JumpEvent, simulate_hawkes
export Trajectory, simulate_ode, SDEResult, simulate_sde

end # module
