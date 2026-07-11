"""
社会動態シミュレータ「未来予報」

仕様は docs/SPEC.md(Single Source of Truth)。実装判断は docs/DECISIONS.md。
"""
module MiraiYohou

using Distributions
using LinearAlgebra
using Random
using TOML

# Phase B(§10): SciML スタック(名前空間つきで使用)
import JumpProcesses
import SciMLBase
import StochasticDiffEq
import Symbolics
using SciMLBase: SDEProblem

include("coordinates.jl")
include("parameters.jl")
include("diagnostics.jl")
include("drift.jl")
include("diffusion.jl")
include("jumps.jl")
include("integrator.jl")
include("observation.jl")
include("enkf.jl")
include("weights.jl")
include("assimilation.jl")
include("phaseb.jl")

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
export EnsembleResult, simulate_ensemble, member_seed

# observation / enkf / weights (§9)
export ObservationSpec, ObservationRecord
export standard_observations, observation_times, synthesize_observations
export enkf_analysis!, enks_analysis!, postprocess_analysis!, rtps!, ensemble_spread
export poisson_logweights, normalize_weights, ess, systematic_resample,
       resample_if_needed!

# assimilation driver (§9.2/§9.3/§13)
export AssimConfig, AssimResult, run_assimilation, free_ensemble, with_theta_sig
export AugmentedParam, augment_ensemble, build_member_params, simulate_sde_augmented

# Phase B (§10)
export simulate_sde_phaseb, drift_jacobian_sparsity

end # module
