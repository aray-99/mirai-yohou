"""
社会動態シミュレータ「未来予報」

仕様は docs/SPEC.md(Single Source of Truth)。実装判断は docs/DECISIONS.md。
"""
module MiraiYohou

using Distributions
using Random
using TOML

include("coordinates.jl")

# coordinates (§2/§3)
export logit, sigmoid, softplus, pluspart
export to_state, from_state, to_state_var, from_state_var
export N_STATE, STATE_NAMES
export IX_P, IX_W, IX_H, IX_K, IX_G, IX_T, IX_PHI, IX_V,
       IX_TAU, IX_TAUA, IX_SIG, IX_PP, IX_LAME

end # module
