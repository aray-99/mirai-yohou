# 座標変換(SPEC §2・§3)
#
# 有界変数(0〜1)は logit 座標、正値変数は log 座標で状態を保持する。
# xi_tauA は「logit 空間上のアンカー値」で非有界(変換なし)、
# lam_e は自然座標(≥0)のまま保持する。

"""
状態ベクトルのインデックス(SPEC §3。固定であり変更してはならない)
"""
const IX_P = 1      # 総人口 P            (log)
const IX_W = 2      # 生産年齢人口比率 w   (logit)
const IX_H = 3      # 人的資本指数 h       (log)
const IX_K = 4      # 一人当たり資本 k     (log)
const IX_G = 5      # 格差指標 g           (logit)
const IX_T = 6      # 技術フロンティア T   (log)
const IX_PHI = 7    # 技術普及率 φ         (logit)
const IX_V = 8      # 情報伝播速度 v       (log)
const IX_TAU = 9    # 制度信頼 τ           (logit)
const IX_TAUA = 10  # 信頼アンカー τ*      (そのまま: logit 空間の値)
const IX_SIG = 11   # 社会的応力 σ_s       (log)
const IX_PP = 12    # 分極度 p             (logit)
const IX_LAME = 13  # Hawkes 興奮度 λ_e    (自然座標 ≥0)

const N_STATE = 13

const LOG_INDICES = (IX_P, IX_H, IX_K, IX_T, IX_V, IX_SIG)
const LOGIT_INDICES = (IX_W, IX_G, IX_PHI, IX_TAU, IX_PP)
const IDENTITY_INDICES = (IX_TAUA, IX_LAME)

const STATE_NAMES = (:xi_P, :xi_w, :xi_h, :xi_k, :xi_g, :xi_T, :xi_phi,
                     :xi_v, :xi_tau, :xi_tauA, :xi_sig, :xi_p, :lam_e)

# --- 基本変換(SPEC §2.2) ---

logit(x) = log(x / (1 - x))

sigmoid(z) = 1 / (1 + exp(-z))

"数値安定な softplus(常に > 0。単調性の保証に使う)。
分岐レス形はシンボリックトレース(Symbolics、#0020)のため。"
softplus(z) = max(z, zero(z)) + log1p(exp(-abs(z)))

"(a)_+ = max(a, 0)"
pluspart(a) = max(a, zero(a))

# --- 変数別変換 ---

"""
    to_state_var(i, x)

自然座標の第 `i` 変数 `x` を保持座標(状態座標)へ変換する。
"""
function to_state_var(i::Integer, x::Real)
    i in LOG_INDICES && return log(x)
    i in LOGIT_INDICES && return logit(x)
    i in IDENTITY_INDICES && return float(x)
    throw(ArgumentError("invalid state index $i (1..$N_STATE)"))
end

"""
    from_state_var(i, xi)

保持座標の第 `i` 変数 `xi` を自然座標へ戻す。
"""
function from_state_var(i::Integer, xi::Real)
    i in LOG_INDICES && return exp(xi)
    i in LOGIT_INDICES && return sigmoid(xi)
    i in IDENTITY_INDICES && return float(xi)
    throw(ArgumentError("invalid state index $i (1..$N_STATE)"))
end

# --- 一括変換 ---

"""
    to_state(x_nat) -> Vector{Float64}

全13変数を自然座標から保持座標へ一括変換する(SPEC §2.2)。
"""
function to_state(x_nat::AbstractVector{<:Real})
    length(x_nat) == N_STATE ||
        throw(DimensionMismatch("expected $N_STATE variables, got $(length(x_nat))"))
    return [to_state_var(i, x_nat[i]) for i in 1:N_STATE]
end

"""
    from_state(xi) -> Vector{Float64}

全13変数を保持座標から自然座標へ一括変換する。
"""
function from_state(xi::AbstractVector{<:Real})
    length(xi) == N_STATE ||
        throw(DimensionMismatch("expected $N_STATE variables, got $(length(xi))"))
    return [from_state_var(i, xi[i]) for i in 1:N_STATE]
end
