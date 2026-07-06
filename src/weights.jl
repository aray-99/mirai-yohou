# ポアソン重み・ESS・系統リサンプリング(SPEC §9.3)
#
# 各解析窓 [t1, t2] でメンバー i の重み w_i ∝ Poisson(N_obs | Λ_i)、
# Λ_i = ∫ lam_i dt(数値積分でトラッキング)。
# ESS = 1 / Σ w_i² が N/2 を下回ったら系統リサンプリング。

"""
    poisson_logweights(N_obs, Lambdas) -> Vector{Float64}

log w_i = N_obs * log(Λ_i) - Λ_i(共通定数 log N_obs! は省略)。
"""
poisson_logweights(N_obs::Integer, Lambdas::AbstractVector{<:Real}) =
    [N_obs * log(L) - L for L in Lambdas]

"""
    normalize_weights(logw) -> Vector{Float64}

log-sum-exp で正規化した重み(Σ w = 1)。
"""
function normalize_weights(logw::AbstractVector{<:Real})
    mx = maximum(logw)
    w = exp.(logw .- mx)
    return w ./ sum(w)
end

"有効サンプルサイズ ESS = 1 / Σ w_i²(§9.3)"
ess(w::AbstractVector{<:Real}) = 1 / sum(abs2, w)

"""
    systematic_resample(rng, w) -> Vector{Int}

系統リサンプリング。正規化済み重み `w` から N 個のインデックスを返す
(期待複製数 N*w_i、分散最小のグリッドサンプリング)。
"""
function systematic_resample(rng::AbstractRNG, w::AbstractVector{<:Real})
    N = length(w)
    u0 = rand(rng) / N
    cum = cumsum(w)
    idx = Vector{Int}(undef, N)
    j = 1
    for k in 1:N
        u = u0 + (k - 1) / N
        while cum[j] < u && j < N
            j += 1
        end
        idx[k] = j
    end
    return idx
end

"""
    resample_if_needed!(E, Lambdas, N_obs; rng, threshold_ratio=0.5)
        -> (w, essval, resampled::Bool)

ポアソン重みを計算し、ESS < N * threshold_ratio なら `E`(列=メンバー)を
系統リサンプリングで置き換える(§9.3)。
"""
function resample_if_needed!(E::AbstractMatrix{Float64},
                             Lambdas::AbstractVector{<:Real}, N_obs::Integer;
                             rng::AbstractRNG, threshold_ratio::Float64 = 0.5)
    N = size(E, 2)
    length(Lambdas) == N ||
        throw(DimensionMismatch("Lambdas length must equal ensemble size"))
    w = normalize_weights(poisson_logweights(N_obs, Lambdas))
    essval = ess(w)
    resampled = essval < N * threshold_ratio
    if resampled
        idx = systematic_resample(rng, w)
        E .= E[:, idx]
    end
    return (w, essval, resampled)
end
