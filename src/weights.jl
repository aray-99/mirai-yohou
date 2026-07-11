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
    poisson_logpmf(N_obs, mu) -> Float64

素のポアソン log pmf(log N! 定数込み。予測 log スコア用 — 重み計算には
定数省略版の `poisson_logweights` を使う)。
"""
poisson_logpmf(N_obs::Integer, mu::Real) =
    logpdf(Poisson(max(float(mu), 1e-10)), N_obs)

"""
    negbin_logpmf(N_obs, mu, r) -> Float64

負の二項分布の log pmf(DECISIONS #0054)。平均 `mu`、サイズ(分散)
パラメータ `r` の (mu, r) パラメタライズ: Var = mu + mu²/r。
`r → ∞` でポアソンに収束する。Distributions の NegativeBinomial(r, p)
(p = r/(r+mu))に委譲。
"""
negbin_logpmf(N_obs::Integer, mu::Real, r::Real) =
    logpdf(NegativeBinomial(float(r), float(r) / (float(r) + max(float(mu), 1e-10))),
           N_obs)

"""
    negbin_logweights(N_obs, mus, r) -> Vector{Float64}

NegBin 尤度のメンバー別 log 重み(#0054)。r・N_obs はメンバー共通なので
定数項は省略可能だが、log スコアとの一貫性のため完全な log pmf を返す
(重みの正規化には影響しない)。
"""
negbin_logweights(N_obs::Integer, mus::AbstractVector{<:Real}, r::Real) =
    [negbin_logpmf(N_obs, mu, r) for mu in mus]

"""
    negbin_profile_r(counts, mus; lo=1e-2, hi=1e4, iters=200) -> Float64

サイズパラメータ r のプロファイル最尤(#0054): 平均列 `mus`(= ν·Λ_w、
較正窓のみ)を固定して Σ_w log NegBin(N_w | mu_w, r) を r について最大化
する。log r 空間の黄金分割探索(単峰。境界解 hi 到達 ≈ ポアソンで十分の
サイン)。ν は現行どおり別途プロファイル(ν* = ΣN/ΣΛ、#0031-1/#0033)し、
本関数では動かさない。
"""
function negbin_profile_r(counts::AbstractVector{<:Integer},
                          mus::AbstractVector{<:Real};
                          lo::Float64 = 1e-2, hi::Float64 = 1e4, iters::Int = 200)
    length(counts) == length(mus) ||
        throw(DimensionMismatch("counts and mus must have equal length"))
    isempty(counts) && throw(ArgumentError("no count windows to profile r"))
    ll(logr) = sum(negbin_logpmf(c, mu, exp(logr)) for (c, mu) in zip(counts, mus))
    a, b = log(lo), log(hi)
    phi = (sqrt(5.0) - 1) / 2
    c = b - phi * (b - a); d = a + phi * (b - a)
    fc, fd = ll(c), ll(d)
    for _ in 1:iters
        if fc > fd
            b, d, fd = d, c, fc
            c = b - phi * (b - a); fc = ll(c)
        else
            a, c, fc = c, d, fd
            d = a + phi * (b - a); fd = ll(d)
        end
        (b - a) < 1e-8 && break
    end
    return exp((a + b) / 2)
end

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
