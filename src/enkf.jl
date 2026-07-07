# EnKF 解析ステップ(SPEC §9.2/§9.4)
#
# 摂動観測型 stochastic EnKF(自作)。接線行列は使わず、H は関数として
# 与え、アンサンブル共分散で線形化を代用する。解析後は乗法的
# インフレーション(既定 rho_inf = 1.02)を常時適用する(§9.2)。

"""
    enkf_analysis!(E, yobs, hfun, R; rng, rho_inf=1.02) -> E

摂動観測 EnKF の解析更新(in-place)。

- `E`: 状態アンサンブル(n × N。列がメンバー)
- `yobs`: 観測ベクトル(m)
- `hfun`: 状態列 → 観測空間ベクトル(m)の関数
- `R`: 観測誤差共分散(m × m)
- `rho_inf`: 乗法的インフレーション(解析後、平均まわりの偏差を拡大)

K = Cxy (Cyy + R)^{-1} をアンサンブル共分散から構成し、各メンバーに
独立な観測摂動 eps_j ~ N(0, R) を加えて更新する。
"""
function enkf_analysis!(E::AbstractMatrix{Float64}, yobs::AbstractVector{Float64},
                        hfun, R::AbstractMatrix{Float64};
                        rng::AbstractRNG, rho_inf::Float64 = 1.02)
    N = size(E, 2)
    m = length(yobs)
    N > 1 || throw(ArgumentError("ensemble size must be > 1"))

    Yf = Matrix{Float64}(undef, m, N)
    for j in 1:N
        Yf[:, j] = hfun(view(E, :, j))
    end

    xbar = sum(E; dims = 2) ./ N
    ybar = sum(Yf; dims = 2) ./ N
    Xp = E .- xbar
    Yp = Yf .- ybar

    Cyy = (Yp * Yp') ./ (N - 1) .+ R
    Cxy = (Xp * Yp') ./ (N - 1)
    K = Cxy / Cyy                                   # n × m

    Rchol = cholesky(Symmetric(R)).L
    innov = Vector{Float64}(undef, m)
    for j in 1:N
        innov .= yobs .+ Rchol * randn(rng, m) .- @view Yf[:, j]
        E[:, j] .+= K * innov
    end

    # 乗法的インフレーション(§9.2。解析後の平均まわり)
    if rho_inf != 1.0
        xbar2 = sum(E; dims = 2) ./ N
        @. E = xbar2 + rho_inf * (E - xbar2)
    end
    return E
end

"""
    rtps!(E, spread_prior; alpha) -> E

RTPS(relaxation to prior spread、§9.2 改訂・DECISIONS #0013)。
解析後のアンサンブル偏差の成分別スプレッド σ_post を、解析前の
スプレッド σ_prior へ割合 `alpha` だけ緩和する:
新偏差 = 偏差 × (σ_post + α(σ_prior − σ_post)) / σ_post。
観測に拘束されない成分は σ_post = σ_prior のため変化せず、
乗法的インフレーションのような弱観測部分空間の複利膨張が起きない。
"""
function rtps!(E::AbstractMatrix{Float64}, spread_prior::AbstractVector{Float64};
               alpha::Float64)
    alpha == 0.0 && return E
    N = size(E, 2)
    xbar = sum(E; dims = 2) ./ N
    for i in axes(E, 1)
        s_post = sqrt(sum(abs2, @view(E[i, :]) .- xbar[i]) / (N - 1))
        s_post > 0 || continue
        c = (s_post + alpha * (spread_prior[i] - s_post)) / s_post
        for j in axes(E, 2)
            E[i, j] = xbar[i] + c * (E[i, j] - xbar[i])
        end
    end
    return E
end

"成分別アンサンブルスプレッド(RTPS の事前保存用)"
function ensemble_spread(E::AbstractMatrix{Float64})
    N = size(E, 2)
    xbar = sum(E; dims = 2) ./ N
    return [sqrt(sum(abs2, @view(E[i, :]) .- xbar[i]) / (N - 1)) for i in axes(E, 1)]
end

"スカラー観測での解析更新(逐次・非同期同化 §9.2 の1ステップ)"
enkf_analysis!(E::AbstractMatrix{Float64}, yobs::Real, hfun_scalar, sd::Real;
               rng::AbstractRNG, rho_inf::Float64 = 1.02) =
    enkf_analysis!(E, [float(yobs)], xi -> [hfun_scalar(xi)],
                   fill(sd^2, 1, 1); rng, rho_inf)

"""
    enks_analysis!(E, snapshots, yobs, hfun, R; rng, rho_inf=1.0) -> E

固定ラグ EnKS の解析更新(SPEC 追補・DECISIONS #0024)。現在状態 `E` を
摂動観測 EnKF で更新し、**同一の観測摂動・イノベーション**を用いて過去の
アンサンブルスナップショット `snapshots`(各 n × N、列=メンバー)も
クロス共分散 K_s = Cov(x_past, y) (Cyy + R)^{-1} で更新する。
`snapshots` が空なら enkf_analysis! と同一の更新になる。
`smooth_rows` を与えると、スナップショットの該当行のみを更新する
(変数局所化。実情報のない変数への標本雑音蓄積を防ぐ、#0025)。
"""
function enks_analysis!(E::AbstractMatrix{Float64},
                        snapshots::AbstractVector{<:AbstractMatrix{Float64}},
                        yobs::AbstractVector{Float64}, hfun,
                        R::AbstractMatrix{Float64};
                        rng::AbstractRNG, rho_inf::Float64 = 1.0,
                        smooth_rows::Union{Colon,Vector{Int}} = Colon())
    N = size(E, 2)
    m = length(yobs)
    N > 1 || throw(ArgumentError("ensemble size must be > 1"))

    Yf = Matrix{Float64}(undef, m, N)
    for j in 1:N
        Yf[:, j] = hfun(view(E, :, j))
    end
    ybar = sum(Yf; dims = 2) ./ N
    Yp = Yf .- ybar
    Cyy = (Yp * Yp') ./ (N - 1) .+ R

    # メンバー別イノベーション(観測摂動込み)を一度だけ引き、全ブロックで共有
    Rchol = cholesky(Symmetric(R)).L
    innov = Matrix{Float64}(undef, m, N)
    for j in 1:N
        innov[:, j] = yobs .+ Rchol * randn(rng, m) .- @view Yf[:, j]
    end

    # 現在状態の更新
    xbar = sum(E; dims = 2) ./ N
    K = ((E .- xbar) * Yp') ./ (N - 1) / Cyy
    E .+= K * innov
    if rho_inf != 1.0
        xbar2 = sum(E; dims = 2) ./ N
        @. E = xbar2 + rho_inf * (E - xbar2)
    end

    # 過去スナップショットの平滑化更新(インフレーションなし、行は smooth_rows)
    for S in snapshots
        Ssel = view(S, smooth_rows, :)
        sbar = sum(Ssel; dims = 2) ./ N
        Ks = ((Ssel .- sbar) * Yp') ./ (N - 1) / Cyy
        Ssel .+= Ks * innov
    end
    return E
end

"""
    postprocess_analysis!(E) -> E

解析ステップ直後の後処理(SPEC §3/§9.4):
lam_e = max(lam_e, 0) のクランプ、|logit 座標| > 10 の警告ログ。
状態が13次元モデル状態(+拡大パラメータ)の場合に用いる。
"""
function postprocess_analysis!(E::AbstractMatrix{Float64})
    for j in axes(E, 2)
        E[IX_LAME, j] = max(E[IX_LAME, j], 0.0)
    end
    for i in (IX_W, IX_G, IX_PHI, IX_TAU, IX_PP)
        mx = maximum(abs, @view E[i, :])
        if mx > 10
            @warn "logit coordinate out of range after analysis" index = i maxabs = mx
        end
    end
    return E
end
