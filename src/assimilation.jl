# 逐次・非同期同化ドライバ(SPEC §9.2/§9.3、双子実験 §13 用)
#
# - 観測は届いた時刻に、届いた観測だけで解析(グリッドに丸めた tstops)。
# - イベントは ExogenousEvents 相当: カタログ時刻で全メンバー強制ジャンプ
#   (マーク rho は各メンバー独立)。内生発火は同化中オフ(§9.3)。
# - 無イベント尤度: 週次窓(0.02年 ≈ 7.3日で近似、DECISIONS #0009)で
#   Λ_i = ∫ lam_i dt をトラッキングし、ポアソン重み+ESS<N/2 で系統再抽選。
# - 乗法的インフレーション 1.02 常時、強制ジャンプ直後の解析は 1.05(§9.3)。
# - E1b: theta_sig を log 座標で状態拡大(行 14。d(param)=0+微小ノイズ、§8.3)。

"""
同化の設定(既定値は SPEC §9/§13)。

`inflation_mode`(DECISIONS #0012/#0013):
- `:per_analysis` — 解析毎に偏差 × rho_inf(SPEC v1.0 の原記述。多レート観測
  では弱観測部分空間が複利膨張しフィルタ崩壊する — #0012)
- `:per_time` — 解析時に rho_inf^(Δt_前回解析から / tau_ref)。単位時間あたり
  注入率を解析頻度から切り離す
- `:rtps` — relaxation to prior spread(rtps_alpha)。観測に拘束されない
  成分への注入がゼロで、複利膨張が原理的に起きない
強制ジャンプ直後の「一時的に強める」(§9.3)は rho_inf_jump / rtps_alpha_jump。
"""
Base.@kwdef struct AssimConfig
    t0::Float64 = 0.0
    t1::Float64 = 45.0
    dt::Float64 = 0.01
    inflation_mode::Symbol = :rtps   # 採用方式(#0013)。:per_analysis は SPEC v1.0 原案
    rho_inf::Float64 = 1.02          # 乗法モードの基礎レート(§9.2)
    rho_inf_jump::Float64 = 1.05     # 強制ジャンプ直後(§9.3)
    tau_ref::Float64 = 0.25          # :per_time の正規化時定数(四半期)
    rtps_alpha::Float64 = 0.7        # :rtps の緩和係数(#0013 の診断マトリクスで選定)
    rtps_alpha_jump::Float64 = 0.8   # 強制ジャンプ直後の :rtps 係数
    event_window::Float64 = 0.02     # 週次バッチ近似(#0009)
    ess_ratio::Float64 = 0.5         # ESS < N * ratio で再抽選(§9.3)
    param_noise_sd::Float64 = 0.01   # 状態拡大パラメータの微小ノイズ(/√年)
end

"同化ランの結果(X は状態行 × 時刻 × メンバー。拡大時は最終行がパラメータ)"
struct AssimResult
    t::Vector{Float64}
    X::Array{Float64,3}
    ranks::Dict{Symbol,Vector{Int}}      # 解析直前の順位(ランクヒストグラム用)
    ess::Vector{Float64}                 # 各週次窓の ESS
    nresample::Int
end

"""
    pathological(xi) -> Bool

メンバーが数値的に病的な領域にいるか(DECISIONS #0011)。
§9.4 の警告水準(|logit ξ| > 10)を大きく超えた |logit ξ| > 15、
または σ_s > e³ ≈ 20(降伏応力 σ_Y = 1 の20倍)、または非有限値。
"""
function pathological(xi::AbstractVector{Float64})
    all(isfinite, xi) || return true
    for i in (IX_W, IX_G, IX_PHI, IX_TAU, IX_PP)
        abs(xi[i]) > 15 && return true
    end
    return xi[IX_SIG] > 3
end

"theta_sig を差し替えた ModelParameters(状態拡大メンバー用)"
with_theta_sig(p::ModelParameters, theta::Real) =
    ModelParameters(p.regime, p.l1, p.l2, L3Params(theta_sig = float(theta)),
                    p.exo, p.x0_nat, p.x0)

"""
    run_assimilation(params, E0, obs, event_times; cfg, seed,
                     augmented=false) -> AssimResult

初期アンサンブル `E0`(n × N。augmented なら n = 14 で最終行 = log theta_sig)
から §9 のハイブリッド同化(EnKF + ポアソン重み + イベント同期)を実行する。
週次イベントカウントは `event_times`(真値カタログ)を窓に集計して用いる。
"""
function run_assimilation(params::ModelParameters, E0::Matrix{Float64},
                          obs::Vector{ObservationRecord},
                          event_times::Vector{Float64};
                          cfg::AssimConfig = AssimConfig(), seed::Integer,
                          augmented::Bool = false)
    n, N = size(E0)
    n == (augmented ? N_STATE + 1 : N_STATE) ||
        throw(DimensionMismatch("E0 has $n rows, augmented=$augmented"))

    nsteps = round(Int, (cfg.t1 - cfg.t0) / cfg.dt)
    ts = collect(range(cfg.t0; step = cfg.dt, length = nsteps + 1))
    grid_index(t) = clamp(round(Int, (t - cfg.t0) / cfg.dt) + 1, 1, nsteps + 1)

    # 観測をグリッド点にグループ化(§9.2: 届いた時刻に届いた観測だけで解析)
    obs_at = Dict{Int,Vector{ObservationRecord}}()
    for o in obs
        push!(get!(obs_at, grid_index(o.t), ObservationRecord[]), o)
    end

    events = sort(event_times)
    next_ev = 1

    # 週次窓の境界グリッドと観測カウント
    wsteps = max(1, round(Int, cfg.event_window / cfg.dt))

    E = copy(E0)
    X = Array{Float64,3}(undef, n, nsteps + 1, N)
    X[:, 1, :] = E
    rngs = [Xoshiro(member_seed(seed, i)) for i in 1:N]
    f = Vector{Float64}(undef, N_STATE)
    sig = Vector{Float64}(undef, N_STATE)
    dW = Vector{Float64}(undef, N_STATE)
    Lambda = zeros(N)
    logw = zeros(N)          # 累積 log 重み(再抽選までウィンドウ間で持ち越す)
    window_count = 0
    ess_hist = Float64[]
    nresample = 0
    jump_since_analysis = false
    t_last_analysis = cfg.t0
    ranks = Dict{Symbol,Vector{Int}}()
    sqdt = sqrt(cfg.dt)

    member_params(i) = augmented ? with_theta_sig(params, exp(E[end, i])) : params

    for step in 1:nsteps
        t = ts[step]
        t_next = ts[step + 1]

        # (a) 強制ジャンプ(イベント同期、§9.3)
        while next_ev <= length(events) && events[next_ev] < t_next
            for i in 1:N
                xi = @view E[1:N_STATE, i]
                rho = draw_mark(rngs[i], params)
                apply_jump!(xi, rho, member_params(i))
            end
            window_count += 1
            jump_since_analysis = true
            next_ev += 1
        end

        # (b) Λ トラッキング(§9.3。直前状態で近似)と EM ステップ
        for i in 1:N
            p_i = member_params(i)
            xi = @view E[1:N_STATE, i]
            Lambda[i] += intensity(xi, p_i) * cfg.dt
            drift!(f, xi, p_i, t)
            diffusion!(sig, xi, p_i, t)
            randn!(rngs[i], dW)
            @. xi += cfg.dt * f + sqdt * sig * dW
            if augmented   # d(param) = 0 + 微小ノイズ(§8.3)
                E[end, i] += cfg.param_noise_sd * sqdt * randn(rngs[i])
            end
        end

        # (c) 週次窓の終端: ポアソン重みを累積し、ESS < N/2 で系統再抽選(§9.3)
        if step % wsteps == 0
            logw .+= poisson_logweights(window_count, Lambda)
            # 病的メンバーは重みゼロ化して強制再抽選(#0011)。ESS は単一
            # 外れ値では下がらないため、暴走メンバーが強制ジャンプ
            # (m ∝ sigma_s^-)で数値爆発する前に淘汰する必要がある。
            npath = 0
            for i in 1:N
                if pathological(view(E, 1:N_STATE, i))
                    logw[i] = -Inf
                    npath += 1
                end
            end
            npath < N || error("filter diverged: all members pathological")
            w = normalize_weights(logw)
            essval = ess(w)
            push!(ess_hist, essval)
            if npath > 0 || essval < N * cfg.ess_ratio
                idx = systematic_resample(rngs[1], w)
                E .= E[:, idx]
                fill!(logw, 0.0)
                nresample += 1
            end
            fill!(Lambda, 0.0)
            window_count = 0
        end

        # (d) 解析ステップ(この時刻に届いた観測のみ、§9.2)
        if haskey(obs_at, step + 1)
            batch = obs_at[step + 1]
            # ランク(解析直前の事前アンサンブルに対する観測の順位)
            for o in batch
                yj = [o.spec.h(view(E, 1:N_STATE, j)) for j in 1:N]
                push!(get!(ranks, o.spec.name, Int[]),
                      count(<(o.value), yj) + 1)
            end
            yobs = [o.value for o in batch]
            R = Diagonal([o.spec.sd^2 for o in batch]) |> Matrix
            hfun = col -> [o.spec.h(view(col, 1:N_STATE)) for o in batch]

            # スプレッド注入(inflation_mode、DECISIONS #0013)
            if cfg.inflation_mode === :rtps
                spread_prior = ensemble_spread(E)
                enkf_analysis!(E, yobs, hfun, R; rng = rngs[1], rho_inf = 1.0)
                alpha = jump_since_analysis ? cfg.rtps_alpha_jump : cfg.rtps_alpha
                rtps!(E, spread_prior; alpha)
            else
                rho_base = jump_since_analysis ? cfg.rho_inf_jump : cfg.rho_inf
                rho = cfg.inflation_mode === :per_time ?
                    rho_base^((t_next - t_last_analysis) / cfg.tau_ref) : rho_base
                enkf_analysis!(E, yobs, hfun, R; rng = rngs[1], rho_inf = rho)
            end
            postprocess_analysis!(E)
            jump_since_analysis = false
            t_last_analysis = t_next
        end

        X[:, step + 1, :] = E
    end
    return AssimResult(ts, X, ranks, ess_hist, nresample)
end

"""
    free_ensemble(params, E0; cfg, seed) -> Array{Float64,3}

同化オフの自由ラン(同じ初期アンサンブル、内生 Hawkes、§13 手順5)。
戻り値は N_STATE × 時刻 × メンバー。
"""
function free_ensemble(params::ModelParameters, E0::Matrix{Float64};
                       cfg::AssimConfig = AssimConfig(), seed::Integer)
    N = size(E0, 2)
    nsteps = round(Int, (cfg.t1 - cfg.t0) / cfg.dt)
    X = Array{Float64,3}(undef, N_STATE, nsteps + 1, N)
    Threads.@threads for i in 1:N
        r = simulate_sde(params; seed = member_seed(seed, i),
                         t0 = cfg.t0, t1 = cfg.t1, dt = cfg.dt,
                         xi0 = E0[1:N_STATE, i])
        X[:, :, i] = r.traj.X
    end
    return X
end
