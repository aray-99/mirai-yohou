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
    smoother_lag::Float64 = 0.0      # 固定ラグ EnKS のラグ(年)。0 = 平滑化オフ(#0024)
    smoother_dt::Float64 = 0.1       # 平滑化スナップショットの間隔(年)
    # 平滑化更新は period ≥ この値の観測を含む解析のみで行う。高頻度観測は
    # 過去状態への実情報が乏しく、クロス共分散の標本雑音だけが累積するため(#0025)
    smoother_min_period::Float64 = 2.0
    # 平滑化で更新する状態行(変数局所化、#0025)。実情報のない変数
    # (k 等)への雑音蓄積を防ぐ。既定は制度ブロック + 格差。
    smoother_vars::Vector{Int} = [IX_G, IX_TAU, IX_TAUA, IX_SIG, IX_PP]
    # tauA(IX_TAUA)への緩い擬似観測の倍率(DECISIONS #0036)。tau 観測と
    # 同時刻・同値の擬似観測を sd = このスカラー × tau 観測 sd で追加する。
    # 既定 0.0 = オフ(従来動作。E1・既存テストの記録結果を保護)。
    tauA_pseudo_sd_mult::Float64 = 0.0
    # 現在時刻解析の変数局所化(DECISIONS #0040-(α))。ここに列挙した状態行は、
    # 解析バッチが analysis_unmask_names のいずれの観測名も含まない場合、
    # 現在時刻の EnKF 更新(K の該当行)をマスクする。EnKS の smoother_vars/
    # smooth_rows(過去平滑化の局所化、#0025)と対になる現在時刻側の局所化。
    # 既定は空 = 従来動作(後方互換)。
    analysis_masked_vars::Vector{Int} = Int[]
    # analysis_masked_vars のマスクを解除する観測名(#0040-(α))。解析バッチに
    # この name の観測が1つでも含まれればマスクを解除する。
    analysis_unmask_names::Vector{Symbol} = Symbol[]
end

"同化ランの結果(X は状態行 × 時刻 × メンバー。拡大時は最終行がパラメータ)"
struct AssimResult
    t::Vector{Float64}
    X::Array{Float64,3}
    ranks::Dict{Symbol,Vector{Int}}      # 解析直前の順位(ランクヒストグラム用)
    ess::Vector{Float64}                 # 各週次窓の ESS
    nresample::Int
    ts_snap::Vector{Float64}             # EnKS スナップショット時刻(平滑化オフなら空)
    Xs::Array{Float64,3}                 # 平滑化アンサンブル(行 × スナップ × メンバー)
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

"""
    select_masked_rows(cfg, batch) -> Vector{Int}

現在時刻解析の変数局所化(DECISIONS #0040-(α))の対象行を選ぶ。
`cfg.analysis_masked_vars` が空なら常に `Int[]`(既定・従来動作)。
非空でも、`batch` に `cfg.analysis_unmask_names` のいずれかの観測名が
1つでも含まれればマスク解除(`Int[]`)。それ以外は `cfg.analysis_masked_vars`
をそのまま返す。
"""
function select_masked_rows(cfg::AssimConfig, batch::AbstractVector{ObservationRecord})
    isempty(cfg.analysis_masked_vars) && return Int[]
    any(o.spec.name in cfg.analysis_unmask_names for o in batch) && return Int[]
    return cfg.analysis_masked_vars
end

"theta_sig を差し替えた ModelParameters(状態拡大メンバー用)"
with_theta_sig(p::ModelParameters, theta::Real) =
    ModelParameters(p.regime, p.l1, p.l2, L3Params(theta_sig = float(theta)),
                    p.exo, p.x0_nat, p.x0)

"""
    run_assimilation(params, E0, obs, event_times; cfg, seed,
                     augmented=false, obs_counts=nothing, count_scale=1.0)
        -> AssimResult

初期アンサンブル `E0`(n × N。augmented なら n = 14 で最終行 = log theta_sig)
から §9 のハイブリッド同化(EnKF + ポアソン重み + イベント同期)を実行する。

週次イベントカウントの扱い(DECISIONS #0031):
- 既定(`obs_counts = nothing`): E1 と同じく `event_times`(真値カタログ)を
  窓に集計して観測カウントとする(モデルジャンプ = 観測イベントが1対1)。
- 実データ(M8): `obs_counts` に窓別の観測カウント列(窓 k は区間
  [t0+(k−1)·event_window, t0+k·event_window))を渡す。**負値はデータなし**を
  意味し、その窓のポアソン重み更新をスキップする(#0031-3)。
  `count_scale` は報告率 ν(N_w 〜 Poisson(ν·Λ)、#0031-1)。`count_temper` は
  過分散カウントの尤度テンパリング係数(1/ν 推奨、#0033。既定 1 = 素のポアソン)。
"""
function run_assimilation(params::ModelParameters, E0::Matrix{Float64},
                          obs::Vector{ObservationRecord},
                          event_times::Vector{Float64};
                          cfg::AssimConfig = AssimConfig(), seed::Integer,
                          augmented::Bool = false,
                          obs_counts::Union{Nothing, Vector{Int}} = nothing,
                          count_scale::Float64 = 1.0,
                          count_temper::Float64 = 1.0)
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
    # ランク計算専用の独立ストリーム(観測ノイズ抽選が力学の乱数列を乱さないように)
    rank_rng = Xoshiro(member_seed(seed, 10_000_019))
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

    # 固定ラグ EnKS(#0024): smoother_dt 刻みでスナップショットを保持し、
    # 解析のたびに現在時刻から smoother_lag 以内のものを同時更新する。
    smoothing = cfg.smoother_lag > 0
    snap_steps = max(1, round(Int, cfg.smoother_dt / cfg.dt))
    snap_ts = Float64[]
    snaps = Matrix{Float64}[]
    lag_start = 1                          # ラグ窓内の最初のスナップショット index
    if smoothing
        push!(snap_ts, ts[1])
        push!(snaps, copy(E))
    end

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
            guard_sigma_drift!(f)                    # σ_s ガード(#0032)
            diffusion!(sig, xi, p_i, t)
            randn!(rngs[i], dW)
            @. xi += cfg.dt * f + sqdt * sig * dW
            guard_sigma_state!(xi)
            if augmented   # d(param) = 0 + 微小ノイズ(§8.3)
                E[end, i] += cfg.param_noise_sd * sqdt * randn(rngs[i])
            end
        end

        # (c) 週次窓の終端: ポアソン重みを累積し、ESS < N/2 で系統再抽選(§9.3)
        if step % wsteps == 0
            # 観測カウント: 既定はカタログ集計(E1)、実データでは窓別列(#0031)。
            # 負値 = データなし窓 → 重み更新スキップ(病的ガードは常時)。
            widx = step ÷ wsteps
            observed = obs_counts === nothing ? window_count :
                       (widx <= length(obs_counts) ? obs_counts[widx] : -1)
            if observed >= 0
                # count_temper: 過分散カウントの情報量換算(#0033。既定1 = 素のポアソン)
                logw .+= count_temper .*
                         poisson_logweights(observed, count_scale .* Lambda)
            end
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
                # メンバー対応を保つため、ラグ窓内のスナップショットも同じ
                # インデックスで再抽選する(EnKS、#0024)
                for s in lag_start:length(snaps)
                    snaps[s] .= snaps[s][:, idx]
                end
                fill!(logw, 0.0)
                nresample += 1
            end
            fill!(Lambda, 0.0)
            window_count = 0
        end

        # (d) 解析ステップ(この時刻に届いた観測のみ、§9.2)
        if haskey(obs_at, step + 1)
            batch = obs_at[step + 1]
            # tauA への緩い擬似観測(DECISIONS #0036、既定オフ)。batch を
            # コピーして追加するため obs_at 由来の元配列は変更しない。
            if cfg.tauA_pseudo_sd_mult > 0
                batch = augment_tauA_pseudo(batch, cfg.tauA_pseudo_sd_mult)
            end
            # ランク(解析直前の事前アンサンブルに対する観測の順位)。
            # 観測 = 真値 + ノイズ のため、メンバー側にも観測ノイズ抽選を
            # 加えるのがランクヒストグラムの標準定義(Hamill 2001、#0017)。
            # 省くと高頻度観測変数で見かけの過小分散が生じる。
            for o in batch
                yj = [o.spec.h(view(E, 1:N_STATE, j)) + o.spec.sd * randn(rank_rng)
                      for j in 1:N]
                push!(get!(ranks, o.spec.name, Int[]),
                      count(<(o.value), yj) + 1)
            end
            yobs = [o.value for o in batch]
            R = Diagonal([o.spec.sd^2 for o in batch]) |> Matrix
            hfun = col -> [o.spec.h(view(col, 1:N_STATE)) for o in batch]

            # 現在時刻解析の変数局所化(#0040-(α))
            masked_rows = select_masked_rows(cfg, batch)

            # ラグ窓の前進(EnKS。窓外のスナップショットは確定)
            while smoothing && lag_start <= length(snap_ts) &&
                  snap_ts[lag_start] < t_next - cfg.smoother_lag
                lag_start += 1
            end
            # 平滑化は疎な観測を含む解析のみ(smoother_min_period)
            do_smooth = smoothing &&
                any(o.spec.period >= cfg.smoother_min_period for o in batch)
            window_snaps = do_smooth ? view(snaps, lag_start:length(snaps)) :
                           Matrix{Float64}[]

            # スプレッド注入(inflation_mode、DECISIONS #0013)。
            # 平滑化更新は現在状態と同一のイノベーションで行う(#0024)。
            if cfg.inflation_mode === :rtps
                spread_prior = ensemble_spread(E)
                enks_analysis!(E, window_snaps, yobs, hfun, R;
                               rng = rngs[1], rho_inf = 1.0,
                               smooth_rows = cfg.smoother_vars,
                               masked_rows)
                alpha = jump_since_analysis ? cfg.rtps_alpha_jump : cfg.rtps_alpha
                rtps!(E, spread_prior; alpha)
            else
                rho_base = jump_since_analysis ? cfg.rho_inf_jump : cfg.rho_inf
                rho = cfg.inflation_mode === :per_time ?
                    rho_base^((t_next - t_last_analysis) / cfg.tau_ref) : rho_base
                enks_analysis!(E, window_snaps, yobs, hfun, R;
                               rng = rngs[1], rho_inf = rho,
                               smooth_rows = cfg.smoother_vars,
                               masked_rows)
            end
            postprocess_analysis!(E)
            jump_since_analysis = false
            t_last_analysis = t_next
        end

        # EnKS スナップショットの追加(解析後の状態、snap_dt 刻み)
        if smoothing && step % snap_steps == 0
            push!(snap_ts, t_next)
            push!(snaps, copy(E))
        end

        X[:, step + 1, :] = E
    end

    if smoothing
        Xs = Array{Float64,3}(undef, n, length(snaps), N)
        for (s, S) in enumerate(snaps)
            Xs[:, s, :] = S
        end
        return AssimResult(ts, X, ranks, ess_hist, nresample, snap_ts, Xs)
    end
    return AssimResult(ts, X, ranks, ess_hist, nresample,
                       Float64[], Array{Float64,3}(undef, 0, 0, 0))
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
