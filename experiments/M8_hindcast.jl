# M8 単国ヒンドキャスト ドライバ(DECISIONS #0030/#0031)
#
# 実データ(experiments/data/raw)による日本・タイのヒンドキャスト:
#   初期アンサンブル(1990 年観測 + 摂動)→ ハイブリッド同化
#   (EnKF + ポアソン重み ν·Λ + 強制ジャンプ)→ 自由ラン対照 → 検証指標。
#
# 時間座標: t = 0 が 1990-01-01(#0030-2)。較正/検証分割は #0030-4。
# 本スクリプトの数値(R・スプレッド・ν 初期値等)は較正段階の暫定値であり、
# 検証ラン前に DECISIONS で凍結する(#0030-7)。
#
# 実行: julia --project=experiments experiments/M8_hindcast.jl [JPN|THA] [--smoke]

using Dates
using Random
using Statistics
using TOML
using JSON3
using MiraiYohou
using MiraiYohou: build_params, run_assimilation, free_ensemble, AssimConfig,
                  ObservationRecord, member_seed, N_STATE, intensity, simulate_sde,
                  IX_P, IX_W, IX_H, IX_K, IX_G, IX_T, IX_PHI, IX_V,
                  IX_TAU, IX_TAUA, IX_SIG, IX_PP, IX_LAME, logit,
                  AugmentedParam, augment_ensemble

include(joinpath(@__DIR__, "data", "build_observations.jl"))
include(joinpath(@__DIR__, "data", "prepare_events.jl"))

const T_EPOCH = Date(1990, 1, 1)
date_to_t(d::Date) = Dates.value(d - T_EPOCH) / 365.25

"国別設定(#0030-3/-4)。時刻は 1990 年起点の年"
const COUNTRY_CFG = Dict(
    "JPN" => (regime = :stable, calib = (5.0, 26.0), verif = (26.0, 35.0),
              exclude_admin1 = String[], acled_from = date_to_t(Date(2018, 1, 1))),
    "THA" => (regime = :volatile, calib = (20.0, 28.0), verif = (28.0, 35.0),
              exclude_admin1 = DEEP_SOUTH_THA, acled_from = date_to_t(Date(2010, 1, 1))),
)

const T1 = 35.0                      # 2024 年末まで(ACLED エンバーゴ内)

"""
駆動パラメータの L3 状態拡大(DECISIONS #0046)の初期スプレッド・RW sd。
較正窓 smoke でのみ感度確認する暫定値であり、検証ラン前に凍結エントリで
確定する(#0030-7 の流儀)。**この1箇所を差し替えれば済む**ようにまとめる。
theta_sig は既存の状態拡大機構(#0010系。E1 と同型)を M8 にも適用したもの、
残り4個は #0046-1 の被覆崩壊変数対応(P/T_proxy/phi/v/THA g,y の駆動自由度)。
"""
const M8_AUG_SETTINGS = (
    theta_sig = (link = :log,      init_sd = 0.5,   rw_sd = 0.01),
    netgrowth = (link = :identity, init_sd = 0.003, rw_sd = 0.002),
    a_T       = (link = :log,      init_sd = 0.1,   rw_sd = 0.05),
    r_phi     = (link = :log,      init_sd = 0.1,   rw_sd = 0.05),
    mu_gbar   = (link = :identity, init_sd = 0.1,   rw_sd = 0.05),
    c_v0      = (link = :identity, init_sd = 0.1,   rw_sd = 0.05),
)

"""
M8 の拡大パラメータ記述子(DECISIONS #0046/#0049)を `params`(fit_exogenous・
較正値適用済み)と国コード `country` から構築する。初期値は theta_sig(非較正の
L3 既定値)を除き較正済み中心(mu_gbar/c_v0 は EKI 較正値、netgrowth は
fit_exogenous 値、a_T/r_phi は較正対象外なので §8.2 のレジーム既定値)を使う。

**国条件付き除外(#0049-2)**: `country == "JPN"` では theta_sig を拡大集合
から除外する(イベント情報のほぼ無い国で RW が無拘束に漂流し λ を過大化する
副作用 — #0048-3)。THA は #0046 どおり theta_sig を含める。除外時、theta_sig
は augmented_params に現れないため member_params が触れず、較正済み定数と
しての従来(非拡大)動作に自動的に戻る。
"""
function build_m8_augmented_params(params::ModelParameters, country::AbstractString)
    s = M8_AUG_SETTINGS
    aug = AugmentedParam[]
    if country != "JPN"
        push!(aug, AugmentedParam(name = :theta_sig, link = s.theta_sig.link,
                       init = params.l3.theta_sig,
                       init_sd = s.theta_sig.init_sd, rw_sd = s.theta_sig.rw_sd))
    end
    append!(aug, AugmentedParam[
        AugmentedParam(name = :netgrowth, link = s.netgrowth.link,
                       init = params.exo.netgrowth,
                       init_sd = s.netgrowth.init_sd, rw_sd = s.netgrowth.rw_sd),
        AugmentedParam(name = :a_T, link = s.a_T.link,
                       init = params.l2.a_T,
                       init_sd = s.a_T.init_sd, rw_sd = s.a_T.rw_sd),
        AugmentedParam(name = :r_phi, link = s.r_phi.link,
                       init = params.l2.r_phi,
                       init_sd = s.r_phi.init_sd, rw_sd = s.r_phi.rw_sd),
        AugmentedParam(name = :mu_gbar, link = s.mu_gbar.link,
                       init = params.l2.mu_gbar,
                       init_sd = s.mu_gbar.init_sd, rw_sd = s.mu_gbar.rw_sd),
        AugmentedParam(name = :c_v0, link = s.c_v0.link,
                       init = params.l2.c_v0,
                       init_sd = s.c_v0.init_sd, rw_sd = s.c_v0.rw_sd),
    ])
    return aug
end

"""
初期アンサンブル: 中心 = レジーム既定値(§8.4)を 1990 年観測で上書きし、
全座標に sd `spread` の独立ガウス摂動(E1 手順3と同型。暫定値)。
"""
function initial_ensemble(country, params, recs::Vector{ObservationRecord};
                          N::Int, spread::Float64 = 0.3, seed::Integer)
    center = copy(params.x0)
    first_obs = Dict{Symbol, Float64}()
    for r in recs
        r.t <= 1.0 || continue
        haskey(first_obs, r.spec.name) || (first_obs[r.spec.name] = r.value)
    end
    haskey(first_obs, :P) && (center[IX_P] = first_obs[:P])
    haskey(first_obs, :w) && (center[IX_W] = first_obs[:w])
    haskey(first_obs, :g_swiid) && (center[IX_G] = first_obs[:g_swiid])
    haskey(first_obs, :T_proxy) && (center[IX_T] = first_obs[:T_proxy])
    haskey(first_obs, :phi) && (center[IX_PHI] = first_obs[:phi])
    haskey(first_obs, :v) && (center[IX_V] = first_obs[:v])
    haskey(first_obs, :tau) && (center[IX_TAU] = first_obs[:tau];
                                center[IX_TAUA] = first_obs[:tau])
    haskey(first_obs, :p) && (center[IX_PP] = first_obs[:p])
    rng = Xoshiro(seed)
    E0 = center .+ spread .* randn(rng, N_STATE, N)
    E0[IX_LAME, :] .= 0.0                            # λ_e は非負・初期0
    return E0
end

"窓別観測カウント列(#0031)。負値 = ACLED カバレッジ外"
function build_obs_counts(country, cfg::AssimConfig)
    ccfg = COUNTRY_CFG[country]
    nwin = floor(Int, round((cfg.t1 - cfg.t0) / cfg.dt) /
                      max(1, round(Int, cfg.event_window / cfg.dt)))
    counts = fill(-1, nwin)
    ev = political_events(load_events(country); exclude_admin1 = ccfg.exclude_admin1)
    tev = sort([date_to_t(e.date) for e in ev])
    w = cfg.event_window
    for k in 1:nwin
        lo, hi = cfg.t0 + (k - 1) * w, cfg.t0 + k * w
        lo >= ccfg.acled_from || continue
        hi <= T1 || continue
        counts[k] = searchsortedfirst(tev, hi) - searchsortedfirst(tev, lo)
    end
    return counts
end

"強制ジャンプ時刻(#0030-3: 較正期間分位のカタログ、週央に配置)"
function build_forced_jumps(country)
    ccfg = COUNTRY_CFG[country]
    ev = political_events(load_events(country); exclude_admin1 = ccfg.exclude_admin1)
    isempty(ev) && return Float64[]
    ws, c, f = weekly_counts(ev)
    lo = max(T_EPOCH + Day(round(Int, ccfg.calib[1] * 365.25)), ws[1])
    hi = T_EPOCH + Day(round(Int, ccfg.calib[2] * 365.25))
    if lo > hi   # ACLED カバレッジが較正期間と重ならない(JPN)→ 強制ジャンプなし
        println("  強制ジャンプ: なし(較正期間にイベントデータなし)")
        return Float64[]
    end
    jw, _, thr = jump_catalog(ws, c, f; calib_from = lo, calib_to = hi)
    times = [date_to_t(w) + 0.01 for w in jw if date_to_t(w) < T1]
    println("  強制ジャンプ: $(length(times)) 週(閾値 $(round(thr, digits=1)))")
    return times
end

# 評価は「事前(解析前)予測」で行う(#0034)。X[k] は観測時刻の解析後
# 状態なので、1グリッド前(k−1、O(dt)=3.65日前)の状態を prior の近似に使う。
_prior_index(ts, t) = clamp(searchsortedlast(ts, t) - 1, 1, length(ts))

"観測空間の被覆率(事前90%区間、Hamill 定義 #0017 と同じノイズ付加)"
function coverage(recs, ts, X, window; rng)
    n_in = 0; n_tot = 0
    for r in recs
        window[1] <= r.t < window[2] || continue
        k = _prior_index(ts, r.t)
        y = [r.spec.h(view(X, :, k, j)) + r.spec.sd * randn(rng)
             for j in 1:size(X, 3)]
        q05, q95 = quantile(y, 0.05), quantile(y, 0.95)
        n_in += (q05 <= r.value <= q95)
        n_tot += 1
    end
    return n_in / max(n_tot, 1), n_tot
end

"""
変数別診断(DECISIONS #0040-(γ)。判定はプール被覆値のみで行い、本関数は
診断能の追加に限る): 観測名ごとの被覆内訳(hit/n。coverage() と同一の
90%区間判定を独立乱数で再計算)、事前アンサンブルのスプレッド(sd)平均、
イノベーション |obs − 事前平均| の平均。
"""
function variable_diagnostics(recs, ts, X, window; rng)
    hits = Dict{Symbol,Int}(); ns = Dict{Symbol,Int}()
    spreads = Dict{Symbol,Vector{Float64}}(); innovs = Dict{Symbol,Vector{Float64}}()
    for r in recs
        window[1] <= r.t < window[2] || continue
        k = _prior_index(ts, r.t)
        raw = [r.spec.h(view(X, :, k, j)) for j in 1:size(X, 3)]
        y = raw .+ r.spec.sd .* randn(rng, length(raw))
        q05, q95 = quantile(y, 0.05), quantile(y, 0.95)
        name = r.spec.name
        hits[name] = get(hits, name, 0) + (q05 <= r.value <= q95)
        ns[name] = get(ns, name, 0) + 1
        push!(get!(spreads, name, Float64[]), std(raw))
        push!(get!(innovs, name, Float64[]), abs(r.value - mean(raw)))
    end
    names = union(keys(ns), keys(spreads))
    return Dict(string(name) => Dict(
                    "hit" => get(hits, name, 0), "n" => get(ns, name, 0),
                    "prior_spread_mean" => mean(get(spreads, name, [NaN])),
                    "prior_innov_mean" => mean(get(innovs, name, [NaN])))
                for name in names)
end

"凍結 TOML の meta.decisions(来歴サイドカー用。ファイルが無ければ nothing)"
function frozen_decisions_string()
    path = joinpath(@__DIR__, "M8_frozen_config.toml")
    isfile(path) || return nothing
    meta = get(TOML.parsefile(path), "meta", Dict())
    return get(meta, "decisions", nothing)
end

"""
1年先予測 RMSE(#0034 基準2): 各観測時刻 s について、同化ランの s−1 時点の
アンサンブルから純予測(内生 Hawkes・同化なし)で1年間発射し、H(x(s)) の
アンサンブル平均と観測を比較する。
"""
function forecast_rmse(recs, res, params, window; horizon::Float64 = 1.0,
                       seed::Integer)
    err = Dict{Symbol, Vector{Float64}}()
    N = size(res.X, 3)
    vals = Vector{Float64}(undef, N)
    for (ridx, r) in enumerate(recs)
        window[1] <= r.t < window[2] || continue
        t0 = r.t - horizon
        t0 >= res.t[1] || continue
        k = clamp(searchsortedlast(res.t, t0), 1, length(res.t))
        Threads.@threads for j in 1:N
            sim = simulate_sde(params; seed = member_seed(seed + 1000 * ridx, j),
                               t0 = res.t[k], t1 = r.t + 0.011, dt = 0.01,
                               xi0 = collect(view(res.X, 1:N_STATE, k, j)))
            vals[j] = r.spec.h(view(sim.traj.X, :, size(sim.traj.X, 2)))
        end
        push!(get!(err, r.spec.name, Float64[]), mean(vals) - r.value)
    end
    return Dict(k => sqrt(mean(abs2, v)) for (k, v) in err)
end

"""
週次カウントの予測ポアソン対数尤度(#0034 基準3): Σ_w log Poisson(N_w | ν·Λ_w)。
Λ_w は状態軌道 X から窓ごとに積分(アンサンブル平均強度)。素の尤度
(テンパリングなし)で同化ラン・自由ランを同一比較する。
"""
function count_loglik(ts, X, obs_counts, params, nu, cfg, window)
    w = cfg.event_window
    wsteps = max(1, round(Int, w / cfg.dt))
    N = size(X, 3)
    ll = 0.0; nwin = 0
    for (k, c) in enumerate(obs_counts)
        c >= 0 || continue
        lo = cfg.t0 + (k - 1) * w
        window[1] <= lo < window[2] || continue
        lam = 0.0
        i0 = (k - 1) * wsteps + 1
        for i in i0:min(i0 + wsteps - 1, length(ts))
            lam += mean(intensity(view(X, 1:N_STATE, i, j), params)
                        for j in 1:N) * cfg.dt
        end
        e = max(nu * lam, 1e-10)
        ll += c * log(e) - e            # log N_w! は両ラン共通の定数なので省略
        nwin += 1
    end
    return ll, nwin
end

"観測変数ごとの RMSE(事前アンサンブル平均 vs 観測値)"
function obs_rmse(recs, ts, X, window)
    err = Dict{Symbol, Vector{Float64}}()
    for r in recs
        window[1] <= r.t < window[2] || continue
        k = _prior_index(ts, r.t)
        m = mean(r.spec.h(view(X, :, k, j)) for j in 1:size(X, 3))
        push!(get!(err, r.spec.name, Float64[]), m - r.value)
    end
    return Dict(k => sqrt(mean(abs2, v)) for (k, v) in err)
end

"""
外生入力を実データから推定(§5: b−d+mig と wbar は既知の時間関数。
ConstantExogenous の枠内で、較正期間の平均成長率・平均 w に固定する)。
netgrowth = 較正期間の平均 d(log P)/dt、wbar = 較正期間の平均 w。
"""
function fit_exogenous(params, recs, calib)
    zP = sort([(r.t, r.value) for r in recs if r.spec.name === :P &&
               calib[1] <= r.t < calib[2]])
    netgrowth = length(zP) >= 2 ?
        (zP[end][2] - zP[1][2]) / (zP[end][1] - zP[1][1]) : params.exo.netgrowth
    zw = [r.value for r in recs if r.spec.name === :w && calib[1] <= r.t < calib[2]]
    wbar = isempty(zw) ? params.exo.wbar : 1 / (1 + exp(-mean(zw)))
    exo = MiraiYohou.ConstantExogenous(; netgrowth, wbar)
    return MiraiYohou.ModelParameters(params.regime, params.l1, params.l2,
                                      params.l3, exo, params.x0_nat, params.x0)
end

"""
c_v0 のモーメント初期化(#0031-2): 式(8)の平衡 xi_v ≈ xi_T + c_v0 から、
較正期間の観測平均で c_v0 ≈ mean(z_v) − mean(z_T)。較正の初期値であり
凍結値ではない。
"""
function init_c_v0(recs, calib)
    zv = [r.value for r in recs if r.spec.name === :v && calib[1] <= r.t < calib[2]]
    zT = [r.value for r in recs if r.spec.name === :T_proxy && calib[1] <= r.t < calib[2]]
    (isempty(zv) || isempty(zT)) && return 0.0
    return mean(zv) - mean(zT)
end

"ν のモーメント初期化(#0031-1): 較正期間の平均週次カウント / (lam0×窓幅)"
function init_nu(obs_counts, cfg, params, calib)
    w = cfg.event_window
    tot = 0; nwin = 0
    for (k, c) in enumerate(obs_counts)
        c >= 0 || continue
        lo = cfg.t0 + (k - 1) * w
        calib[1] <= lo < calib[2] || continue
        tot += c; nwin += 1
    end
    nwin == 0 && return 1.0
    return max(tot / (nwin * w * params.l2.lam0), 1.0)
end

"""
較正結果を読む。凍結 TOML(M8_frozen_config.toml、#0034)があればそれを優先し、
なければ EKI 出力 JSON(較正作業中の未凍結値)。nothing = 未較正の初期値運用。
"""
function load_calibrated(country)
    frozen = joinpath(@__DIR__, "M8_frozen_config.toml")
    if isfile(frozen)
        cfg = TOML.parsefile(frozen)
        haskey(cfg, country) || return nothing
        sec = cfg[country]
        return (; (Symbol(k) => Float64(v) for (k, v) in sec
                   if v isa Real && k != "calib_window")...)
    end
    path = joinpath(@__DIR__, "output", "M8_calib_$(country).json")
    isfile(path) || return nothing
    return JSON3.read(read(path, String)).params
end

function run_country(country::String; N::Int = 100, seed::Integer = 20260708,
                     smoke::Bool = false, calibrated::Bool = false)
    ccfg = COUNTRY_CFG[country]
    params0 = build_params(ccfg.regime)
    println("== $country ($(ccfg.regime)) ==")
    recs = build_observations(country, params0; t1 = T1)
    println("  観測: $(length(recs)) 点")
    cal = calibrated ? load_calibrated(country) : nothing
    nu_cal = nothing
    if cal !== nothing
        kw = Dict(Symbol(k) => Float64(v) for (k, v) in pairs(cal) if k !== :nu)
        params = fit_exogenous(build_params(ccfg.regime; kw...), recs, ccfg.calib)
        nu_cal = max(Float64(cal.nu), 1.0)
        println("  較正済み θ̂ を使用: ", kw, "  ν = ", round(nu_cal, digits = 1))
    else
        c_v0 = init_c_v0(recs, ccfg.calib)
        params = fit_exogenous(build_params(ccfg.regime; c_v0), recs, ccfg.calib)
        println("  c_v0 初期値(モーメント法) = ", round(c_v0, digits = 2))
    end
    println("  外生: netgrowth = ", round(params.exo.netgrowth, digits = 4),
            " wbar = ", round(params.exo.wbar, digits = 3))

    cfg = AssimConfig(t0 = 0.0, t1 = smoke ? ccfg.calib[2] : T1,
                      smoother_lag = 5.0,
                      smoother_vars = [IX_G, IX_TAU, IX_SIG, IX_PP],  # #0036: tauA除外
                      tauA_pseudo_sd_mult = 3.0,                     # #0036
                      analysis_masked_vars = [IX_TAUA],               # #0040-(α)
                      analysis_unmask_names = [:tau, :tauA_pseudo],   # #0040-(α)
                      rtps_alpha = 0.85,                              # #0040-(β)
                      obs_spread_floor_frac = 0.5)                    # #0043
    E0_state = initial_ensemble(country, params, recs; N, seed = seed + 1)
    aug = build_m8_augmented_params(params, country)        # #0046/#0049
    E0 = augment_ensemble(E0_state, aug; rng = Xoshiro(seed + 6))
    obs_counts = build_obs_counts(country, cfg)
    ncov = count(>=(0), obs_counts)
    println("  週次窓: $(length(obs_counts))(データあり $ncov)")
    event_times = build_forced_jumps(country)
    filter!(t -> t < cfg.t1, event_times)

    nu = nu_cal !== nothing ? nu_cal :
         init_nu(obs_counts, cfg, params, ccfg.calib)   # 報告率 ν(#0031-1)
    println("  ν = ", round(nu, digits = 1),
            nu_cal === nothing ? "(モーメント初期値)" : "(較正プロファイル値)")
    recs_run = [r for r in recs if r.t <= cfg.t1]
    res = run_assimilation(params, E0, recs_run, event_times;
                           cfg, seed, obs_counts, count_scale = nu,
                           count_temper = 1 / nu,           # #0033
                           augmented_params = aug)          # #0046
    println("  同化完了: 再抽選 $(res.nresample) 回, ESS範囲 ",
            round.(extrema(res.ess), digits = 1))
    for (k, ap) in enumerate(aug)
        final = res.X[N_STATE + k, end, :]
        natural = ap.link === :log ? exp.(final) : final
        println("  拡大パラメータ $(ap.name): 事後平均 ", round(mean(natural), digits = 5),
                "(初期 ", round(ap.init, digits = 5), ")")
    end
    Xf = free_ensemble(params, E0; cfg, seed = seed + 2)

    rng = Xoshiro(seed + 3)
    win = smoke ? ccfg.calib : ccfg.verif
    cov, ntot = coverage(recs_run, res.t, res.X, win; rng)
    println("  [基準1] 被覆率($(win)) = ", round(cov, digits = 3), " (n=$ntot)",
            0.80 <= cov <= 0.98 ? "  PASS" : "  FAIL")
    ra = smoke ? obs_rmse(recs_run, res.t, res.X, win) :
         forecast_rmse(recs_run, res, params, win; seed = seed + 4)
    rf = obs_rmse(recs_run, res.t, Xf, win)
    nbetter = 0; nvars = 0
    for k in sort(collect(keys(ra)); by = string)
        haskey(rf, k) || continue
        nvars += 1; nbetter += ra[k] < rf[k]
        println("  RMSE$(smoke ? "" : "(1年先)") $(rpad(k, 9)) 同化 ",
                round(ra[k], digits = 4), "  自由 ", round(rf[k], digits = 4),
                "  比 ", round(ra[k] / rf[k], digits = 3))
    end
    println("  [基準2] 改善 $nbetter/$nvars",
            nbetter > nvars ÷ 2 ? "  PASS" : "  FAIL")
    lla, nwin_ll = count_loglik(res.t, res.X, obs_counts, params, nu, cfg, win)
    if nwin_ll > 0
        llf, _ = count_loglik(res.t, Xf, obs_counts, params, nu, cfg, win)
        println("  [基準3] カウント予測 logL 同化 ", round(lla, digits = 1),
                "  自由 ", round(llf, digits = 1),
                lla > llf ? "  PASS" : "  FAIL")
    end

    # 変数別診断の成果物化(DECISIONS #0040-(γ)。判定・表示はプール被覆値の
    # ままで変更しない。検証ラン(--calib かつ非 smoke)のみ保存し、smoke
    # ランでは保存をスキップする)。
    if !smoke && calibrated
        diag_vars = variable_diagnostics(recs_run, res.t, res.X, win;
                                         rng = Xoshiro(seed + 5))
        diag_out = Dict(
            "country" => country,
            "window" => collect(win),
            "variables" => diag_vars,
            "provenance" => Dict(
                "commit" => strip(read(`git -C $(dirname(@__DIR__)) rev-parse HEAD`,
                                        String)),
                "seed" => seed,
                "generated_at" => Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"),
                "frozen_decisions" => frozen_decisions_string()))
        diag_path = joinpath(@__DIR__, "output", "M8_verify_diag_$(country).json")
        mkpath(dirname(diag_path))
        write(diag_path, JSON3.write(diag_out))
        println("  変数別診断: $diag_path")
    end

    return (; res, Xf, recs = recs_run, cfg)
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    smoke = "--smoke" in ARGS
    calibrated = "--calib" in ARGS
    countries = [c for c in ARGS if haskey(COUNTRY_CFG, c)]
    isempty(countries) && (countries = ["JPN", "THA"])
    for c in countries
        run_country(c; smoke, calibrated)
    end
end
