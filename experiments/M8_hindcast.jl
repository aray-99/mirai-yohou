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
using MiraiYohou
using MiraiYohou: build_params, run_assimilation, free_ensemble, AssimConfig,
                  ObservationRecord, member_seed, N_STATE,
                  IX_P, IX_W, IX_H, IX_K, IX_G, IX_T, IX_PHI, IX_V,
                  IX_TAU, IX_TAUA, IX_SIG, IX_PP, IX_LAME, logit

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

"観測空間の被覆率(事前90%区間、Hamill 定義 #0017 と同じノイズ付加)"
function coverage(recs, ts, X, window; rng)
    n_in = 0; n_tot = 0
    for r in recs
        window[1] <= r.t < window[2] || continue
        k = clamp(searchsortedlast(ts, r.t), 1, length(ts))
        y = [r.spec.h(view(X, :, k, j)) + r.spec.sd * randn(rng)
             for j in 1:size(X, 3)]
        q05, q95 = quantile(y, 0.05), quantile(y, 0.95)
        n_in += (q05 <= r.value <= q95)
        n_tot += 1
    end
    return n_in / max(n_tot, 1), n_tot
end

"観測変数ごとの RMSE(アンサンブル平均 vs 観測値)"
function obs_rmse(recs, ts, X, window)
    err = Dict{Symbol, Vector{Float64}}()
    for r in recs
        window[1] <= r.t < window[2] || continue
        k = clamp(searchsortedlast(ts, r.t), 1, length(ts))
        m = mean(r.spec.h(view(X, :, k, j)) for j in 1:size(X, 3))
        push!(get!(err, r.spec.name, Float64[]), m - r.value)
    end
    return Dict(k => sqrt(mean(abs2, v)) for (k, v) in err)
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

function run_country(country::String; N::Int = 100, seed::Integer = 20260708,
                     smoke::Bool = false)
    ccfg = COUNTRY_CFG[country]
    params0 = build_params(ccfg.regime)
    println("== $country ($(ccfg.regime)) ==")
    recs = build_observations(country, params0; t1 = T1)
    println("  観測: $(length(recs)) 点")
    c_v0 = init_c_v0(recs, ccfg.calib)
    params = build_params(ccfg.regime; c_v0)
    println("  c_v0 初期値(モーメント法) = ", round(c_v0, digits = 2))

    cfg = AssimConfig(t0 = 0.0, t1 = smoke ? ccfg.calib[2] : T1,
                      smoother_lag = 5.0)
    E0 = initial_ensemble(country, params, recs; N, seed = seed + 1)
    obs_counts = build_obs_counts(country, cfg)
    ncov = count(>=(0), obs_counts)
    println("  週次窓: $(length(obs_counts))(データあり $ncov)")
    event_times = build_forced_jumps(country)
    filter!(t -> t < cfg.t1, event_times)

    nu = init_nu(obs_counts, cfg, params, ccfg.calib)   # 報告率 ν の初期値(#0031-1)
    println("  ν 初期値(モーメント法) = ", round(nu, digits = 1))
    recs_run = [r for r in recs if r.t <= cfg.t1]
    res = run_assimilation(params, E0, recs_run, event_times;
                           cfg, seed, obs_counts, count_scale = nu,
                           count_temper = 1 / nu)   # #0033
    println("  同化完了: 再抽選 $(res.nresample) 回, ESS範囲 ",
            round.(extrema(res.ess), digits = 1))
    Xf = free_ensemble(params, E0; cfg, seed = seed + 2)

    rng = Xoshiro(seed + 3)
    win = smoke ? ccfg.calib : ccfg.verif
    cov, ntot = coverage(recs_run, res.t, res.X, win; rng)
    println("  被覆率($(win)) = ", round(cov, digits = 3), " (n=$ntot)")
    ra = obs_rmse(recs_run, res.t, res.X, win)
    rf = obs_rmse(recs_run, res.t, Xf, win)
    for k in sort(collect(keys(ra)); by = string)
        println("  RMSE $(rpad(k, 9)) 同化 ", round(ra[k], digits = 4),
                "  自由 ", round(rf[k], digits = 4),
                "  比 ", round(ra[k] / rf[k], digits = 3))
    end
    return (; res, Xf, recs = recs_run, cfg)
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    smoke = "--smoke" in ARGS
    countries = [c for c in ARGS if haskey(COUNTRY_CFG, c)]
    isempty(countries) && (countries = ["JPN", "THA"])
    for c in countries
        run_country(c; smoke)
    end
end
