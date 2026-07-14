# M11 30年ホライズン数値安定性スモーク(Issue #2、docs/PHASE3_DESIGN.md §2)
#
# 将来予報(M11)のホライズンは 30 年だが、現行の予報アンサンブルは最長 6 年
# でしか回していない。実装前に、長期積分での数値安定性(特に σ_s の log 座標
# 剛性に対する tame+floor ガード #0032、src/integrator.jl)を確認する。
#
# 手順(国ごと):
#   1. M8 凍結較正値(M8_frozen_config.toml)で params を構築(fit_exogenous は
#      較正窓、run_country --calib と同じ流儀)
#   2. 較正窓 [win_start, t_k](t_k = 較正窓末尾 = M9/M10 の初回オリジン)まで
#      M10 と同一の AssimConfig・NegBin 手続きで同化し、窓末尾状態を得る
#   3. その末尾アンサンブル状態から forecast_ensemble(M10_walkforward.jl 経由、
#      内生ジャンプ + 拡大パラメータ RW 継続 #0053、Γ シンニング p_ex #0068 は
#      walk-forward と同じ規則)で 30 年 × N=100(M10 と同サイズ)を積分
#   4. 全状態の NaN/Inf 検査、tame/floor ガード発火の集計、変数別 5-95% 分位幅
#      時系列の保存(experiments/output/M11_horizon_smoke_{ISO3}.json +
#      .meta.json 来歴サイドカー)
#
# ガード発火の検出は保存軌道からの事後再構成(コア変更なしの読み取り専用):
#   - floor: 積分器は EM ステップ後に xi_sig = max(xi_sig, XI_SIG_FLOOR) を適用
#     する。拡散ノイズが加わった直後の値が床と厳密一致する確率は 0 なので、
#     保存グリッド値 == XI_SIG_FLOOR ⟺ そのステップで floor が発火。厳密。
#   - tame: 各保存グリッド状態(拡大行から build_member_params でメンバー
#     パラメータを再構築)で drift を再評価し |f_sig| > SIG_DRIFT_MAX を数える。
#     積分器は「区間内ジャンプ適用後」の状態で drift を評価するため、ジャンプ
#     発生ステップでは保存状態(区間開始時点)との差が O(1ジャンプ) だけ
#     生じ得る近似(ジャンプは希少なのでスモーク診断としては十分)。
#
# 凍結値・合格基準の変更なし。シード 20260711(#0063)を来歴に記録。
#
# 実行: julia --project=experiments -t 8 experiments/M11_horizon_smoke.jl JPN THA
#         [--horizon 30] [--N 100]

include(joinpath(@__DIR__, "M10_walkforward.jl"))   # M9/M8_calibrate/M8_hindcast も連鎖 include 済み

using MiraiYohou: drift!, build_member_params
const SIG_DRIFT_MAX_SMOKE = MiraiYohou.SIG_DRIFT_MAX   # 50.0(#0032)
const XI_SIG_FLOOR_SMOKE = MiraiYohou.XI_SIG_FLOOR     # -12.0(#0032)

"""
    guard_stats(fe, params, aug) -> Dict

予報アンサンブル `fe`(forecast_ensemble の戻り値)の保存軌道から tame/floor
ガード(#0032)の発火を事後再構成で集計する(冒頭コメント参照。floor は厳密、
tame はジャンプ発生ステップで O(1ジャンプ) の近似)。発火時刻分布は予報開始
からの経過年(床関数)ごとのヒストグラムで返す。読み取り専用・乱数消費なし。
"""
function guard_stats(fe, params, aug::Vector{AugmentedParam})
    nsteps = length(fe.t) - 1
    N = size(fe.X, 3)
    t0 = fe.t[1]
    horizon_years = ceil(Int, fe.t[end] - t0)
    tame_hist = zeros(Int, horizon_years)
    floor_hist = zeros(Int, horizon_years)
    tame_total = 0; floor_total = 0
    tame_members = falses(N); floor_members = falses(N)
    f = Vector{Float64}(undef, N_STATE)
    for j in 1:N
        Ej = view(fe.X, :, :, j)
        for i in 1:nsteps
            yr = clamp(floor(Int, fe.t[i] - t0) + 1, 1, horizon_years)
            # floor 発火: EM 後の保存値が床と厳密一致(i+1 列 = ステップ i の結果)
            if fe.X[IX_SIG, i + 1, j] == XI_SIG_FLOOR_SMOKE
                floor_total += 1; floor_hist[yr] += 1; floor_members[j] = true
            end
            # tame 発火: 保存状態で drift を再評価し |f_sig| > 上限
            p = build_member_params(params, aug, view(Ej, :, i:i), N_STATE, 1)
            drift!(f, view(Ej, 1:N_STATE, i), p, fe.t[i])
            if abs(f[IX_SIG]) > SIG_DRIFT_MAX_SMOKE
                tame_total += 1; tame_hist[yr] += 1; tame_members[j] = true
            end
        end
    end
    total_steps = nsteps * N
    return Dict(
        "detection_note" => "post-hoc from saved trajectories: floor exact (grid value == XI_SIG_FLOOR after clamp), tame approximate at jump steps (drift re-evaluated at interval-start state; integrator evaluates after intra-step jumps)",
        "sig_drift_max" => SIG_DRIFT_MAX_SMOKE,
        "xi_sig_floor" => XI_SIG_FLOOR_SMOKE,
        "member_steps_total" => total_steps,
        "tame" => Dict("count" => tame_total,
                       "rate" => tame_total / total_steps,
                       "members_affected" => count(tame_members),
                       "hist_by_forecast_year" => tame_hist),
        "floor" => Dict("count" => floor_total,
                        "rate" => floor_total / total_steps,
                        "members_affected" => count(floor_members),
                        "hist_by_forecast_year" => floor_hist),
    )
end

"""
    quantile_width_series(fe, aug; stride = 25) -> Dict

予報アンサンブルの各行(状態 13 変数 + 拡大パラメータ行、保持座標のまま)に
ついて 5-95% 分位幅の時系列を `stride` グリッド刻み(dt=0.01 で既定 0.25 年)
で集計する。予報円の発散が単調か飽和かの定性確認用(Issue #2)。
"""
function quantile_width_series(fe, aug::Vector{AugmentedParam}; stride::Int = 25)
    idx = 1:stride:length(fe.t)
    names = vcat([string(s) for s in STATE_NAMES],
                 ["aug_" * string(ap.name) for ap in aug])
    series = Dict{String,Vector{Float64}}()
    for (r, name) in enumerate(names)
        series[name] = [quantile(view(fe.X, r, i, :), 0.95) -
                        quantile(view(fe.X, r, i, :), 0.05) for i in idx]
    end
    return Dict("t" => fe.t[idx], "coordinates_note" =>
                    "retained (internal) coordinates: log/logit per src/coordinates.jl; augmented rows in link coordinates",
                "width_5_95" => series)
end

"""
    run_horizon_smoke(country; N = 100, seed = 20260711, horizon = 30.0) -> Dict

冒頭コメントの手順 1〜4 を1国分実行する。同化部は M10 の `run_origin`
(M9_walkforward.jl)の (b) と同一の構成(AssimConfig・NegBin プロファイル・
若返り a=0.95・シードオフセット)だが、EKI 再較正 (a) は行わず M8 凍結値を
そのまま使う(スモークの目的は数値安定性の確認であり較正ではない)。
"""
function run_horizon_smoke(country::String; N::Int = 100, seed::Integer = 20260711,
                           horizon::Float64 = 30.0)
    ccfg = COUNTRY_CFG[country]
    win_start = ccfg.calib[1]
    t_k = ccfg.calib[2]          # 較正窓末尾 = M9/M10 の初回オリジン(JPN 26 / THA 28)
    window = (win_start, t_k)
    seed_o = seed + Int(t_k)     # run_origin 呼び出し規約(seed + t_k)と同じ
    println("== $country M11 horizon smoke: 較正窓 $window → 予報 $horizon 年 (N=$N, seed=$seed) ==")

    frozen = load_calibrated(country)
    frozen === nothing &&
        error("run_horizon_smoke: $country の M8 凍結値がありません(M8_frozen_config.toml)")
    frozen_center = Dict(CAL_PARAMS[k].name => Float64(getproperty(frozen, CAL_PARAMS[k].name))
                         for k in eachindex(CAL_PARAMS))

    # theta_sig 規則(#0052/#0054): 較正窓のフィルタ後 ΣN しきい値判定。
    # t_k = 初回オリジンなので M9/M10 の判定と同一値になる。
    cfg0 = AssimConfig(t0 = 0.0, t1 = Float64(t_k))
    total_counts = windowed_count_total(build_obs_counts(country, cfg0),
                                        cfg0.t0, cfg0.event_window, window)
    include_theta_sig = total_counts >= THETA_SIG_COUNT_MIN
    println("  theta_sig 規則: ΣN = $total_counts (しきい値 $THETA_SIG_COUNT_MIN) → ",
            include_theta_sig ? "拡大集合に含める" : "除外")

    # params 構築(run_country --calib / run_origin と同じ流儀)
    params0 = build_params(ccfg.regime)
    recs_all = build_observations(country, params0; t1 = Float64(t_k))
    params = fit_exogenous(build_params(ccfg.regime; frozen_center...), recs_all, window)
    recs_calib = [r for r in recs_all if r.t <= t_k]
    nu_frozen = max(Float64(frozen.nu), 1.0)

    # 較正窓末尾までの同化(run_origin (b) と同一の AssimConfig・手続き)
    cfg = AssimConfig(t0 = 0.0, t1 = Float64(t_k), smoother_lag = 5.0,
                      smoother_vars = [IX_G, IX_TAU, IX_SIG, IX_PP],
                      tauA_pseudo_sd_mult = 3.0,
                      analysis_masked_vars = [IX_TAUA],
                      analysis_unmask_names = [:tau, :tauA_pseudo],
                      rtps_alpha = 0.85,
                      obs_spread_floor_frac = 0.5,
                      rejuvenation_a = REJUVENATION_A)
    E0_state = initial_ensemble(country, params, recs_all; N, seed = seed_o + 1)
    aug = build_m8_augmented_params(params, country; include_theta_sig)
    E0 = augment_ensemble(E0_state, aug; rng = Xoshiro(seed_o + 6))
    obs_counts = build_obs_counts(country, cfg)
    event_times = filter(t -> t < cfg.t1,
                         build_forced_jumps(country; calib_window = window))
    res = run_assimilation(params, E0, recs_calib, event_times;
                           cfg, seed = seed_o, obs_counts, count_scale = nu_frozen,
                           count_temper = 1 / nu_frozen, augmented_params = aug)
    ks = count_windows_in(obs_counts, cfg, window)
    if isempty(ks)
        nu_star, r_hat = nu_frozen, Inf
        println("  カウント窓なし — (ν*, r̂) プロファイルをスキップ(実質ポアソン)")
    else
        lams = window_lambdas(res.t, res.X, params, cfg, ks)
        prof = profile_count_dispersion([obs_counts[k] for k in ks], lams)
        nu_star, r_hat = prof.nu_star, prof.r_hat
        println("  プロファイル ν* = ", round(nu_star, digits = 3),
                "  r̂ = ", round(r_hat, digits = 4), "(カウント窓 ", length(ks), ")")
        res = run_assimilation(params, E0, recs_calib, event_times;
                               cfg, seed = seed_o, obs_counts, count_scale = nu_star,
                               count_model = :negbin, count_dispersion = r_hat,
                               augmented_params = aug)
    end
    println("  同化完了(t=0→$t_k): 再抽選 $(res.nresample) 回, ESS範囲 ",
            round.(extrema(res.ess), digits = 1))

    # p_ex(#0068): walk-forward と同じ規則(強制ジャンプ週数/カウントデータ週数)
    n_forced_window = count(t -> window[1] <= t < window[2], event_times)
    p_ex = isempty(ks) ? 1.0 : n_forced_window / length(ks)
    println("  p_ex(Γ シンニング, #0068) = ", round(p_ex, digits = 4))

    # 30 年予報アンサンブル(内生ジャンプ + 拡大 RW 継続 #0053)
    t_fore = time()
    fe = forecast_ensemble(params, aug, res; horizon, seed = seed_o + 7,
                           gamma_thinning_p = p_ex)
    elapsed_fore = time() - t_fore
    println("  予報積分完了: $(round(elapsed_fore, digits = 1)) 秒 ",
            "($(size(fe.X, 1)) 行 × $(length(fe.t)) 時刻 × $N メンバー)")

    # (1) NaN/Inf 検査(拡大行を含む全状態)
    nonfinite = count(!isfinite, fe.X)
    finite_ok = nonfinite == 0
    println("  NaN/Inf 検査: ", finite_ok ? "PASS(全値有限)" : "FAIL($nonfinite 個)")

    # (2) ガード発火集計
    gstats = guard_stats(fe, params, aug)
    println("  ガード発火: tame ", gstats["tame"]["count"], " 回(",
            gstats["tame"]["members_affected"], "/$N メンバー) floor ",
            gstats["floor"]["count"], " 回(",
            gstats["floor"]["members_affected"], "/$N メンバー) / 総ステップ ",
            gstats["member_steps_total"])

    # (3) 変数別 5-95% 分位幅時系列(0.25 年刻み)
    widths = quantile_width_series(fe, aug)

    njumps = sum(length, fe.jumps)
    println("  内生ジャンプ: 合計 $njumps 発(全メンバー・$horizon 年)")

    return Dict(
        "country" => country,
        "regime" => string(ccfg.regime),
        "calib_window" => collect(window),
        "forecast_start" => Float64(t_k),
        "horizon" => horizon,
        "N" => N,
        "include_theta_sig" => include_theta_sig,
        "nu_star" => nu_star,
        "r_hat" => isfinite(r_hat) ? r_hat : nothing,
        "gamma_thinning_p" => p_ex,
        "assim_nresample" => res.nresample,
        "finite_check" => Dict("pass" => finite_ok, "nonfinite_values" => nonfinite),
        "guard_stats" => gstats,
        "quantile_widths" => widths,
        "endogenous_jumps_total" => njumps,
        "forecast_elapsed_sec" => elapsed_fore,
    )
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    horizon = 30.0
    N = 100
    for (i, a) in enumerate(ARGS)
        a == "--horizon" && (global horizon = parse(Float64, ARGS[i + 1]))
        a == "--N" && (global N = parse(Int, ARGS[i + 1]))
    end
    seed = 20260711   # #0063 凍結シード
    countries = [c for c in ARGS if haskey(COUNTRY_CFG, c)]
    isempty(countries) && (countries = ["JPN", "THA"])
    for c in countries
        t0 = time()
        out = run_horizon_smoke(c; N, seed, horizon)
        out["elapsed_sec"] = time() - t0
        provenance = Dict(
            "commit" => strip(read(`git -C $(dirname(@__DIR__)) rev-parse HEAD`, String)),
            "seed" => seed,
            "generated_at" => Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"),
            "frozen_decisions" => frozen_decisions_string(),
            "design_decision" => "Issue #2 (M11 horizon smoke); guards #0032; forecast rules #0053/#0068; frozen seed #0063",
            "script" => "experiments/M11_horizon_smoke.jl",
        )
        out["provenance"] = provenance
        path = joinpath(@__DIR__, "output", "M11_horizon_smoke_$(c).json")
        mkpath(dirname(path))
        write(path, JSON3.write(out))
        write(path * ".meta.json", JSON3.write(provenance))   # 来歴サイドカー
        println("保存: $path(+ .meta.json)")
    end
end
