# M10 mu_gbar アンカリング prior sd 事前選定(DECISIONS #0062 決定3)
#
# THA 初回オリジン(t=28)の較正窓(t=20..28、検証窓未参照)のみで、mu_gbar
# アンカリング prior の sd 候補 {0.1, 0.2, 0.3} それぞれについて EKI 較正 →
# 同化ラン(M10_walkforward.jl の run_origin と同一の較正・同化手続き、
# #0054 NegBin プロファイル込み)を実行し、窓内 g_swiid 観測に対する
#   (a) 事前平均バイアス(モデル予測 logit − 観測 logit の平均)
#   (b) 事前90%被覆(Hamill 定義、lag=1 = 解析更新前グリッド)
# を算出して JSON に保存する。**選定判定(どの sd を採用するか)は行わない**
# (#0062: 数値の算出まで。判定はオーケストレーター側)。
#
# 実行: julia --project=experiments -t 8 experiments/M10_prior_select.jl
#         [--smoke] [--candidates 0.1,0.2,0.3]
#   --smoke: N=40, J=12, iters=2(動作確認のみ)。
#   非 smoke 実行(本番の事前選定)は長時間ランのため、本スクリプト単体では
#   実行しない(オーケストレーターが切り離しプロトコルで別途実行する)。

include(joinpath(@__DIR__, "M10_walkforward.jl"))

const PRIOR_SELECT_COUNTRY = "THA"
const PRIOR_SELECT_ORIGIN = 28.0

"""
    calibrate_assimilate_window(country, t_k, prior_center; N, seed, J, iters,
                                N_eki, include_theta_sig, prior_sd_override)
        -> (; res, params, recs_all, window, nu_star, r_hat, theta_hat)

`M10_walkforward.jl` の `run_origin`(継承元 `M9_walkforward.jl`)の
(a) EKI 較正 + (b) t_k までの同化(#0054 NegBin プロファイル込み)と同一の
手続きを較正窓のみに限定して実行する。予報・自由ラン・評価は行わない
(事前選定は較正窓内の事前診断のみ必要、#0062 決定3)。
"""
function calibrate_assimilate_window(country::String, t_k::Real,
                                     prior_center::Dict;
                                     N::Int, seed::Integer, J::Int, iters::Int,
                                     N_eki::Int, include_theta_sig::Bool,
                                     prior_sd_override)
    ccfg = COUNTRY_CFG[country]
    win_start = ccfg.calib[1]
    window = (win_start, Float64(t_k))

    calib = calibrate(country; J, iters, N = N_eki, seed = seed + 101,
                      window, prior_center, prior_sd = 0.5, prior_sd_override,
                      save = false, include_theta_sig)
    theta_hat, nu_eki = calib.theta_hat, calib.nu_star
    kw = Dict(CAL_PARAMS[k].name => theta_hat[k] for k in eachindex(CAL_PARAMS))

    params0 = build_params(ccfg.regime)
    recs_all = build_observations(country, params0; t1 = Float64(t_k))
    params = fit_exogenous(build_params(ccfg.regime; kw...), recs_all, window)
    recs_calib = [r for r in recs_all if r.t <= t_k]

    cfg = AssimConfig(t0 = 0.0, t1 = Float64(t_k), smoother_lag = 5.0,
                      smoother_vars = [IX_G, IX_TAU, IX_SIG, IX_PP],
                      tauA_pseudo_sd_mult = 3.0,
                      analysis_masked_vars = [IX_TAUA],
                      analysis_unmask_names = [:tau, :tauA_pseudo],
                      rtps_alpha = 0.85,
                      obs_spread_floor_frac = 0.5,
                      rejuvenation_a = REJUVENATION_A)
    E0_state = initial_ensemble(country, params, recs_all; N, seed = seed + 1)
    aug = build_m8_augmented_params(params, country; include_theta_sig)
    E0 = augment_ensemble(E0_state, aug; rng = Xoshiro(seed + 6))
    obs_counts = build_obs_counts(country, cfg)
    event_times = filter(t -> t < cfg.t1,
                         build_forced_jumps(country; calib_window = window))
    # 予備ラン(現行ポアソン + 1/ν): (ν*, r̂) プロファイルの Λ̄_w 供給源
    res = run_assimilation(params, E0, recs_calib, event_times;
                           cfg, seed, obs_counts, count_scale = nu_eki,
                           count_temper = 1 / nu_eki, augmented_params = aug)
    ks = count_windows_in(obs_counts, cfg, window)
    nu_star, r_hat = nu_eki, Inf
    if !isempty(ks)
        lams = window_lambdas(res.t, res.X, params, cfg, ks)
        prof = profile_count_dispersion([obs_counts[k] for k in ks], lams)
        nu_star, r_hat = prof.nu_star, prof.r_hat
        # 本同化ラン: NegBin 尤度(#0054。テンパリングなし、同一シード)
        res = run_assimilation(params, E0, recs_calib, event_times;
                               cfg, seed, obs_counts, count_scale = nu_star,
                               count_model = :negbin, count_dispersion = r_hat,
                               augmented_params = aug)
    end
    return (; res, params, recs_all, window, nu_star, r_hat, theta_hat)
end

"""
    g_swiid_prior_diagnostics(recs, ts, X, window; rng, lag = 1) -> (; bias, coverage, n)

窓内 g_swiid 観測に対する事前(解析更新前グリッド、`_prior_index` 既定
lag=1、M8 `variable_diagnostics` と同じ規約)アンサンブル平均・90%区間から
(a) 事前平均バイアス(モデル予測 logit − 観測 logit の平均。符号: 正なら
モデルが観測より高い/過大評価)、(b) 事前90%被覆(Hamill)を算出する。
g_swiid は恒等観測(target = IX_G)なので h(xi) = xi[IX_G] がそのまま
logit 座標(#0062 のアンカリング対象と同じ座標系)。
"""
function g_swiid_prior_diagnostics(recs, ts, X, window; rng, lag::Int = 1)
    biases = Float64[]
    n_in = 0; n_tot = 0
    for r in recs
        r.spec.name === :g_swiid || continue
        window[1] <= r.t <= window[2] || continue
        k = _prior_index(ts, r.t; lag)
        raw = [r.spec.h(view(X, :, k, j)) for j in 1:size(X, 3)]
        y = raw .+ r.spec.sd .* randn(rng, length(raw))
        q05, q95 = quantile(y, 0.05), quantile(y, 0.95)
        push!(biases, mean(raw) - r.value)
        n_in += (q05 <= r.value <= q95)
        n_tot += 1
    end
    return (; bias = isempty(biases) ? NaN : mean(biases),
            coverage = n_tot > 0 ? n_in / n_tot : NaN,
            n = n_tot)
end

"""
    run_prior_select(; candidates, smoke, N, seed, J, iters, N_eki) -> Dict

候補 sd それぞれで THA 初回オリジンの較正窓のみ EKI 較正 → 同化し、窓内
g_swiid の事前バイアス・事前90%被覆を算出する(#0062 決定3)。
"""
function run_prior_select(; candidates::Vector{Float64} = [0.1, 0.2, 0.3],
                          smoke::Bool = false, N::Int = 100,
                          seed::Integer = 20260713, J::Int = 24, iters::Int = 4,
                          N_eki::Int = 100)
    if smoke
        N = 40; J = 12; iters = 2; N_eki = 40
    end
    country = PRIOR_SELECT_COUNTRY
    t_k = PRIOR_SELECT_ORIGIN
    ccfg = COUNTRY_CFG[country]
    win_start = ccfg.calib[1]
    window = (win_start, t_k)

    cfg0 = AssimConfig(t0 = 0.0, t1 = t_k)
    total_counts = windowed_count_total(build_obs_counts(country, cfg0),
                                        cfg0.t0, cfg0.event_window, window)
    include_theta_sig = total_counts >= THETA_SIG_COUNT_MIN

    frozen = load_calibrated(country)
    frozen === nothing && error("run_prior_select: $country の M8 凍結値がありません")
    frozen_center = Dict(CAL_PARAMS[k].name => Float64(getproperty(frozen, CAL_PARAMS[k].name))
                         for k in eachindex(CAL_PARAMS))
    anchor = gbar_anchor(country, window)
    prior_center = merge(frozen_center, Dict(:mu_gbar => anchor))

    println("== M10 prior_select: $country t_k=$t_k window=$window ==")
    println("  mu_gbar prior 中心(アンカリング) = ", round(anchor, digits = 3))

    results = Dict{String,Any}[]
    t_start = time()
    for sd in candidates
        println("-- 候補 mu_gbar_sd = $sd --")
        ca = calibrate_assimilate_window(country, t_k, prior_center;
                                         N, seed = seed + round(Int, sd * 1000),
                                         J, iters, N_eki, include_theta_sig,
                                         prior_sd_override = Dict(:mu_gbar => sd))
        diag = g_swiid_prior_diagnostics(ca.recs_all, ca.res.t, ca.res.X, window;
                                         rng = Xoshiro(seed + 999))
        println("  事前平均バイアス = ", round(diag.bias, digits = 4),
                "  事前90%被覆 = ", round(diag.coverage, digits = 3),
                " (n=$(diag.n))")
        push!(results, Dict(
            "mu_gbar_sd" => sd,
            "prior_mean_bias" => diag.bias,
            "prior_coverage" => diag.coverage,
            "n" => diag.n,
            "theta_hat" => Dict(string(CAL_PARAMS[k].name) => ca.theta_hat[k]
                                for k in eachindex(CAL_PARAMS)),
            "nu_star" => ca.nu_star,
            "r_hat" => isfinite(ca.r_hat) ? ca.r_hat : nothing))
    end
    elapsed = time() - t_start
    println("  所要時間 $(round(elapsed, digits = 1)) 秒($(length(candidates)) 候補)")

    return Dict(
        "country" => country,
        "t_k" => t_k,
        "window" => collect(window),
        "smoke" => smoke,
        "mu_gbar_anchor" => anchor,
        "candidates" => candidates,
        "results" => results,
        "elapsed_sec" => elapsed,
        "provenance" => Dict(
            "commit" => strip(read(`git -C $(dirname(@__DIR__)) rev-parse HEAD`, String)),
            "seed" => seed,
            "generated_at" => Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"),
            "frozen_decisions" => frozen_decisions_string(),
            "design_decision" => "#0062"))
end

function parse_candidates(spec::AbstractString)
    return [parse(Float64, strip(s)) for s in split(spec, ",") if !isempty(strip(s))]
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    smoke = "--smoke" in ARGS
    candidates = [0.1, 0.2, 0.3]
    for (i, a) in enumerate(ARGS)
        a == "--candidates" && (global candidates = parse_candidates(ARGS[i + 1]))
    end
    out = run_prior_select(; candidates, smoke)
    suffix = smoke ? "_smoke" : ""
    path = joinpath(@__DIR__, "output", "M10_prior_select$(suffix).json")
    mkpath(dirname(path))
    write(path, JSON3.write(out))
    println("保存: $path")
end
