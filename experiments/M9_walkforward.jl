# M9 expanding-window walk-forward ドライバ(DECISIONS #0052)
#
# JPN/THA について年次オリジン t_k ごとに:
#   (a) 窓開始(JPN t=5 / THA t=20)〜t_k のデータのみで EKI 再較正
#       (warm-start: 初回は M8 凍結値、以後は前オリジンの較正値。prior sd
#        は毎回既定値 0.5 に復元)
#   (b) t_k まで同化(M8 と同一の AssimConfig・拡大集合)
#   (c) 1年先予報アンサンブルを発行し、[t_k, t_k+1) の観測で評価
#       (事前90%被覆・変数別 RMSE・THA 週次カウント予測尤度)。
#       自由ラン基準(同一較正値・同化なし・同区間)も同時に計算する。
# 検証区間の観測は較正・同化とも未参照(疑似リアルタイム、#0052)。
# 全オリジンをプールして国別の基準1〜3判定値を集計する。
#
# 実行: julia --project=experiments -t 8 experiments/M9_walkforward.jl JPN THA
#         [--smoke] [--origins a,b]
#   --smoke: 各国先頭2オリジン・N=40・EKI J=12/iters=2 の動作確認モード。
#   --origins a,b: 既定のオリジン年列を上書き(国ごとに有効な年のみ採用)。

using Dates
using Random
using Statistics
using JSON3
using MiraiYohou: N_STATE, member_seed, run_assimilation, free_ensemble,
                  AssimConfig, augment_ensemble, simulate_sde, intensity

include(joinpath(@__DIR__, "M8_calibrate.jl"))   # M8_hindcast.jl も連鎖 include 済み

"""
国別オリジン列(年次、#0052): JPN t=26..33(8本)、THA t=28..33(6本)。
窓開始(expanding の左端)は `COUNTRY_CFG[country].calib[1]` を流用する
(JPN t=5 / THA t=20、M8 較正窓の開始と同じ)。
"""
const M9_ORIGINS = Dict(
    "JPN" => collect(26:33),
    "THA" => collect(28:33),
)

const M9_HORIZON = 1.0    # 1年先予報(#0052)

"""
    forecast_ensemble(params, res; horizon, seed, dt) -> (; t, X)

同化結果 `res` の最終時刻(オリジン t_k)のアンサンブル状態から、内生 Hawkes
のみ(強制ジャンプなし・同化なし)で `horizon` 年だけ純粋に前進させる
1年先予報アンサンブル(#0052)。`forecast_rmse`(M8_hindcast.jl)と同じ
`simulate_sde` の既定モード(EndogenousHawkes)を使う。
"""
function forecast_ensemble(params, res; horizon::Float64 = M9_HORIZON,
                           seed::Integer, dt::Float64 = 0.01)
    N = size(res.X, 3)
    t0 = res.t[end]
    t1 = t0 + horizon
    nsteps = round(Int, (t1 - t0) / dt)
    X = Array{Float64,3}(undef, N_STATE, nsteps + 1, N)
    xi0_all = view(res.X, 1:N_STATE, size(res.X, 2), :)
    Threads.@threads for j in 1:N
        sim = simulate_sde(params; seed = member_seed(seed, j),
                           t0 = t0, t1 = t1, dt = dt,
                           xi0 = collect(view(xi0_all, :, j)))
        X[:, :, j] = sim.traj.X
    end
    ts = collect(range(t0; step = dt, length = nsteps + 1))
    return (; t = ts, X)
end

"""
    run_origin(country, t_k, prior_center; N, seed, J, iters, N_eki) -> NamedTuple

1オリジン分の (較正 → 同化 → 予報 → 自由ラン → 評価)。戻り値に次オリジンへ
の warm-start 用 `theta_center`(Dict)と評価指標一式を含む。
"""
function run_origin(country::String, t_k::Real,
                    prior_center::Union{Nothing,Dict};
                    N::Int, seed::Integer, J::Int, iters::Int, N_eki::Int)
    ccfg = COUNTRY_CFG[country]
    win_start = ccfg.calib[1]
    window = (win_start, Float64(t_k))
    t_fore_end = min(Float64(t_k) + M9_HORIZON, T1)
    eval_win = (Float64(t_k), t_fore_end)

    println("-- $country origin t=$(t_k) (calib window $window, eval $eval_win) --")

    # (a) EKI 再較正(expanding window, warm-start, prior sd 復元)
    calib = calibrate(country; J, iters, N = N_eki, seed = seed + 101,
                      window, prior_center, prior_sd = 0.5, save = false)
    theta_hat, nu_star = calib.theta_hat, calib.nu_star
    theta_center = Dict(CAL_PARAMS[k].name => theta_hat[k] for k in eachindex(CAL_PARAMS))
    kw = theta_center
    println("  較正値 θ̂ = ", round.(theta_hat, digits = 3), "  ν* = ", round(nu_star, digits = 2))

    # 評価区間終端までの観測(較正・同化には window[2]=t_k 以前のみ使用)
    params0 = build_params(ccfg.regime)
    recs_all = build_observations(country, params0; t1 = t_fore_end)
    params = fit_exogenous(build_params(ccfg.regime; kw...), recs_all, window)
    recs_calib = [r for r in recs_all if r.t <= t_k]

    # (b) t_k までの同化(M8 検証ランと同一の AssimConfig)
    cfg = AssimConfig(t0 = 0.0, t1 = Float64(t_k), smoother_lag = 5.0,
                      smoother_vars = [IX_G, IX_TAU, IX_SIG, IX_PP],
                      tauA_pseudo_sd_mult = 3.0,
                      analysis_masked_vars = [IX_TAUA],
                      analysis_unmask_names = [:tau, :tauA_pseudo],
                      rtps_alpha = 0.85,
                      obs_spread_floor_frac = 0.5)
    E0_state = initial_ensemble(country, params, recs_all; N, seed = seed + 1)
    aug = build_m8_augmented_params(params, country)
    E0 = augment_ensemble(E0_state, aug; rng = Xoshiro(seed + 6))
    obs_counts = build_obs_counts(country, cfg)
    event_times = filter(t -> t < cfg.t1,
                         build_forced_jumps(country; calib_window = window))
    res = run_assimilation(params, E0, recs_calib, event_times;
                           cfg, seed, obs_counts, count_scale = nu_star,
                           count_temper = 1 / nu_star, augmented_params = aug)

    # (c) 1年先予報アンサンブル
    fe = forecast_ensemble(params, res; horizon = t_fore_end - t_k, seed = seed + 7)

    # 自由ラン基準(同一較正値・同化なし・同区間、#0052)
    cfg_free = AssimConfig(t0 = 0.0, t1 = t_fore_end)
    Xf = free_ensemble(params, E0; cfg = cfg_free, seed = seed + 8)
    ts_free = collect(range(cfg_free.t0; step = cfg_free.dt, length = size(Xf, 2)))

    # 評価(forecast は事後の再解析を含まないため lag=0。M8 の lag=1(事前近似)
    # とは異なり、予報・自由ランのいずれも解析更新を経ないので直近グリッド点
    # をそのまま使ってよい)
    rng_cov = Xoshiro(seed + 9)
    cov_f, n_f, hit_f = coverage(recs_all, fe.t, fe.X, eval_win; rng = rng_cov, lag = 0)
    rng_cov2 = Xoshiro(seed + 10)
    cov_r, n_r, hit_r = coverage(recs_all, ts_free, Xf, eval_win; rng = rng_cov2, lag = 0)

    err_fore = obs_errors(recs_all, fe.t, fe.X, eval_win; lag = 0)
    err_free = obs_errors(recs_all, ts_free, Xf, eval_win; lag = 0)

    # count_loglik の窓インデックスは各系列自身の t0(fe: t_k、free: 0)基準
    # なので、obs_counts も系列ごとに別グリッドで構築する(同一カレンダー
    # 窓だが index 対応が異なる — 取り違えると窓ずれで logL が壊れる)。
    cfg_fore_eval = AssimConfig(t0 = t_k, t1 = t_fore_end)
    obs_counts_fore = build_obs_counts(country, cfg_fore_eval)
    ll_fore, nwin_fore = count_loglik(fe.t, fe.X, obs_counts_fore, params, nu_star,
                                      cfg_fore_eval, eval_win)
    obs_counts_free = build_obs_counts(country, cfg_free)
    ll_free, _ = count_loglik(ts_free, Xf, obs_counts_free, params, nu_star,
                              cfg_free, eval_win)

    println("  被覆(予報) = ", round(cov_f, digits = 3), " (n=$n_f)  ",
            "被覆(自由) = ", round(cov_r, digits = 3), " (n=$n_r)")
    for k in sort(collect(keys(err_fore)); by = string)
        haskey(err_free, k) || continue
        rmse_f = sqrt(mean(abs2, err_fore[k]))
        rmse_r = sqrt(mean(abs2, err_free[k]))
        println("  RMSE $(rpad(k, 9)) 予報 ", round(rmse_f, digits = 4),
                "  自由 ", round(rmse_r, digits = 4))
    end
    if nwin_fore > 0
        println("  カウント予測 logL 予報 ", round(ll_fore, digits = 2),
                "  自由 ", round(ll_free, digits = 2))
    end

    return (; t_k = Float64(t_k), theta_center, nu_star,
            cov_fore = (hit = hit_f, n = n_f), cov_free = (hit = hit_r, n = n_r),
            err_fore, err_free, ll_fore, ll_free, nwin = nwin_fore,
            resample = res.nresample, ess = extrema(res.ess))
end

"""
    run_walkforward(country; N, seed, J, iters, N_eki, smoke, origins) -> Dict

国1つ分の walk-forward 全オリジンを実行し、プール判定値を含む診断 Dict を
返す(JSON 化前の中間形)。
"""
function run_walkforward(country::String; N::Int = 100, seed::Integer = 20260711,
                         J::Int = 24, iters::Int = 4, N_eki::Int = 100,
                         smoke::Bool = false, origins::Union{Nothing,Vector{Int}} = nothing)
    orig_list = origins !== nothing ? [t for t in origins if t in M9_ORIGINS[country]] :
                copy(M9_ORIGINS[country])
    if smoke
        orig_list = first(orig_list, 2)
        N = 40; J = 12; iters = 2; N_eki = 40
    end
    println("== $country M9 walk-forward: origins = $orig_list (N=$N, J=$J, iters=$iters) ==")

    # 初回オリジンの prior 中心 = M8 凍結値(#0050)で warm-start(#0052)。
    # 凍結値が無い場合(未凍結国)は nothing(EKI 既定のモーメント初期化)。
    frozen = load_calibrated(country)
    prior_center = frozen === nothing ? nothing :
        Dict(CAL_PARAMS[k].name => Float64(getproperty(frozen, CAL_PARAMS[k].name))
             for k in eachindex(CAL_PARAMS))

    origin_results = []
    t_start = time()
    for t_k in orig_list
        r = run_origin(country, t_k, prior_center; N, seed = seed + t_k, J, iters, N_eki)
        push!(origin_results, r)
        prior_center = r.theta_center    # 次オリジンへ warm-start(#0052)
    end
    elapsed = time() - t_start

    # プール集計(#0052: 全オリジンをプールした国別判定値)
    hit_f = sum(r.cov_fore.hit for r in origin_results)
    n_f = sum(r.cov_fore.n for r in origin_results)
    hit_r = sum(r.cov_free.hit for r in origin_results)
    n_r = sum(r.cov_free.n for r in origin_results)
    cov_pooled = hit_f / max(n_f, 1)

    err_fore_pooled = Dict{Symbol, Vector{Float64}}()
    err_free_pooled = Dict{Symbol, Vector{Float64}}()
    for r in origin_results
        for (k, v) in r.err_fore; append!(get!(err_fore_pooled, k, Float64[]), v); end
        for (k, v) in r.err_free; append!(get!(err_free_pooled, k, Float64[]), v); end
    end
    rmse_fore = Dict(k => sqrt(mean(abs2, v)) for (k, v) in err_fore_pooled)
    rmse_free = Dict(k => sqrt(mean(abs2, v)) for (k, v) in err_free_pooled)
    nbetter = 0; nvars = 0
    var_detail = Dict{String, Any}()
    for k in sort(collect(keys(rmse_fore)); by = string)
        haskey(rmse_free, k) || continue
        nvars += 1
        better = rmse_fore[k] < rmse_free[k]
        nbetter += better
        var_detail[string(k)] = Dict("rmse_forecast" => rmse_fore[k],
                                     "rmse_free" => rmse_free[k], "better" => better,
                                     "n" => length(err_fore_pooled[k]))
    end

    ll_fore_total = sum(r.ll_fore for r in origin_results)
    ll_free_total = sum(r.ll_free for r in origin_results)
    nwin_total = sum(r.nwin for r in origin_results)

    crit1 = 0.80 <= cov_pooled <= 0.98
    crit2 = nvars > 0 && nbetter > nvars ÷ 2
    crit3 = nwin_total > 0 ? ll_fore_total > ll_free_total : nothing

    println("== $country プール判定 ==")
    println("  [基準1] 被覆率 = ", round(cov_pooled, digits = 3), " (n=$n_f, hit=$hit_f)",
            crit1 ? "  PASS" : "  FAIL")
    println("  [基準2] RMSE改善 $nbetter/$nvars", crit2 ? "  PASS" : "  FAIL")
    if crit3 !== nothing
        println("  [基準3] カウント予測logL 予報 ", round(ll_fore_total, digits = 1),
                "  自由 ", round(ll_free_total, digits = 1), crit3 ? "  PASS" : "  FAIL")
    else
        println("  [基準3] N/A(カウント窓なし)")
    end
    println("  所要時間 $(round(elapsed, digits = 1)) 秒($(length(orig_list)) オリジン)")

    return Dict(
        "country" => country,
        "smoke" => smoke,
        "origins" => orig_list,
        "window_start" => COUNTRY_CFG[country].calib[1],
        "horizon" => M9_HORIZON,
        "criteria" => Dict(
            "coverage" => Dict("value" => cov_pooled, "n" => n_f, "hit" => hit_f,
                               "pass" => crit1, "band" => [0.80, 0.98]),
            "rmse_majority_improve" => Dict("nbetter" => nbetter, "nvars" => nvars,
                                            "pass" => crit2),
            "count_loglik" => Dict("forecast" => ll_fore_total, "free" => ll_free_total,
                                   "nwin" => nwin_total,
                                   "pass" => crit3 === nothing ? nothing : crit3)),
        "variables" => var_detail,
        "free_run_coverage" => Dict("value" => hit_r / max(n_r, 1), "n" => n_r, "hit" => hit_r),
        "per_origin" => [Dict(
            "t_k" => r.t_k,
            "theta" => Dict(string(k) => v for (k, v) in r.theta_center),
            "nu_star" => r.nu_star,
            "coverage_forecast" => Dict("hit" => r.cov_fore.hit, "n" => r.cov_fore.n),
            "coverage_free" => Dict("hit" => r.cov_free.hit, "n" => r.cov_free.n),
            "ll_forecast" => r.ll_fore, "ll_free" => r.ll_free, "nwin" => r.nwin,
            "resample" => r.resample, "ess_range" => collect(r.ess))
            for r in origin_results],
        "elapsed_sec" => elapsed,
        "provenance" => Dict(
            "commit" => strip(read(`git -C $(dirname(@__DIR__)) rev-parse HEAD`, String)),
            "seed" => seed,
            "generated_at" => Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"),
            "frozen_decisions" => frozen_decisions_string(),
            "design_decision" => "#0052"))
end

function parse_origins(spec::AbstractString)
    return [parse(Int, strip(s)) for s in split(spec, ",") if !isempty(strip(s))]
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    smoke = "--smoke" in ARGS
    origins = nothing
    for (i, a) in enumerate(ARGS)
        a == "--origins" && (global origins = parse_origins(ARGS[i + 1]))
    end
    countries = [c for c in ARGS if haskey(COUNTRY_CFG, c)]
    isempty(countries) && (countries = ["JPN", "THA"])
    for c in countries
        out = run_walkforward(c; smoke, origins)
        suffix = smoke ? "_smoke" : ""
        path = joinpath(@__DIR__, "output", "M9_walkforward_$(c)$(suffix).json")
        mkpath(dirname(path))
        write(path, JSON3.write(out))
        println("保存: $path")
    end
end
