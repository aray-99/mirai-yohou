# M11 真の将来予報モード(Issue #4)
#
# データ末尾(t_k)を起点に 30 年予報アンサンブルを生成し、M12 ダッシュボード
# 用 JSON(experiments/output/M11_forecast_{ISO3}.json)を出力する。
#
# 設計根拠(docs/DECISIONS.md):
#   - #0052: expanding-window walk-forward の較正/同化手続きの土台(run_origin)。
#   - #0053: 予報アンサンブルにおける拡大パラメータ(#0046)の扱い — メンバー別
#     事後値を初期値として持ち込み、予報区間中も同一 rw_sd の RW を継続する。
#   - #0062/#0063: EKI prior の mu_gbar アンカリング(窓内 g_swiid 観測の logit
#     平均)と mu_gbar prior sd = 0.3 の凍結値(M10 既定 0.2 は使わない)。
#   - #0068: 予報窓の内生 Γ ジャンプへの超過確率シンニング p_ex = (同化窓内の
#     強制ジャンプ週数)/(同化窓内のカウントデータ週数)。
#   - #0069: M10 検証ラン合格(walk-forward 回帰基準として不変更のまま踏襲)。
#   - #0071: 将来予報モード専用の p_ex 規則。カウント週ゼロ(0/0)のフォール
#     バックは 1.0 でなく 0.0(JPN の λ_e カスケード・アーティファクト判定)。
#     発動時は provenance と jump_thinning の両方に p_ex_fallback を記録する。
#
# 参考実装(必読): experiments/M11_horizon_smoke.jl(同化→予報の骨格)、
# experiments/M9_walkforward.jl の run_origin(EKI→NegBin 同化→p_ex→予報の
# 完全な手続き)、experiments/M10_walkforward.jl(mu_gbar アンカリング)。
#
# 実行例:
#   julia --project=experiments -t 8 experiments/M11_forecast.jl JPN THA
#   julia --project=experiments -t 8 experiments/M11_forecast.jl THA --horizon 30 --N 100
#   julia --project=experiments -t 8 experiments/M11_forecast.jl JPN --smoke

include(joinpath(@__DIR__, "M10_walkforward.jl"))   # M9/M8_calibrate/M8_hindcast も連鎖 include 済み

using MiraiYohou: sigmoid, deborah_number, intensity

const M11_SEED = 20260711   # #0063 凍結シード(将来予報モードでも固定)

"""
    calendar_year(t) -> Int

保持座標の時刻 `t`(build_observations.jl の t = (year − 1990) + 0.5 の逆変換、
1990 年起点・年央近似)から西暦年へ戻す。
"""
calendar_year(t::Real) = Int(round(1990 + t - 0.5))

"""
    year_idx(y) -> Int

予報開始(y=0)からの経過年 `y`(整数)に対応する予報グリッドのインデックス
(dt=0.01 固定。forecast_ensemble の既定 dt と一致)。y=0 は起点断面
(初期条件 = 同化窓末尾 t_k のアンサンブル)。
"""
year_idx(y::Integer) = 1 + round(Int, y / 0.01)

"1990 年値(なければ最初の観測年値)を基準値として返す(build_observations.jl の _load_series/_baseline を再利用)"
_series_base(country, file) = _baseline(_load_series(_series_path(country, file))...)

"観測レコード列 `recs` 中の変数 `name` の最終観測年(西暦)。観測が無ければ nothing"
function _last_obs_year(recs, name::Symbol)
    ts = [r.t for r in recs if r.spec.name === name]
    isempty(ts) && return nothing
    return calendar_year(maximum(ts))
end

"変数別の自然単位・変換説明(FORECAST_JSON.md §variables と対応。REAL_SERIES の var をキーとする)"
const VAR_META = Dict(
    "P"       => (unit = "persons",
                  transform = "base_P * exp(xi_P)(1990年人口を基準とした対数正規変換の逆変換)"),
    "w"       => (unit = "fraction (0-1)",
                  transform = "sigmoid(xi_w)(logit座標の逆変換)"),
    "y"       => (unit = "USD (2015, per capita)",
                  transform = "base_y * exp(h_logy(xi))(合成観測 y の観測演算子 _h_logy の逆合成、SPEC §9.1)"),
    "g_swiid" => (unit = "fraction (0-1)",
                  transform = "sigmoid(xi_g)(logit座標の逆変換)"),
    "T_proxy" => (unit = "patent applications (resident)",
                  transform = "base_T * exp(xi_T)(1990年値を基準とした対数正規変換の逆変換)"),
    "phi"     => (unit = "fraction (0-1)",
                  transform = "sigmoid(xi_phi)(logit座標の逆変換)"),
    "v"       => (unit = "subscriptions per 100 people",
                  transform = "base_v * exp(xi_v)(1990年値を基準とした対数正規変換の逆変換)"),
    "tau"     => (unit = "fraction (0-1)",
                  transform = "sigmoid(xi_tau)(logit座標の逆変換)"),
    "p"       => (unit = "fraction (0-1)",
                  transform = "sigmoid(xi_p)(logit座標の逆変換)"),
)

"""
    variable_series(country, fe, params, recs_all, s, year_idxs) -> Dict

REAL_SERIES の1エントリ `s` について、予報アンサンブル `fe` の各年次グリッド点
(`year_idxs`)でメンバー横断の自然単位分位点 q05/q25/q50/q75/q95 を計算する
(#0071 実装仕様: 保持座標→自然単位はモデル座標変換の規約 src/coordinates.jl
に従う。log 系列は 1990 年基準値 × exp、logit 系列は sigmoid、y は合成観測
h_logy の指数変換)。
"""
function variable_series(country, fe, params, recs_all, s, year_idxs)
    kind = s.kind
    N = size(fe.X, 3)
    nyears = length(year_idxs)
    h_logy = kind === :log_y ? _h_logy(params) : nothing
    base = kind in (:log_norm, :log_y) ? _series_base(country, s.file) : nothing
    q = Dict(k => Vector{Float64}(undef, nyears) for k in ("q05", "q25", "q50", "q75", "q95"))
    vals = Vector{Float64}(undef, N)
    for (yi, idx) in enumerate(year_idxs)
        for j in 1:N
            if kind === :log_norm
                vals[j] = base * exp(fe.X[s.target, idx, j])
            elseif kind === :log_y
                vals[j] = base * exp(h_logy(view(fe.X, 1:N_STATE, idx, j)))
            else   # :logit_pct または :logit(いずれも状態は logit 座標。sigmoid で 0-1 の自然値)
                vals[j] = sigmoid(fe.X[s.target, idx, j])
            end
        end
        q["q05"][yi] = quantile(vals, 0.05)
        q["q25"][yi] = quantile(vals, 0.25)
        q["q50"][yi] = quantile(vals, 0.5)
        q["q75"][yi] = quantile(vals, 0.75)
        q["q95"][yi] = quantile(vals, 0.95)
    end
    meta = VAR_META[s.var]
    return merge(Dict("unit" => meta.unit, "transform" => meta.transform,
                      "last_observation_year" => _last_obs_year(recs_all, Symbol(s.var))), q)
end

"""
    de_series(fe, params, year_idxs) -> Dict

社会的デボラ数 De(診断量、SPEC §7)の年次分位点。De は状態変数ではないので
`last_observation_year` は null。De の構成パラメータ(eta_g/g_c/eta_p/delta_sig)
は L3 拡大の対象外なので較正後 `params.l2` の共有値でよい。
"""
function de_series(fe, params, year_idxs)
    N = size(fe.X, 3)
    nyears = length(year_idxs)
    q = Dict(k => Vector{Float64}(undef, nyears) for k in ("q05", "q25", "q50", "q75", "q95"))
    vals = Vector{Float64}(undef, N)
    for (yi, idx) in enumerate(year_idxs)
        for j in 1:N
            vals[j] = deborah_number(params.l2, sigmoid(fe.X[IX_G, idx, j]), sigmoid(fe.X[IX_PP, idx, j]))
        end
        q["q05"][yi] = quantile(vals, 0.05)
        q["q25"][yi] = quantile(vals, 0.25)
        q["q50"][yi] = quantile(vals, 0.5)
        q["q75"][yi] = quantile(vals, 0.75)
        q["q95"][yi] = quantile(vals, 0.95)
    end
    return merge(Dict("unit" => "dimensionless",
                      "transform" => "deborah_number(params.l2, sigmoid(xi_g), sigmoid(xi_p))(SPEC §7 診断量)",
                      "last_observation_year" => nothing), q)
end

"""
    count_forecast_series(fe, params, nu_star, r_hat, H) -> Dict

各予報年 y=1..H について、年次期待イベント数 nu_star·Λ_j(intensity の年次
窓積分、stride なし)のメンバー横断分位点を計算する(騒乱リスクのカウント
予報。実現カウントの分位ではなく期待値の分位である点に注意 — NegBin 分散は
`r_hat` 参照)。
"""
function count_forecast_series(fe, params, nu_star, r_hat, H::Int; dt::Float64 = 0.01)
    N = size(fe.X, 3)
    wsteps = round(Int, 1.0 / dt)
    q = Dict(k => Vector{Float64}(undef, H) for k in ("q05", "q25", "q50", "q75", "q95"))
    events = Vector{Float64}(undef, N)
    for y in 1:H
        i0 = (y - 1) * wsteps + 1
        for j in 1:N
            acc = 0.0
            for i in i0:(i0 + wsteps - 1)
                acc += intensity(view(fe.X, 1:N_STATE, i, j), params) * dt
            end
            events[j] = nu_star * acc
        end
        q["q05"][y] = quantile(events, 0.05)
        q["q25"][y] = quantile(events, 0.25)
        q["q50"][y] = quantile(events, 0.5)
        q["q75"][y] = quantile(events, 0.75)
        q["q95"][y] = quantile(events, 0.95)
    end
    return merge(Dict(
        "unit" => "expected political disorder events per year, nu* × ∫intensity",
        "nu_star" => nu_star,
        "r_hat" => isfinite(r_hat) ? r_hat : nothing,
        "note" => "NegBin 分散は r_hat 参照(実現カウントの分位ではなく期待値の分位)"),
        q)
end

"REAL_SERIES 各 file + \"events\" の raw meta.json から fetched_at を集める(#0067/#0068 流儀の来歴収集)。
手動配置系列(tau 等)は fetched_at の代わりに placed_at を持つのでフォールバックで拾う。"
function data_fetched_at_map(country)
    m = Dict{String, Any}()
    files = vcat([(s.var, s.file) for s in REAL_SERIES], [("events", "events")])
    for (var, file) in files
        path = _series_path(country, file) * ".meta.json"
        isfile(path) || continue
        meta = JSON3.read(read(path, String))
        if haskey(meta, :fetched_at)
            m[var] = meta[:fetched_at]
        elseif haskey(meta, :placed_at)
            m[var] = meta[:placed_at]
        end
    end
    return m
end

"""
    run_forecast(country; N, seed, horizon, J, iters, N_eki) -> Dict

将来予報モードの1国分手続き(冒頭コメント参照)。データ末尾を起点に
EKI 最終オリジン再較正 → 同化 → p_ex(#0071)→ 30 年予報アンサンブルを実行し、
M12 ダッシュボード用の JSON(化前 Dict)を返す。
"""
function run_forecast(country::String; N::Int = 100, seed::Integer = M11_SEED,
                      horizon::Float64 = 30.0, J::Int = 24, iters::Int = 4,
                      N_eki::Int = 100)
    ccfg = COUNTRY_CFG[country]

    # 1. 予報起点 t_k = データ末尾
    params0 = build_params(ccfg.regime)
    recs_probe = build_observations(country, params0; t1 = 1e6)
    t_k = maximum(r.t for r in recs_probe)
    window = (ccfg.calib[1], t_k)
    seed_o = seed + Int(floor(t_k))
    println("== $country M11 将来予報: 予報起点 t_k = $t_k (較正窓 $window, 予報 $horizon 年, N=$N, seed=$seed) ==")

    # 2. theta_sig 規則(#0052/#0054)。窓は本ランの同化窓
    cfg0 = AssimConfig(t0 = 0.0, t1 = t_k)
    total_counts = windowed_count_total(build_obs_counts(country, cfg0),
                                        cfg0.t0, cfg0.event_window, window)
    include_theta_sig = total_counts >= THETA_SIG_COUNT_MIN
    println("  theta_sig 規則: ΣN = $total_counts (しきい値 $THETA_SIG_COUNT_MIN) → ",
            include_theta_sig ? "拡大集合に含める" : "除外")

    # 3. EKI 最終オリジン再較正(M10 #0062 規則。mu_gbar prior sd = 0.3 #0063)
    frozen = load_calibrated(country)
    frozen === nothing &&
        error("run_forecast: $country の M8 凍結値がありません(M8_frozen_config.toml)")
    frozen_center = Dict(CAL_PARAMS[k].name => Float64(getproperty(frozen, CAL_PARAMS[k].name))
                         for k in eachindex(CAL_PARAMS))
    anchor = gbar_anchor(country, window)
    prior_center = merge(frozen_center, Dict(:mu_gbar => anchor))
    println("  EKI 再較正: mu_gbar prior 中心(アンカリング) = ", round(anchor, digits = 3),
            "  sd = 0.3(#0063)")
    calib = calibrate(country; J, iters, N = N_eki, seed = seed_o + 101, window,
                      prior_center, prior_sd = 0.5,
                      prior_sd_override = Dict(:mu_gbar => 0.3), save = false,
                      include_theta_sig, validate_prior = true)
    theta_hat, nu_eki = calib.theta_hat, calib.nu_star
    theta_center = Dict(CAL_PARAMS[k].name => theta_hat[k] for k in eachindex(CAL_PARAMS))
    println("  較正値 θ̂ = ", round.(theta_hat, digits = 3), "  ν(EKI) = ", round(nu_eki, digits = 2))

    # 4. params 構築(データ末尾 t_k までの全観測を使用)
    recs_all = build_observations(country, params0; t1 = t_k)
    params = fit_exogenous(build_params(ccfg.regime; theta_center...), recs_all, window)
    recs_calib = recs_all   # t1 = t_k なので全件が較正・同化対象

    # 5. 同化(run_origin (b) と同一、M11_horizon_smoke.jl と同一 AssimConfig)
    cfg = AssimConfig(t0 = 0.0, t1 = t_k, smoother_lag = 5.0,
                      smoother_vars = [IX_G, IX_TAU, IX_SIG, IX_PP],
                      tauA_pseudo_sd_mult = 3.0,
                      analysis_masked_vars = [IX_TAUA],
                      analysis_unmask_names = [:tau, :tauA_pseudo],
                      rtps_alpha = 0.85,
                      obs_spread_floor_frac = 0.5,
                      rejuvenation_a = REJUVENATION_A)   # #0057
    E0_state = initial_ensemble(country, params, recs_all; N, seed = seed_o + 1)
    aug = build_m8_augmented_params(params, country; include_theta_sig)
    E0 = augment_ensemble(E0_state, aug; rng = Xoshiro(seed_o + 6))
    obs_counts = build_obs_counts(country, cfg)
    event_times = filter(t -> t < cfg.t1,
                         build_forced_jumps(country; calib_window = window))
    res = run_assimilation(params, E0, recs_calib, event_times;
                           cfg, seed = seed_o, obs_counts, count_scale = nu_eki,
                           count_temper = 1 / nu_eki, augmented_params = aug)
    ks = count_windows_in(obs_counts, cfg, window)
    if isempty(ks)
        nu_star, r_hat = nu_eki, Inf
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

    # 6. p_ex(#0071 凍結規則 — #0068 と異なるのはフォールバックのみ)
    n_forced_window = count(t -> window[1] <= t < window[2], event_times)
    p_ex_fallback = isempty(ks) ? "no_count_data" : nothing
    p_ex = isempty(ks) ? 0.0 : n_forced_window / length(ks)
    println("  p_ex(Γ シンニング, #0071) = ", round(p_ex, digits = 4),
            "  (強制ジャンプ週 $n_forced_window / カウントデータ週 $(length(ks)))",
            p_ex_fallback === nothing ? "" : "  [フォールバック: $p_ex_fallback]")

    # 7. 30年予報
    t_fore = time()
    fe = forecast_ensemble(params, aug, res; horizon, seed = seed_o + 7,
                           gamma_thinning_p = p_ex)
    elapsed_fore = time() - t_fore
    nonfinite = count(!isfinite, fe.X)
    finite_ok = nonfinite == 0
    println("  予報積分完了: $(round(elapsed_fore, digits = 1)) 秒  NaN/Inf 検査: ",
            finite_ok ? "PASS" : "FAIL($nonfinite 個)")

    # --- 出力 JSON の組み立て ---
    H = round(Int, horizon)
    year_idxs = [year_idx(y) for y in 0:H]
    t_grid = [t_k + y for y in 0:H]
    years_cal = [calendar_year(t) for t in t_grid]

    variables = Dict{String, Any}()
    for s in REAL_SERIES
        variables[s.var] = variable_series(country, fe, params, recs_all, s, year_idxs)
    end
    variables["De"] = de_series(fe, params, year_idxs)

    count_forecast = count_forecast_series(fe, params, nu_star, r_hat, H)

    njumps = sum(length, fe.jumps)
    println("  内生ジャンプ: 合計 $njumps 発(全メンバー・$horizon 年)")

    verified_until = t_k + 6
    out = Dict{String, Any}(
        "country" => country,
        "regime" => string(ccfg.regime),
        "N" => N,
        "horizon_years" => H,
        "forecast_start" => Dict("t" => t_k, "calendar_year" => calendar_year(t_k)),
        "verified_horizon" => Dict(
            "years" => 6,
            "until_t" => verified_until,
            "until_calendar_year" => calendar_year(verified_until),
            "note" => "#0052/#0069 の検証範囲(1年先×全オリジン + M10 合格)。これ以遠は外挿領域"),
        "finite_check" => Dict("pass" => finite_ok, "nonfinite_values" => nonfinite),
        "assimilation" => Dict(
            "theta_hat" => Dict(string(k) => v for (k, v) in theta_center),
            "nu_star" => nu_star,
            "r_hat" => isfinite(r_hat) ? r_hat : nothing,
            "include_theta_sig" => include_theta_sig,
            "sigma_n_total_counts" => total_counts,
            "mu_gbar_prior" => Dict("center" => anchor, "sd" => 0.3, "rule" => "#0062"),
            "nresample" => res.nresample,
            "ess_range" => collect(extrema(res.ess)),
            "n_obs" => length(recs_calib),
            "window" => collect(window)),
        "jump_thinning" => Dict(
            "p_ex" => p_ex,
            "n_forced_weeks" => n_forced_window,
            "n_count_weeks" => length(ks),
            "p_ex_fallback" => p_ex_fallback,
            "design_decision" => "#0068 式 / #0071 将来予報フォールバック規則"),
        "variables" => variables,
        "count_forecast" => count_forecast,
        "years" => years_cal,
        "t" => t_grid,
        "endogenous_jumps_total" => njumps,
    )
    return out
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    horizon = 30.0
    horizon_explicit = false
    N = 100
    smoke = "--smoke" in ARGS
    for (i, a) in enumerate(ARGS)
        if a == "--horizon"
            global horizon = parse(Float64, ARGS[i + 1])
            global horizon_explicit = true
        end
        a == "--N" && (global N = parse(Int, ARGS[i + 1]))
    end
    J, iters, N_eki = 24, 4, 100
    if smoke
        N = 40; J = 12; iters = 2; N_eki = 40
        horizon_explicit || (global horizon = 5.0)
    end
    seed = M11_SEED   # #0063 凍結シード
    countries = [c for c in ARGS if haskey(COUNTRY_CFG, c)]
    isempty(countries) && (countries = ["JPN", "THA"])
    for c in countries
        t0 = time()
        out = run_forecast(c; N, seed, horizon, J, iters, N_eki)
        elapsed = time() - t0
        provenance = Dict(
            "commit" => strip(read(`git -C $(dirname(@__DIR__)) rev-parse HEAD`, String)),
            "seed" => seed,
            "seed_offsets" => "seed_o = seed + floor(t_k)(run_origin 規約)。EKI: seed_o+101、初期アンサンブル: seed_o+1、拡大アンサンブル rng: seed_o+6、同化: seed_o、予報: seed_o+7",
            "generated_at" => Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"),
            "frozen_decisions" => frozen_decisions_string(),
            "frozen_config" => "M8_frozen_config.toml",
            "design_decision" => "Issue #4; p_ex rule #0071; forecast rules #0053/#0068; anchoring #0062/#0063; seed #0063",
            "script" => "experiments/M11_forecast.jl",
            "data_fetched_at" => data_fetched_at_map(c),
            "p_ex_fallback" => out["jump_thinning"]["p_ex_fallback"],
        )
        out["provenance"] = provenance
        out["elapsed_sec"] = elapsed
        suffix = smoke ? "_smoke" : ""
        path = joinpath(@__DIR__, "output", "M11_forecast_$(c)$(suffix).json")
        mkpath(dirname(path))
        write(path, JSON3.write(out))
        write(path * ".meta.json", JSON3.write(provenance))   # 来歴サイドカー(M11_horizon_smoke.jl と同じ流儀)
        println("保存: $path(+ .meta.json)  所要時間 $(round(elapsed, digits = 1)) 秒")
    end
end
