# カウント尤度の NegBin 評価(DECISIONS #0054、較正窓のみ)
#
# THA 較正窓(t=20〜28)で、M8 凍結設定(M8_frozen_config.toml、N=100)の
# 較正窓同化ランを2条件で実行して比較する:
#   (a) baseline: 現行 = 複合ポアソン近似 + 1/ν テンパリング(#0033)
#   (b) negbin  : NegBin 尤度(平均 ν·Λ、サイズ r はプロファイル最尤 #0054)
#
# 指標(#0054 で事前凍結):
#   (i)  カウント観測の逐次1ステップ先予測 log スコア
#        (更新前アンサンブルの重み付き混合予測分布。AssimResult.count_logscore)
#   (ii) カウント更新時 ESS の最小・中央値
# 採否規則: NegBin 採用は (i) 改善 かつ (ii) ESS 中央値がベースラインの
# 0.8 倍以上の場合のみ。
#
# r のプロファイル: baseline ランの週次 Λ̄_w(アンサンブル平均強度の窓積分)
# と観測 N_w から、平均 ν·Λ̄_w を固定して r を1次元最尤(negbin_profile_r)。
# ν は M8 凍結値(較正時のプロファイル値)を両条件共通で使う。
# 検証窓 [28, 35] は一切参照しない。シードは両条件共通。
#
# 実行: julia --project=experiments -t 8 experiments/M9_negbin_eval.jl

using Dates
using Random
using Statistics
using JSON3
using MiraiYohou: negbin_profile_r, intensity, N_STATE

include(joinpath(@__DIR__, "M8_hindcast.jl"))

"THA 較正窓の同化ラン入力(M8 検証ランと同一機構・凍結値)"
function negbin_eval_inputs(; N::Int = 100, seed::Integer = 20260708)
    country = "THA"
    ccfg = COUNTRY_CFG[country]
    calib = ccfg.calib                              # (20.0, 28.0)
    cal = load_calibrated(country)
    cal !== nothing || error("M8_frozen_config.toml に THA の凍結値がありません")
    kw = Dict(Symbol(k) => Float64(v) for (k, v) in pairs(cal) if k !== :nu)
    recs = build_observations(country, build_params(ccfg.regime); t1 = calib[2])
    params = fit_exogenous(build_params(ccfg.regime; kw...), recs, calib)
    nu = max(Float64(cal.nu), 1.0)
    cfg = AssimConfig(t0 = 0.0, t1 = calib[2], smoother_lag = 5.0,
                      smoother_vars = [IX_G, IX_TAU, IX_SIG, IX_PP],  # #0036
                      tauA_pseudo_sd_mult = 3.0,                     # #0036
                      analysis_masked_vars = [IX_TAUA],               # #0040-(α)
                      analysis_unmask_names = [:tau, :tauA_pseudo],   # #0040-(α)
                      rtps_alpha = 0.85,                              # #0040-(β)
                      obs_spread_floor_frac = 0.5)                    # #0043
    E0_state = initial_ensemble(country, params, recs; N, seed = seed + 1)
    aug = build_m8_augmented_params(params, country)
    E0 = augment_ensemble(E0_state, aug; rng = Xoshiro(seed + 6))
    obs_counts = build_obs_counts(country, cfg)
    event_times = filter(t -> t < cfg.t1, build_forced_jumps(country))
    recs_run = [r for r in recs if r.t <= cfg.t1]
    return (; country, calib, params, nu, cfg, E0, aug, obs_counts, event_times,
            recs_run, seed)
end

"較正窓内でデータのある週次窓のインデックス集合"
function count_window_indices(inp)
    cfg, calib = inp.cfg, inp.calib
    w = cfg.event_window
    return [k for (k, c) in enumerate(inp.obs_counts)
            if c >= 0 && calib[1] <= cfg.t0 + (k - 1) * w < calib[2]]
end

"1条件の同化ラン + 指標抽出(#0054: log スコア・カウント更新時 ESS)"
function run_condition(inp; count_model::Symbol, count_dispersion::Float64 = Inf,
                       label::String)
    println("== 条件: $label ==")
    t = @elapsed res = run_assimilation(inp.params, copy(inp.E0), inp.recs_run,
                                        inp.event_times;
                                        cfg = inp.cfg, seed = inp.seed,
                                        obs_counts = inp.obs_counts,
                                        count_scale = inp.nu,
                                        count_temper = 1 / inp.nu,   # :negbin では無視
                                        augmented_params = inp.aug,
                                        count_model, count_dispersion)
    ks = count_window_indices(inp)
    scores = [res.count_logscore[k] for k in ks]
    esss = [res.ess[k] for k in ks]
    m = Dict(
        "label" => label,
        "count_model" => string(count_model),
        # Inf(= :poisson の未使用値)は JSON 非対応なので nothing にする
        "count_dispersion" => isfinite(count_dispersion) ? count_dispersion : nothing,
        "n_count_windows" => length(ks),
        "logscore_total" => sum(scores),
        "logscore_mean" => mean(scores),
        "ess_min" => minimum(esss),
        "ess_median" => median(esss),
        "nresample" => res.nresample,
        "elapsed_sec" => t)
    println("  カウント窓 $(length(ks))  logスコア合計 ", round(sum(scores), digits = 2),
            "  ESS min/med ", round(minimum(esss), digits = 1), "/",
            round(median(esss), digits = 1),
            "  再抽選 ", res.nresample, " 回  (", round(t, digits = 1), "s)")
    return res, m, ks
end

"baseline ランの週次 Λ̄_w(アンサンブル平均強度の窓積分。count_loglik と同一)"
function window_lambdas(res, inp, ks)
    cfg = inp.cfg
    wsteps = max(1, round(Int, cfg.event_window / cfg.dt))
    N = size(res.X, 3)
    lams = Float64[]
    for k in ks
        acc = 0.0
        i0 = (k - 1) * wsteps + 1
        for i in i0:min(i0 + wsteps - 1, length(res.t))
            acc += mean(intensity(view(res.X, 1:N_STATE, i, j), inp.params)
                        for j in 1:N) * cfg.dt
        end
        push!(lams, acc)
    end
    return lams
end

function main(; N::Int = 100, seed::Integer = 20260708)
    inp = negbin_eval_inputs(; N, seed)
    println("THA 較正窓 $(inp.calib)  N=$N  ν(凍結) = ", round(inp.nu, digits = 3))

    # (a) baseline: 現行ポアソン + 1/ν テンパリング
    res_base, m_base, ks = run_condition(inp; count_model = :poisson,
                                         label = "baseline (Poisson + 1/nu temper)")

    # r のプロファイル最尤(baseline ランの Λ̄_w、平均 ν·Λ̄_w 固定、#0054)
    counts = [inp.obs_counts[k] for k in ks]
    mus = inp.nu .* window_lambdas(res_base, inp, ks)
    r_hat = negbin_profile_r(counts, mus)
    println("プロファイル r̂ = ", round(r_hat, digits = 4),
            "(ΣN = ", sum(counts), ", Σν·Λ̄ = ", round(sum(mus), digits = 1), ")")

    # (b) NegBin(同一シード)
    _, m_nb, _ = run_condition(inp; count_model = :negbin,
                               count_dispersion = r_hat, label = "negbin")

    # 機械的判定(#0054 の採否規則。最終判定はエントリ側)
    logscore_improved = m_nb["logscore_total"] > m_base["logscore_total"]
    ess_ok = m_nb["ess_median"] >= 0.8 * m_base["ess_median"]
    adopt = logscore_improved && ess_ok
    println("== #0054 採否規則の機械的判定 ==")
    println("  (i) logスコア改善: ", logscore_improved,
            " (negbin ", round(m_nb["logscore_total"], digits = 2),
            " vs baseline ", round(m_base["logscore_total"], digits = 2), ")")
    println("  (ii) ESS中央値 ≥ 0.8×baseline: ", ess_ok,
            " (", round(m_nb["ess_median"], digits = 1), " vs 0.8×",
            round(m_base["ess_median"], digits = 1), " = ",
            round(0.8 * m_base["ess_median"], digits = 1), ")")
    println("  → NegBin 採用条件: ", adopt ? "満たす" : "満たさない")

    out = Dict(
        "country" => inp.country,
        "calib_window" => collect(inp.calib),
        "N" => N,
        "nu_frozen" => inp.nu,
        "r_profiled" => r_hat,
        "conditions" => [m_base, m_nb],
        "verdict" => Dict(
            "logscore_improved" => logscore_improved,
            "ess_median_ratio" => m_nb["ess_median"] / m_base["ess_median"],
            "adopt_negbin" => adopt,
            "rule" => "#0054: adopt iff logscore improves AND ess_median >= 0.8x baseline"),
        "provenance" => Dict(
            "commit" => strip(read(`git -C $(dirname(@__DIR__)) rev-parse HEAD`, String)),
            "seed" => seed,
            "generated_at" => Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"),
            "frozen_decisions" => frozen_decisions_string(),
            "design_decision" => "#0054"))
    path = joinpath(@__DIR__, "output", "M9_negbin_eval_THA.json")
    mkpath(dirname(path))
    write(path, JSON3.write(out))
    println("保存: $path")
    return out
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    main()
end
