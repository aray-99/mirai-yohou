# M10 walk-forward ドライバ(DECISIONS #0062): M9(#0052)の expanding-window
# walk-forward ハーネスを土台に、EKI prior の扱いを2点変更する。
#
#   1. warm-start 連鎖の廃止: M9 は初回オリジンのみ M8 凍結値、以後は前
#      オリジンの較正値を prior 中心にしていた(lam0 のオリジン間ドリフト
#      蓄積チャネル、#0056 診断2)。M10 は **全オリジンで prior 中心 = M8
#      凍結値**(experiments/M8_frozen_config.toml)に固定する。
#   2. mu_gbar のアンカリング prior: mu_gbar だけは prior 中心を「当該
#      オリジン較正窓内(窓開始〜t_k)の g_swiid 観測の logit 平均」に係留
#      する(build_observations.jl が g_swiid を既に logit 座標で保持する
#      ので追加変換は不要)。mu_gbar の prior sd も他パラメータ(0.5 維持)
#      と独立に指定できる(--mu-gbar-sd、既定 0.2。確定値は事前選定
#      M10_prior_select.jl の結果を見て設定凍結エントリで決める)。
#
# mu_gbar は EKI 較正の変換空間で `:id`(恒等)なので、l2.mu_gbar 自体が
# 既に logit 座標の平衡値(#0031、M8_calibrate.jl の CAL_PARAMS)。したがって
# アンカリング中心は「窓内 g_swiid observation の logit 値の単純平均」を
# そのまま prior_center[:mu_gbar] に渡せばよい(空間変換不要)。
#
# 検証形式・判定基準・オリジン列・NegBin・theta_sig 規則・若返り a=0.95・
# 予報生成則は #0062 により M9(#0052/#0055/#0057)から不変で引き継ぐ。
#
# 実行: julia --project=experiments -t 8 experiments/M10_walkforward.jl JPN THA
#         [--smoke] [--origins a,b] [--mu-gbar-sd 0.2]
#   --smoke: 各国先頭2オリジン・N=40・EKI J=12/iters=2 の動作確認モード。
#   --mu-gbar-sd: mu_gbar アンカリング prior の sd(変換空間、恒等なので
#     logit 座標の sd と同義)。既定 0.2(仮。#0062 の事前選定プロトコル
#     M10_prior_select.jl の結果を見て設定凍結で確定)。

include(joinpath(@__DIR__, "M9_walkforward.jl"))   # M8_calibrate.jl 等も連鎖 include 済み

"""
    gbar_anchor(country, window) -> Float64

`window = (窓開始, t_k)` 内の g_swiid 観測値(logit 座標、build_observations.jl
の REAL_SERIES 定義により既に logit 変換済み)の単純平均を返す(#0062 規則1)。
`recs_calib`(run_origin 内)と同じ包含規約(`r.t <= t_k`)に合わせ
`window[1] <= r.t <= window[2]` で絞る。観測が1点もない場合はエラー
(#0062 は両国とも g_swiid が較正窓内に存在する前提。JPN/THA いずれも
較正窓開始が SWIID カバレッジ開始以後なので通常発生しない)。
"""
function gbar_anchor(country::String, window)
    ccfg = COUNTRY_CFG[country]
    params0 = build_params(ccfg.regime)
    recs = build_observations(country, params0; t1 = window[2])
    vals = [r.value for r in recs
            if r.spec.name === :g_swiid && window[1] <= r.t <= window[2]]
    isempty(vals) &&
        error("gbar_anchor: $country window=$window に g_swiid 観測がありません")
    return mean(vals)
end

"""
    run_walkforward_m10(country; N, seed, J, iters, N_eki, smoke, origins,
                        mu_gbar_sd) -> Dict

M9 の `run_walkforward` と同じ流れ(オリジンごとに較正 → 同化 → 予報 →
自由ラン → 評価、プール集計)だが、prior 規則を #0062 のとおりに変更する:
warm-start を行わず、毎オリジン `prior_center` を M8 凍結値から再構築し、
`mu_gbar` のみ `gbar_anchor` で上書き・sd を `mu_gbar_sd` で独立指定する。
出力 JSON にはオリジンごとのアンカリング中心・sd を追加で記録する。
"""
function run_walkforward_m10(country::String; N::Int = 100, seed::Integer = 20260711,
                             J::Int = 24, iters::Int = 4, N_eki::Int = 100,
                             smoke::Bool = false,
                             origins::Union{Nothing,Vector{Int}} = nothing,
                             mu_gbar_sd::Float64 = 0.2)
    orig_list = origins !== nothing ? [t for t in origins if t in M9_ORIGINS[country]] :
                copy(M9_ORIGINS[country])
    if smoke
        orig_list = first(orig_list, 2)
        N = 40; J = 12; iters = 2; N_eki = 40
    end
    println("== $country M10 walk-forward: origins = $orig_list (N=$N, J=$J, iters=$iters, mu_gbar_sd=$mu_gbar_sd) ==")

    # theta_sig 拡大の適用規則(#0052/#0054、M9 と同一。M10 で変更しない)
    win_start = COUNTRY_CFG[country].calib[1]
    t_k1 = Float64(first(orig_list))
    cfg0 = AssimConfig(t0 = 0.0, t1 = t_k1)
    total_counts = windowed_count_total(build_obs_counts(country, cfg0),
                                        cfg0.t0, cfg0.event_window, (win_start, t_k1))
    include_theta_sig = total_counts >= THETA_SIG_COUNT_MIN
    println("  theta_sig 規則: 較正窓 [$win_start, $t_k1) のフィルタ後 ΣN = ",
            total_counts, " (しきい値 $THETA_SIG_COUNT_MIN) → ",
            include_theta_sig ? "拡大集合に含める" : "除外")

    # M8 凍結値(#0050)。M10 は全オリジンでこれを prior 中心とし、warm-start
    # 連鎖を行わない(#0062 規則2)。未凍結国は非対応(M10 は THA/JPN 前提)。
    frozen = load_calibrated(country)
    frozen === nothing &&
        error("run_walkforward_m10: $country の M8 凍結値がありません(M8_frozen_config.toml)")
    frozen_center = Dict(CAL_PARAMS[k].name => Float64(getproperty(frozen, CAL_PARAMS[k].name))
                         for k in eachindex(CAL_PARAMS))

    prior_sd_override = Dict(:mu_gbar => mu_gbar_sd)   # #0062 規則3(他は 0.5 のまま)

    origin_results = []
    anchors = Dict{Float64,Float64}()
    t_start = time()
    for t_k in orig_list
        window = (win_start, Float64(t_k))
        anchor = gbar_anchor(country, window)
        anchors[Float64(t_k)] = anchor
        prior_center = merge(frozen_center, Dict(:mu_gbar => anchor))
        println("  オリジン t=$t_k: mu_gbar prior 中心(アンカリング) = ",
                round(anchor, digits = 3), "  sd = $mu_gbar_sd")
        r = run_origin(country, t_k, prior_center; N, seed = seed + t_k, J, iters,
                       N_eki, include_theta_sig, prior_sd_override,
                       validate_prior = true)   # #0064: prior 妥当性リサンプリング
        push!(origin_results, r)
        # 注意: 次オリジンへの warm-start は行わない(#0062 規則2)。
        # prior_center は次のループ反復で frozen_center から再構築される。
    end
    elapsed = time() - t_start

    out = build_walkforward_output(country, orig_list, origin_results, smoke, seed,
                                   elapsed; design_decision = "#0062",
                                   extra_provenance = Dict(
                                       "mu_gbar_prior_sd" => mu_gbar_sd,
                                       "mu_gbar_prior_rule" =>
                                           "logit mean of g_swiid observations in [window_start, t_k] (#0062 rule 1)",
                                       "prior_center_rule" =>
                                           "M8 frozen config every origin, no warm-start chaining (#0062 rule 2)"))
    for entry in out["per_origin"]
        entry["mu_gbar_prior_center"] = anchors[Float64(entry["t_k"])]
        entry["mu_gbar_prior_sd"] = mu_gbar_sd
    end
    return out
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    smoke = "--smoke" in ARGS
    origins = nothing
    mu_gbar_sd = 0.2
    for (i, a) in enumerate(ARGS)
        a == "--origins" && (global origins = parse_origins(ARGS[i + 1]))
        a == "--mu-gbar-sd" && (global mu_gbar_sd = parse(Float64, ARGS[i + 1]))
    end
    countries = [c for c in ARGS if haskey(COUNTRY_CFG, c)]
    isempty(countries) && (countries = ["JPN", "THA"])
    for c in countries
        out = run_walkforward_m10(c; smoke, origins, mu_gbar_sd)
        suffix = smoke ? "_smoke" : ""
        path = joinpath(@__DIR__, "output", "M10_walkforward_$(c)$(suffix).json")
        mkpath(dirname(path))
        write(path, JSON3.write(out))
        println("保存: $path")
    end
end
