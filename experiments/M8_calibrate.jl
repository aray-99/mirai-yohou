# M8 L2 較正: 反復 EnKF(EKI 方式、DECISIONS #0030-5/#0031)
#
# 較正対象 θ = {eta_g, delta_sig, lam0, mu_gbar, mu_p, c_v0, nu}(7個)。
# 変換空間(正値パラメータは log、mu_* と c_v0 は恒等)で J 本の
# パラメータアンサンブルを事前分布から引き、各 θ_j で較正期間の同化ランを
# 実行して G(θ_j) を評価し、EKI 更新 θ ← θ + C_θG (C_GG + Γ)^{-1}(0 − G_j)
# を反復する。
#
# G(θ) の定義(標準化イノベーション。目標 y* = 0):
#   - 観測系列ごと: mean_t [ (y_t − mean_j H(x_j,t)) / sd_t ](較正窓)
#   - カウント: (総観測数 − ν·Λ̄_total) / √(ν·Λ̄_total)(カバレッジ窓のみ)
# 較正期間のデータのみ使用(#0030-4)。結果は JSON(来歴付き)に保存。
#
# 実行: julia --project=experiments experiments/M8_calibrate.jl [JPN|THA]
#       [--J 16] [--iters 3] [--N 50]

using JSON3
include(joinpath(@__DIR__, "M8_hindcast.jl"))
using MiraiYohou: intensity

# ---- パラメータ変換(#0031: 正値は log) ----
# ν は較正対象に含めない: ポアソン尤度の最尤解 ν* = 総観測数 / Λ̄ で
# 各 θ ごとに解析的にプロファイル化する(#0031-1 の実装細目。λ の水準は
# lam0・θ_sig・σ_s 経路で桁ごと動くため、ν を EKI で回すと縮退する)
const CAL_PARAMS = [
    (name = :eta_g,     trans = :log),
    (name = :delta_sig, trans = :log),
    (name = :lam0,      trans = :log),
    (name = :mu_gbar,   trans = :id),
    (name = :mu_p,      trans = :id),
    (name = :c_v0,      trans = :id),
]

to_eta(vals) = [CAL_PARAMS[k].trans === :log ? log(vals[k]) : vals[k]
                for k in eachindex(CAL_PARAMS)]
from_eta(eta) = [CAL_PARAMS[k].trans === :log ? exp(eta[k]) : eta[k]
                 for k in eachindex(CAL_PARAMS)]

"θ ベクトル → ModelParameters"
function params_from_theta(regime, theta)
    kw = Dict(CAL_PARAMS[k].name => theta[k] for k in eachindex(CAL_PARAMS))
    return build_params(regime; kw...)
end

"""
較正の前処理(θ 非依存の入力を1回だけ構築)。

`window` 省略時は国既定の較正窓(#0030-4、M8 の従来動作)。M9(#0052)の
expanding walk-forward では `(窓開始, オリジン t_k)` を明示的に渡し、
オリジン以前のデータのみで較正する(検証区間は較正・同化とも未参照)。
`include_theta_sig` は theta_sig 拡大の適用可否(既定 = M8 の国条件 #0049-2。
M9 はデータ規則 #0054 の判定結果を渡す)。
"""
function calib_inputs(country; N, seed, window = COUNTRY_CFG[country].calib,
                      include_theta_sig::Bool = country != "JPN")
    ccfg = COUNTRY_CFG[country]
    params0 = build_params(ccfg.regime)
    recs = build_observations(country, params0; t1 = window[2])
    cfg = AssimConfig(t0 = 0.0, t1 = window[2], smoother_lag = 0.0,
                      smoother_vars = [IX_G, IX_TAU, IX_SIG, IX_PP],  # #0036: tauA除外
                      tauA_pseudo_sd_mult = 3.0,                     # #0036
                      analysis_masked_vars = [IX_TAUA],               # #0040-(α)
                      analysis_unmask_names = [:tau, :tauA_pseudo],   # #0040-(α)
                      rtps_alpha = 0.85,                              # #0040-(β)
                      obs_spread_floor_frac = 0.5)                    # #0043
    E0 = initial_ensemble(country, params0, recs; N, seed = seed + 1)
    obs_counts = build_obs_counts(country, cfg)
    event_times = filter(t -> t < cfg.t1,
                         build_forced_jumps(country; calib_window = window))
    nu0 = init_nu(obs_counts, cfg, params0, window)
    return (; country, ccfg, window, recs, cfg, E0, obs_counts, event_times, nu0,
            include_theta_sig)
end

"""
G(θ): 較正窓の系列別平均標準化イノベーション + カウント整合。
h_logy が使う θ_T/θ_phi/α は較正対象外なので観測レコードは θ 非依存。
発散・レジーム違反ランは nothing(EKI 更新から除外)。
"""
function forward_G(inp, theta; N, seed)
    ccfg, recs, cfg, window = inp.ccfg, inp.recs, inp.cfg, inp.window
    params = try
        fit_exogenous(params_from_theta(ccfg.regime, theta), recs, window)
    catch                                  # De レジーム条件等の違反
        return nothing
    end
    nu = inp.nu0    # 重み付け用の ν(モーメント初期値。プロファイル解は下で計算)
    obs_counts = inp.obs_counts
    # 拡大パラメータ(#0046)。M8_hindcast.jl の検証ランと同じ機構・数値
    # (M8_AUG_SETTINGS)を使い、較正と検証の力学を一致させる。初期アンサンブル
    # の乱数は θ_j 間で共通(seed のみに依存)にして CRN 性を保つ(呼び出し側
    # は j に依らず同一の seed = seed_it を渡す設計、上記コメント参照)。
    aug = build_m8_augmented_params(params, inp.country;
                                    include_theta_sig = inp.include_theta_sig)
    E0 = augment_ensemble(inp.E0, aug; rng = Xoshiro(seed))
    res = try
        run_assimilation(params, E0, recs, inp.event_times;
                         cfg, seed, obs_counts, count_scale = nu,
                         count_temper = 1 / nu, augmented_params = aug)
    catch
        return nothing
    end
    # 系列別イノベーション(較正窓のみ)。較正対象 θ が力学的に影響する
    # 系列に限定する(P 等は外生・非依存で、含めると EKI を歪めるだけ)
    G = Float64[]
    for s in REAL_SERIES
        Symbol(s.var) in (:g_swiid, :p, :v, :tau) || continue
        num = Float64[]
        for r in recs
            r.spec.name === Symbol(s.var) || continue
            window[1] <= r.t < window[2] || continue
            k = clamp(searchsortedlast(res.t, r.t), 1, length(res.t))
            m = mean(r.spec.h(view(res.X, 1:N_STATE, k, j)) for j in 1:N)
            push!(num, (r.value - m) / r.spec.sd)
        end
        isempty(num) || push!(G, mean(num))
    end
    # カウント整合(年別パターン): ν はプロファイル最尤解
    # ν* = Σ N_y / Σ Λ_y で吸収し、年別の標準化残差でタイミング情報を取る
    w = cfg.event_window
    wsteps = max(1, round(Int, w / cfg.dt))
    Ny = Dict{Int, Int}(); Ly = Dict{Int, Float64}()
    for (k, c) in enumerate(obs_counts)
        c >= 0 || continue
        lo = cfg.t0 + (k - 1) * w
        window[1] <= lo < window[2] || continue
        yr = floor(Int, lo)
        Ny[yr] = get(Ny, yr, 0) + c
        acc = 0.0
        i0 = (k - 1) * wsteps + 1
        for i in i0:min(i0 + wsteps - 1, length(res.t))
            acc += mean(intensity(view(res.X, 1:N_STATE, i, j), params)
                        for j in 1:N) * cfg.dt
        end
        Ly[yr] = get(Ly, yr, 0.0) + acc
    end
    nu_star = nu
    if !isempty(Ny)
        Ltot = sum(values(Ly))
        nu_star = max(sum(values(Ny)) / max(Ltot, 1e-8), 1.0)
        for yr in sort(collect(keys(Ny)))
            e = nu_star * Ly[yr]
            push!(G, (Ny[yr] - e) / sqrt(max(e, 1.0)))
        end
    end
    return (; G, nu_star)
end

"EKI 1反復: η アンサンブル(d×J)と G(m×J)からカルマン型更新"
function eki_update(H::Matrix{Float64}, G::Matrix{Float64}; gamma::Float64 = 1.0,
                    rng)
    J = size(H, 2); m = size(G, 1)
    Hm = H .- mean(H; dims = 2)
    Gm = G .- mean(G; dims = 2)
    C_hg = (Hm * Gm') ./ (J - 1)
    C_gg = (Gm * Gm') ./ (J - 1)
    K = C_hg / (C_gg + gamma * I(m))
    # 摂動観測: 目標 0 + √γ ノイズ
    return H .+ K * (sqrt(gamma) .* randn(rng, m, J) .- G)
end

using LinearAlgebra: I

"""
    calibrate(country; J, iters, N, seed, window, prior_center, prior_sd, save) -> (; theta_hat, nu_star, out)

EKI 較正(DECISIONS #0030-5/-6/-31)。`window` 省略時は国既定の較正窓
(M8 の従来動作)。M9(#0052)の expanding walk-forward では:

- `window = (窓開始, オリジン t_k)` を渡し、オリジン以前のデータのみで較正。
- `prior_center`: 初回オリジンは `nothing`(既定値からモーメント初期化、
  従来動作と同じ)。以後のオリジンは前オリジンの較正値の Dict
  (`CAL_PARAMS` の名前をキーとする)を渡して warm-start する。
- `prior_sd`: 変換空間の事前 sd。毎オリジンとも既定値 0.5 に復元する
  (#0052 の「prior sd は毎回初期値に復元」)。
- `prior_sd_override`(#0062、M10): パラメータ名 → sd の Dict。指定された
  パラメータのみ `prior_sd` の代わりにこの値を使う(他は `prior_sd` のまま)。
  `nothing`(既定)なら従来どおり全パラメータ一律 `prior_sd`(M8/M9 の動作を
  変えない後方互換の既定)。
- `save`: `true`(既定)なら結果を `output/M8_calib_<country>.json` に保存
  (M8 CLI の従来動作)。M9 は複数オリジンで上書きされるのを避けるため
  `false` を渡し、戻り値の `theta_hat`/`nu_star` を直接使う。
"""
function calibrate(country::String; J::Int = 16, iters::Int = 3, N::Int = 50,
                   seed::Integer = 20260710,
                   window = COUNTRY_CFG[country].calib,
                   prior_center::Union{Nothing,Dict} = nothing,
                   prior_sd::Float64 = 0.5,
                   prior_sd_override::Union{Nothing,AbstractDict} = nothing,
                   save::Bool = true,
                   include_theta_sig::Bool = country != "JPN")
    ccfg = COUNTRY_CFG[country]
    rng = Xoshiro(seed)
    params0 = build_params(ccfg.regime)
    recs0 = build_observations(country, params0; t1 = window[2])
    cfg0 = AssimConfig(t0 = 0.0, t1 = window[2])
    center = prior_center !== nothing ? prior_center : Dict(
        :eta_g => params0.l2.eta_g, :delta_sig => params0.l2.delta_sig,
        :lam0 => params0.l2.lam0, :mu_gbar => params0.l2.mu_gbar,
        :mu_p => params0.l2.mu_p, :c_v0 => init_c_v0(recs0, window))
    theta0 = [center[p.name] for p in CAL_PARAMS]
    eta0 = to_eta(theta0)
    # 事前: 変換空間で sd `prior_sd`(既定 0.5。§8.2 の LogNormal(·, 0.5) と
    # 同型)。mu_* / c_v0 は logit/log 空間の値なので加法。パラメータ別に
    # sd を変えたい場合は `prior_sd_override`(#0062、M10 の mu_gbar 用)。
    sd_vec = prior_sd_override === nothing ? fill(prior_sd, length(CAL_PARAMS)) :
        [Float64(get(prior_sd_override, CAL_PARAMS[k].name, prior_sd))
         for k in eachindex(CAL_PARAMS)]
    H = eta0 .+ sd_vec .* randn(rng, length(eta0), J)
    println("== $country EKI 較正: J=$J iters=$iters N=$N window=$window ==")
    println("  初期中心 θ = ", round.(theta0, digits = 3))
    inp = calib_inputs(country; N, seed, window, include_theta_sig)
    nu_star = inp.nu0
    for it in 1:iters
        Gcols = Vector{Any}(undef, J)
        # 共通乱数(全 θ_j で同一シード)で G の標本雑音を θ 間で相殺する
        Threads.@threads for j in 1:J
            Gcols[j] = forward_G(inp, from_eta(H[:, j]); N, seed = seed + it)
        end
        ok = findall(!isnothing, Gcols)
        length(ok) >= max(4, J ÷ 2) ||
            error("EKI iteration $it: 有効ラン $(length(ok))/$J が少なすぎます")
        m = length(Gcols[ok[1]].G)
        Gm = zeros(m, length(ok))
        for (c, j) in enumerate(ok); Gm[:, c] = Gcols[j].G; end
        nu_star = mean(Gcols[j].nu_star for j in ok)
        Hok = H[:, ok]
        misfit = mean(sqrt.(sum(abs2, Gm; dims = 1))[:])
        # LM 型の適応減衰(Iglesias 流): γ を平均二乗ミスフィットに比例させ、
        # 大ミスフィット時の過大なカルマンステップを抑える
        gamma = max(1.0, mean(sum(abs2, Gm; dims = 1)) / m)
        println("  iter $it: 有効 $(length(ok))/$J  平均ミスフィット ‖G‖ = ",
                round(misfit, digits = 2), "  γ = ", round(gamma, digits = 1))
        Hnew = eki_update(Hok, Gm; rng, gamma)
        # J 本に補充(欠損は更新後アンサンブルからリサンプル)
        H = Hnew[:, mod1.(1:J, size(Hnew, 2))]
    end
    eta_hat = vec(mean(H; dims = 2))
    theta_hat = from_eta(eta_hat)
    println("  較正結果 θ̂ = ", round.(theta_hat, digits = 3),
            "  ν* = ", round(nu_star, digits = 1))
    out = Dict(
        "country" => country,
        "params" => merge(
            Dict(string(CAL_PARAMS[k].name) => theta_hat[k]
                 for k in eachindex(CAL_PARAMS)),
            Dict("nu" => nu_star)),
        "ensemble_final" => [from_eta(H[:, j]) for j in 1:J],
        "J" => J, "iters" => iters, "N" => N, "seed" => seed,
        "calib_window" => collect(window),
        "commit" => strip(read(`git -C $(dirname(@__DIR__)) rev-parse HEAD`, String)),
        "generated_at" => Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))
    if save
        path = joinpath(@__DIR__, "output", "M8_calib_$(country).json")
        mkpath(dirname(path))
        write(path, JSON3.write(out))
        println("  保存: $path")
    end
    return (; theta_hat, nu_star, out)
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    J = 16; iters = 3; N = 50
    for (i, a) in enumerate(ARGS)
        a == "--J" && (global J = parse(Int, ARGS[i + 1]))
        a == "--iters" && (global iters = parse(Int, ARGS[i + 1]))
        a == "--N" && (global N = parse(Int, ARGS[i + 1]))
    end
    countries = [c for c in ARGS if haskey(COUNTRY_CFG, c)]
    isempty(countries) && (countries = ["THA", "JPN"])
    for c in countries
        calibrate(c; J, iters, N)
    end
end
