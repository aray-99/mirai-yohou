# 双子実験 E1(SPEC §13)
#
# 手順: 真値ラン(変動国・固定シード・45年)→ 合成観測 → E1a(状態推定のみ)
# → E1b(theta_sig 状態拡大)→ 対照の自由ラン → 凍結しきい値
# (test/acceptance_thresholds.toml)による合格判定。
# 判定手続きの細部は DECISIONS #0009/#0010 で事前確定済み。
#
# 実行: julia --project=experiments -t auto experiments/E1_twin.jl

using MiraiYohou
using Distributions
using Plots
using Random
using Statistics
using TOML
using Dates

const TRUTH_SEED = 4501
const OBS_SEED = 4502
const ENS_SEED = 4503
const FREE_SEED = 4504
const E1B_SEED = 4505
const N_ENS = 100
const T1 = 45.0
const INIT_SD = 0.3          # 初期アンサンブルの摂動 sd(§13 手順3)
const THRESH = TOML.parsefile(joinpath(dirname(@__DIR__), "test",
                                       "acceptance_thresholds.toml"))

# --- 評価関数(手続きは DECISIONS #0010) ---

rmse(est::AbstractVector, truth::AbstractVector) =
    sqrt(mean(abs2, est .- truth))

"アンサンブル平均の時系列(1変数)"
ens_mean(X, ix) = vec(mean(@view(X[ix, :, :]); dims = 2))

function coverage_sigma_s(X, truth_sig; ci = 0.90)
    lo, hi = (1 - ci) / 2, 1 - (1 - ci) / 2
    nt = size(X, 2)
    hits = 0
    for k in 1:nt
        vals = exp.(@view X[IX_SIG, k, :])
        hits += quantile(vals, lo) <= truth_sig[k] <= quantile(vals, hi)
    end
    return hits / nt
end

"プールした事前ランクの平坦性カイ二乗検定(10ビン、#0010)"
function rank_histogram_chi2(ranks::Dict{Symbol,Vector{Int}}, N::Int)
    pooled = reduce(vcat, values(ranks))
    nbins = 10
    counts = zeros(Int, nbins)
    for r in pooled
        counts[clamp(ceil(Int, r * nbins / (N + 1)), 1, nbins)] += 1
    end
    e = length(pooled) / nbins
    chi2 = sum((c - e)^2 / e for c in counts)
    crit = quantile(Chisq(nbins - 1), 1 - THRESH["E1"]["rank_histogram_alpha"])
    return (chi2 = chi2, critical = crit, counts = counts)
end

# --- 実験本体 ---

function main()
    figdir = joinpath(@__DIR__, "figures")
    mkpath(figdir)
    params = build_params(:volatile)
    println("=== E1 twin experiment: ", dimensionless_numbers(params))

    # 1. 真値ラン(§13 手順1)
    truth = simulate_sde(params; seed = TRUTH_SEED, t1 = T1)
    event_times = [e.t for e in truth.jumps]
    ms = [e.m for e in truth.jumps]
    println("truth: $(length(event_times)) jumps, max m = $(maximum(ms))")

    # 2. 合成観測(§13 手順2)
    obs = synthesize_observations(truth.traj, standard_observations(params);
                                  rng = Xoshiro(OBS_SEED))
    println("synthesized $(length(obs)) observations")

    # 3. 初期アンサンブル(真値から sd 0.3 の摂動、§13 手順3)
    rng_ens = Xoshiro(ENS_SEED)
    E0 = params.x0 .+ INIT_SD .* randn(rng_ens, N_STATE, N_ENS)
    postprocess_analysis!(E0)

    # E1a: 状態推定のみ(パラメータ真値固定)
    cfg = AssimConfig(t1 = T1)
    e1a = run_assimilation(params, copy(E0), obs, event_times;
                           cfg, seed = ENS_SEED)
    println("E1a done: $(e1a.nresample) resamples, min ESS = $(minimum(e1a.ess))")

    # 5. 対照: 同化オフの自由ラン(同じ初期アンサンブル)
    Xfree = free_ensemble(params, E0; cfg, seed = FREE_SEED)

    # --- E1a 合格判定(§13。しきい値は凍結 TOML) ---
    results = Dict{String,Any}()
    pass = Dict{String,Bool}()

    # (i) 観測に近い状態の RMSE 比 < 0.5
    near_obs = [(:xi_k, IX_K), (:xi_g, IX_G), (:xi_tau, IX_TAU), (:xi_p, IX_PP)]
    for (name, ix) in near_obs
        tr = vec(truth.traj.X[ix, :])
        ratio = rmse(ens_mean(e1a.X, ix), tr) / rmse(ens_mean(Xfree, ix), tr)
        results["rmse_ratio_$name"] = ratio
        pass["rmse_$name"] = ratio < THRESH["E1"]["rmse_ratio_max"]
    end

    # (ii) 隠れ状態 sigma_s: 時間相関 > 0.6、90% 区間被覆率 80〜98%
    truth_sig = exp.(vec(truth.traj.X[IX_SIG, :]))
    est_sig = vec(mean(exp.(e1a.X[IX_SIG, :, :]); dims = 2))
    results["sigma_s_correlation"] = cor(est_sig, truth_sig)
    pass["sigma_s_correlation"] =
        results["sigma_s_correlation"] > THRESH["E1"]["sigma_s_time_correlation_min"]
    results["sigma_s_coverage"] =
        coverage_sigma_s(e1a.X, truth_sig; ci = THRESH["E1"]["credible_interval"])
    pass["sigma_s_coverage"] =
        THRESH["E1"]["coverage_min"] <= results["sigma_s_coverage"] <=
        THRESH["E1"]["coverage_max"]

    # (iii) 塑性変数 tauA: 最大ジャンプ後 10 年窓で自由ランより誤差減少(#0010)
    tstar = event_times[argmax(ms)]
    win = findall(t -> tstar <= t <= min(tstar + THRESH["E1"]["tauA_window_years"], T1),
                  truth.traj.t)
    tauA_true = vec(truth.traj.X[IX_TAUA, :])
    err_assim = mean(abs.(ens_mean(e1a.X, IX_TAUA)[win] .- tauA_true[win]))
    err_free = mean(abs.(ens_mean(Xfree, IX_TAUA)[win] .- tauA_true[win]))
    results["tauA_err_assim"] = err_assim
    results["tauA_err_free"] = err_free
    pass["tauA_updated"] = err_assim < err_free

    # (iv) ランクヒストグラム平坦性(プール、#0010)
    rh = rank_histogram_chi2(e1a.ranks, N_ENS)
    results["rank_chi2"] = rh.chi2
    results["rank_chi2_critical"] = rh.critical
    results["rank_counts"] = rh.counts
    pass["rank_histogram"] = rh.chi2 <= rh.critical

    # 4. E1b: theta_sig 状態拡大(§13 手順4。初期値は事前分布からサンプル)
    prior = l3_priors().theta_sig
    E0b = vcat(copy(E0), reshape(log.(rand(Xoshiro(E1B_SEED), prior, N_ENS)), 1, :))
    e1b = run_assimilation(params, E0b, obs, event_times;
                           cfg, seed = E1B_SEED, augmented = true)
    theta_post = mean(exp.(e1b.X[end, end, :]))
    theta_true = params.l3.theta_sig
    results["theta_sig_posterior_mean"] = theta_post
    results["theta_sig_true"] = theta_true
    rel_err = abs(theta_post - theta_true) / theta_true
    results["theta_sig_relative_error"] = rel_err
    pass["e1b_theta_sig"] = rel_err <= THRESH["E1b"]["theta_sig_relative_error_max"]

    # --- 出力 ---
    println("\n=== acceptance results ===")
    for k in sort(collect(keys(pass)))
        println(rpad(k, 24), pass[k] ? "PASS" : "FAIL")
    end
    println("\n=== metrics ===")
    for k in sort(collect(keys(results)))
        v = results[k]
        v isa Vector || println(rpad(k, 28), v)
    end
    e1a_pass = all(v for (k, v) in pass if k != "e1b_theta_sig")
    e1b_pass = pass["e1b_theta_sig"]
    println("\nE1a: ", e1a_pass ? "PASS" : "FAIL")
    println("E1b: ", e1b_pass ? "PASS" : "FAIL")

    # 図: 隠れ状態と塑性変数の復元
    sha = strip(read(`git -C $(dirname(@__DIR__)) rev-parse HEAD`, String))
    panels = [(IX_SIG, "sigma_s (hidden)", x -> exp.(x)),
              (IX_TAUA, "tauA (plastic anchor)", identity),
              (IX_TAU, "tau", x -> sigmoid.(x)),
              (IX_G, "g", x -> sigmoid.(x))]
    plts = map(panels) do (ix, label, tf)
        plt = plot(; title = label, titlefontsize = 9,
                   legend = (ix == IX_SIG ? :topright : false))
        vals = tf(e1a.X[ix, :, :])
        lo = [quantile(@view(vals[k, :]), 0.05) for k in axes(vals, 1)]
        hi = [quantile(@view(vals[k, :]), 0.95) for k in axes(vals, 1)]
        med = vec(mean(vals; dims = 2))
        plot!(plt, e1a.t, med; ribbon = (med .- lo, hi .- med), fillalpha = 0.25,
              color = :steelblue, lw = 1.2, label = "assimilated")
        plot!(plt, truth.traj.t, tf(vec(truth.traj.X[ix, :]));
              color = :black, lw = 1.2, label = "truth")
        plot!(plt, e1a.t, vec(mean(tf(Xfree[ix, :, :]); dims = 2));
              color = :gray, ls = :dash, lw = 1.0, label = "free run")
        vline!(plt, [tstar]; color = :firebrick, alpha = 0.4, label = false)
        plt
    end
    fig = plot(plts...; layout = (2, 2), size = (1000, 700),
               plot_title = "E1a twin experiment (volatile, 45y, N=$N_ENS)")
    figpath = joinpath(figdir, "E1a_twin.png")
    savefig(fig, figpath)

    meta = Dict("experiment" => "E1_twin", "commit_sha" => sha,
                "seeds" => Dict("truth" => TRUTH_SEED, "obs" => OBS_SEED,
                                "ensemble" => ENS_SEED, "free" => FREE_SEED,
                                "e1b" => E1B_SEED),
                "parameter_set" => "volatile", "N" => N_ENS, "t1_years" => T1,
                "generated_at" => Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"),
                "results" => Dict(k => v for (k, v) in results if !(v isa Vector)),
                "pass" => pass, "E1a_pass" => e1a_pass, "E1b_pass" => e1b_pass)
    open(figpath * ".meta.json", "w") do io
        # 依存追加を避けた素朴な JSON 出力
        function j(x)
            x isa Dict ? "{" * join(["\"$k\": $(j(v))" for (k, v) in x], ", ") * "}" :
            x isa Bool ? string(x) :
            x isa Number ? string(x) : "\"$x\""
        end
        write(io, j(meta))
    end
    println("wrote $figpath (+ sidecar metadata)")
    return (; e1a_pass, e1b_pass, results)
end

main()
