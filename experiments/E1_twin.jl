# 双子実験 E1(SPEC §13 プロトコル v1.1、DECISIONS #0015)
#
# 事前凍結した未使用の5 truth シード(acceptance_thresholds.toml の truth_seeds)
# で真値ラン→合成観測→E1a/E1b/自由ランを反復し、数値指標は中央値・
# 二値指標は多数決で判定する。数値しきい値は v1.0 から不変。
# v1.0(単一シード4501)の結果は DECISIONS #0012/#0014 に恒久記録。
#
# 実行: julia --project=experiments -t auto experiments/E1_twin.jl

using MiraiYohou
using Distributions
using Plots
using Random
using Statistics
using TOML
using Dates
using Logging

const N_ENS = 100
const T1 = 45.0
const INIT_SD = 0.3          # 初期アンサンブルの摂動 sd(§13 手順3)
const THRESH = TOML.parsefile(joinpath(dirname(@__DIR__), "test",
                                       "acceptance_thresholds.toml"))
const SEEDS = Int.(THRESH["E1"]["protocol"]["truth_seeds"])
const MAJORITY = THRESH["E1"]["protocol"]["majority_min"]

rmse(est, truth) = sqrt(mean(abs2, est .- truth))
ens_mean(X, ix) = vec(mean(@view(X[ix, :, :]); dims = 2))

"1シード分の E1a/E1b 指標(シード規約: truth=s, obs=s+1, ens=s+2, free=s+3, e1b=s+4)"
function run_seed(s::Int, params::ModelParameters, cfg::AssimConfig)
    truth = simulate_sde(params; seed = s, t1 = T1)
    event_times = [e.t for e in truth.jumps]
    ms = [e.m for e in truth.jumps]
    obs = synthesize_observations(truth.traj, standard_observations(params);
                                  rng = Xoshiro(s + 1))
    E0 = params.x0 .+ INIT_SD .* randn(Xoshiro(s + 2), N_STATE, N_ENS)
    postprocess_analysis!(E0)

    e1a = run_assimilation(params, copy(E0), obs, event_times; cfg, seed = s + 2)
    Xfree = free_ensemble(params, E0; cfg, seed = s + 3)

    m = Dict{Symbol,Any}()
    for (name, ix) in [(:k, IX_K), (:g, IX_G), (:tau, IX_TAU)]
        tr = vec(truth.traj.X[ix, :])
        m[Symbol(:rmse_, name)] =
            rmse(ens_mean(e1a.X, ix), tr) / rmse(ens_mean(Xfree, ix), tr)
    end

    truth_sig = exp.(vec(truth.traj.X[IX_SIG, :]))
    m[:corr_sig] = cor(vec(mean(exp.(e1a.X[IX_SIG, :, :]); dims = 2)), truth_sig)
    lo, hi = (1 - THRESH["E1"]["credible_interval"]) / 2,
             1 - (1 - THRESH["E1"]["credible_interval"]) / 2
    m[:cover_sig] = mean(1:length(truth_sig)) do k
        vals = exp.(@view e1a.X[IX_SIG, k, :])
        quantile(vals, lo) <= truth_sig[k] <= quantile(vals, hi)
    end

    # 較正(v1.1): 全観測の事前90%区間被覆(プールしたランクから計算)
    pooled = reduce(vcat, values(e1a.ranks))
    rlo, rhi = lo * (N_ENS + 1), hi * (N_ENS + 1)
    m[:cover_obs] = mean(r -> rlo <= r <= rhi, pooled)

    # tauA: 最大ジャンプ後 10 年窓(二値、#0010)
    tstar = event_times[argmax(ms)]
    win = findall(t -> tstar <= t <= min(tstar + THRESH["E1"]["tauA_window_years"], T1),
                  truth.traj.t)
    tauA_true = vec(truth.traj.X[IX_TAUA, :])
    m[:tauA_assim] = mean(abs.(ens_mean(e1a.X, IX_TAUA)[win] .- tauA_true[win]))
    m[:tauA_free] = mean(abs.(ens_mean(Xfree, IX_TAUA)[win] .- tauA_true[win]))
    m[:tauA_pass] = m[:tauA_assim] < m[:tauA_free]

    # E1b: theta_sig 状態拡大(§13 手順4)
    prior = l3_priors().theta_sig
    E0b = vcat(copy(E0), reshape(log.(rand(Xoshiro(s + 4), prior, N_ENS)), 1, :))
    e1b = run_assimilation(params, E0b, obs, event_times; cfg, seed = s + 4,
                           augmented = true)
    theta_post = mean(exp.(e1b.X[end, end, :]))
    m[:theta_rel_err] = abs(theta_post - params.l3.theta_sig) / params.l3.theta_sig
    m[:theta_post] = theta_post
    m[:njumps] = length(event_times)
    return m, (; truth, e1a, Xfree, tstar)
end

function main()
    figdir = joinpath(@__DIR__, "figures")
    mkpath(figdir)
    params = build_params(:volatile)
    cfg = AssimConfig(t1 = T1)   # 採用方式: RTPS α=0.7(#0013)
    println("=== E1 twin experiment v1.1: ", dimensionless_numbers(params))
    println("seeds = ", SEEDS)

    all_m = Dict{Symbol,Any}[]
    detail1 = nothing
    for (i, s) in enumerate(SEEDS)
        m, detail = with_logger(NullLogger()) do
            run_seed(s, params, cfg)
        end
        i == 1 && (detail1 = detail)
        push!(all_m, m)
        println("seed $s: ", join(["$k=$(v isa Bool ? v : round(v, digits=3))"
                                   for (k, v) in sort(collect(m))], " "))
    end

    med(key) = median([m[key] for m in all_m])
    pass = Dict{String,Bool}()
    results = Dict{String,Any}()
    for name in (:k, :g, :tau)
        key = Symbol(:rmse_, name)
        results["rmse_ratio_$(name)_median"] = med(key)
        pass["rmse_$name"] = med(key) < THRESH["E1"]["rmse_ratio_max"]
    end
    results["sigma_s_correlation_median"] = med(:corr_sig)
    pass["sigma_s_correlation"] =
        med(:corr_sig) > THRESH["E1"]["sigma_s_time_correlation_min"]
    results["sigma_s_coverage_median"] = med(:cover_sig)
    pass["sigma_s_coverage"] =
        THRESH["E1"]["coverage_min"] <= med(:cover_sig) <= THRESH["E1"]["coverage_max"]
    results["obs_coverage_median"] = med(:cover_obs)
    pass["obs_calibration"] =
        THRESH["E1"]["coverage_min"] <= med(:cover_obs) <= THRESH["E1"]["coverage_max"]
    ntauA = count(m -> m[:tauA_pass], all_m)
    results["tauA_pass_count"] = ntauA
    pass["tauA_updated"] = ntauA >= MAJORITY
    results["theta_sig_relative_error_median"] = med(:theta_rel_err)
    pass["e1b_theta_sig"] =
        med(:theta_rel_err) <= THRESH["E1b"]["theta_sig_relative_error_max"]

    println("\n=== acceptance (v1.1: median / majority over $(length(SEEDS)) seeds) ===")
    for k in sort(collect(keys(pass)))
        println(rpad(k, 24), pass[k] ? "PASS" : "FAIL")
    end
    println("\n=== aggregated metrics ===")
    for k in sort(collect(keys(results)))
        println(rpad(k, 34), results[k])
    end
    e1a_pass = all(v for (k, v) in pass if k != "e1b_theta_sig")
    e1b_pass = pass["e1b_theta_sig"]
    println("\nE1a: ", e1a_pass ? "PASS" : "FAIL")
    println("E1b: ", e1b_pass ? "PASS" : "FAIL")

    # 図(先頭シードの代表4変数)+ 来歴メタデータ(§0.5.6)
    sha = strip(read(`git -C $(dirname(@__DIR__)) rev-parse HEAD`, String))
    truth, e1a, Xfree, tstar = detail1.truth, detail1.e1a, detail1.Xfree, detail1.tstar
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
        med_t = vec(mean(vals; dims = 2))
        plot!(plt, e1a.t, med_t; ribbon = (med_t .- lo, hi .- med_t),
              fillalpha = 0.25, color = :steelblue, lw = 1.2, label = "assimilated")
        plot!(plt, truth.traj.t, tf(vec(truth.traj.X[ix, :]));
              color = :black, lw = 1.2, label = "truth")
        plot!(plt, e1a.t, vec(mean(tf(Xfree[ix, :, :]); dims = 2));
              color = :gray, ls = :dash, lw = 1.0, label = "free run")
        vline!(plt, [tstar]; color = :firebrick, alpha = 0.4, label = false)
        plt
    end
    fig = plot(plts...; layout = (2, 2), size = (1000, 700),
               plot_title = "E1a twin v1.1 (volatile, 45y, N=$N_ENS, seed $(SEEDS[1]))")
    figpath = joinpath(figdir, "E1a_twin.png")
    savefig(fig, figpath)

    meta = Dict("experiment" => "E1_twin_v1.1", "commit_sha" => sha,
                "truth_seeds" => "$(SEEDS)", "parameter_set" => "volatile",
                "N" => N_ENS, "t1_years" => T1,
                "inflation" => "rtps alpha=0.7 (#0013)",
                "generated_at" => Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"),
                "results" => results, "pass" => pass,
                "E1a_pass" => e1a_pass, "E1b_pass" => e1b_pass)
    open(figpath * ".meta.json", "w") do io
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
