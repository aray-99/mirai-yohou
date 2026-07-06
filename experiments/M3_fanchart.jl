# M3 ファンチャート(SPEC §12 M3: 予報円プロット出力)
#
# 両レジームで N=100 のアンサンブルを50年走らせ、主要変数の
# 分位ファン(5–95% / 25–75% / 中央値)を出力する。
# 来歴メタデータは §0.5.6 に従いサイドカー JSON。
#
# 実行: julia --project=experiments -t auto experiments/M3_fanchart.jl

using MiraiYohou
using Plots
using Statistics
using Dates

const SEED = 20260707

function fanchart!(plt, t, vals; color)
    # vals: 時刻 × メンバー
    q(p) = [quantile(@view(vals[k, :]), p) for k in eachindex(t)]
    lo5, lo25, med, hi75, hi95 = q(0.05), q(0.25), q(0.5), q(0.75), q(0.95)
    plot!(plt, t, med; ribbon = (med .- lo5, hi95 .- med),
          fillalpha = 0.15, color, lw = 0, label = false)
    plot!(plt, t, med; ribbon = (med .- lo25, hi75 .- med),
          fillalpha = 0.3, color, lw = 1.5, label = false)
    return plt
end

function run_fancharts()
    figdir = joinpath(@__DIR__, "figures")
    mkpath(figdir)

    panels = [
        (IX_G, "g (Gini)", true),
        (IX_TAU, "tau (institutional trust)", true),
        (IX_SIG, "sigma_s (social stress)", false),
        (IX_PP, "p (polarization)", true),
        (IX_K, "k (capital per capita)", false),
        (IX_PHI, "phi (diffusion)", true),
    ]

    njumps = Dict{Symbol,Float64}()
    for regime in (:stable, :volatile)
        params = build_params(regime)
        println("=== $regime: ", dimensionless_numbers(params))
        ens = simulate_ensemble(params; N = 100, seed = SEED, t1 = 50.0)
        njumps[regime] = sum(length, ens.jumps) / 100

        plts = map(panels) do (ix, label, islogit)
            plt = plot(; title = label, titlefontsize = 9, legend = false)
            vals = islogit ? sigmoid.(ens.X[ix, :, :]) : exp.(ens.X[ix, :, :])
            fanchart!(plt, ens.t, vals;
                      color = regime === :stable ? :steelblue : :firebrick)
        end
        fig = plot(plts...; layout = (3, 2), size = (900, 900),
                   plot_title = "M3 fan chart: $regime, N=100, 50y " *
                                "(mean jumps/run = $(njumps[regime]))")
        figpath = joinpath(figdir, "M3_fanchart_$(regime).png")
        savefig(fig, figpath)

        sha = strip(read(`git -C $(dirname(@__DIR__)) rev-parse HEAD`, String))
        write(figpath * ".meta.json", """
        {
          "experiment": "M3_fanchart",
          "regime": "$regime",
          "commit_sha": "$sha",
          "seed": $SEED,
          "ensemble_size": 100,
          "dt": 0.01,
          "t1_years": 50.0,
          "mean_jumps_per_run": $(njumps[regime]),
          "generated_at": "$(Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))"
        }
        """)
        println("wrote $figpath (+ sidecar metadata)")
    end
end

run_fancharts()
