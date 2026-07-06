# E0 スモークテスト(SPEC §12 M1)
#
# 安定国セットでドリフトのみの決定論 ODE を50年走らせ、主要変数の
# 時系列プロットを experiments/figures/ に出力する。参考として変動国も重ねる。
# 成果物には来歴メタデータ(コミット SHA・シード・パラメータセット・生成日時)を
# サイドカー JSON として付す(§0.5.6)。
#
# 実行: julia --project=experiments experiments/E0_smoke.jl

using MiraiYohou
using Plots
using Dates

function run_smoke()
    figdir = joinpath(@__DIR__, "figures")
    mkpath(figdir)

    trajs = Dict{Symbol,Trajectory}()
    for regime in (:stable, :volatile)
        params = build_params(regime)
        println("=== $regime: ", dimensionless_numbers(params))
        trajs[regime] = simulate_ode(params; t1 = 50.0)
    end

    # 主要変数(自然座標に戻して描画)
    panels = [
        (IX_G, "g (Gini)", true),
        (IX_TAU, "tau (institutional trust)", true),
        (IX_SIG, "sigma_s (social stress)", false),   # log 座標 → exp
        (IX_PP, "p (polarization)", true),
        (IX_T, "T (tech frontier)", false),
        (IX_PHI, "phi (diffusion)", true),
    ]
    plts = map(panels) do (ix, label, islogit)
        plt = plot(; title = label, titlefontsize = 9, legend = (ix == IX_G ? :topleft : false))
        for (regime, color) in ((:stable, :steelblue), (:volatile, :firebrick))
            tr = trajs[regime]
            vals = islogit ? sigmoid.(tr.X[ix, :]) : exp.(tr.X[ix, :])
            plot!(plt, tr.t, vals; label = String(regime), color, lw = 1.5)
        end
        plt
    end
    fig = plot(plts...; layout = (3, 2), size = (900, 900),
               plot_title = "E0 smoke: drift-only ODE, 50y (M1)")
    figpath = joinpath(figdir, "E0_smoke.png")
    savefig(fig, figpath)

    # 来歴メタデータ(§0.5.6)。決定論ランのためシードなし。
    sha = strip(read(`git -C $(dirname(@__DIR__)) rev-parse HEAD`, String))
    meta = """
    {
      "experiment": "E0_smoke",
      "commit_sha": "$sha",
      "seed": null,
      "parameter_sets": ["stable", "volatile"],
      "dt": 0.01,
      "t1_years": 50.0,
      "generated_at": "$(Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))"
    }
    """
    write(figpath * ".meta.json", meta)
    println("wrote $figpath (+ sidecar metadata)")
    return trajs
end

run_smoke()
