# M10 prior 規則(DECISIONS #0062)の単体テスト:
#   (1) gbar_anchor の窓内 logit 平均計算
#   (2) prior_sd_override が calibrate の事前アンサンブル摂動に正しく伝わること
#
# experiments/M8_hindcast.jl 系列と同じ理由(パッケージテスト環境に無い
# Dates/JSON3/TOML 等に依存)で package の test/runtests.jl には組み込まず
# experiments 側の簡易テストとして分離する(test_m8_augmentation.jl と同じ流儀)。
#
# 実行: julia --project=experiments experiments/test_m10_prior.jl

using Test
using Statistics
using Random
using MiraiYohou
using MiraiYohou: logit

include(joinpath(@__DIR__, "M10_walkforward.jl"))

@testset "gbar_anchor 窓内 logit 平均(#0062 規則1)" begin
    country = "THA"
    window = (20.0, 28.0)
    anchor = gbar_anchor(country, window)

    # 独立に CSV を読み、build_observations と同一の logit 変換・時間座標
    # (t = year - 1990 + 0.5)で窓内平均を手計算して照合する
    years = Int[]; values = Float64[]
    for (i, line) in enumerate(eachline(joinpath(@__DIR__, "data", "raw", "THA_g_swiid.csv")))
        i == 1 && continue
        y, v = split(line, ",")
        push!(years, parse(Int, y)); push!(values, parse(Float64, v))
    end
    expected = mean(logit(v) for (y, v) in zip(years, values)
                    if window[1] <= (y - 1990.0 + 0.5) <= window[2])
    @test anchor ≈ expected atol=1e-9

    # 平衡的な範囲(g ∈ (0,1) の妥当域)にあること — logit 座標なので実数値だが
    # 極端(|anchor| > 10、§9.4 警告水準相当)ではないはず
    @test abs(anchor) < 10

    # 窓を変えると値が変わる(較正窓ごとに違う中心を与える設計であることの確認)
    anchor2 = gbar_anchor(country, (20.0, 33.0))
    @test anchor2 != anchor
end

@testset "prior_sd_override が calibrate の事前摂動に伝わる(#0062 規則3)" begin
    # iters=0: EKI 更新ループを回さず初期事前アンサンブル H の摂動のみを見る
    # (forward_G/同化ランを一切呼ばないので高速)。CAL_PARAMS の並びは
    # [eta_g(log), delta_sig(log), lam0(log), mu_gbar(id), mu_p(id), c_v0(id)]
    # で mu_gbar と mu_p はどちらも恒等変換なので prior_sd の効き方を単純比較できる。
    small_sd = 0.02
    default_sd = 0.5
    calib = calibrate("THA"; J = 400, iters = 0, N = 5, seed = 12345,
                      prior_sd = default_sd,
                      prior_sd_override = Dict(:mu_gbar => small_sd),
                      save = false)
    ens = calib.out["ensemble_final"]     # J 本の θ ベクトル(from_eta 適用済み)
    mu_gbar_vals = [e[4] for e in ens]
    mu_p_vals = [e[5] for e in ens]

    @test std(mu_gbar_vals) < 0.3 * std(mu_p_vals)
    @test isapprox(std(mu_gbar_vals), small_sd; rtol = 0.35)
    @test isapprox(std(mu_p_vals), default_sd; rtol = 0.35)

    # override 無し(nothing、既定)なら M8/M9 の従来動作(全パラメータ一律 sd)と
    # 同じ挙動 — mu_gbar と mu_p の広がりが同程度になることを確認
    calib_uniform = calibrate("THA"; J = 400, iters = 0, N = 5, seed = 12345,
                              prior_sd = default_sd, save = false)
    ens_u = calib_uniform.out["ensemble_final"]
    mu_gbar_u = [e[4] for e in ens_u]
    mu_p_u = [e[5] for e in ens_u]
    @test isapprox(std(mu_gbar_u), std(mu_p_u); rtol = 0.35)
end

println("OK: M10 prior 規則テスト green")
