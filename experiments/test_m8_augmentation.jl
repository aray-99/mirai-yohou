# build_m8_augmented_params の国分岐(DECISIONS #0049)の単体テスト。
#
# experiments/M8_hindcast.jl はスクリプト単体(パッケージテスト環境に無い
# Dates/DataFrames/JSON3/TOML 等に依存)のため、test/runtests.jl には
# 組み込まず experiments 側の簡易テストとして分離する(package の 393 本
# には数えない)。
#
# 実行: julia --project=experiments experiments/test_m8_augmentation.jl

using Test
using MiraiYohou

include(joinpath(@__DIR__, "M8_hindcast.jl"))

@testset "build_m8_augmented_params (#0049)" begin
    params_jpn = build_params(COUNTRY_CFG["JPN"].regime)
    params_tha = build_params(COUNTRY_CFG["THA"].regime)

    aug_jpn = build_m8_augmented_params(params_jpn, "JPN")
    aug_tha = build_m8_augmented_params(params_tha, "THA")

    names_jpn = [p.name for p in aug_jpn]
    names_tha = [p.name for p in aug_tha]

    @test length(aug_jpn) == 5
    @test :theta_sig ∉ names_jpn
    @test :c_v0 ∈ names_jpn
    @test Set(names_jpn) == Set([:netgrowth, :a_T, :r_phi, :mu_gbar, :c_v0])

    @test length(aug_tha) == 6
    @test :theta_sig ∈ names_tha
    @test :c_v0 ∈ names_tha
    @test Set(names_tha) == Set([:theta_sig, :netgrowth, :a_T, :r_phi, :mu_gbar, :c_v0])

    # c_v0 の init は params.l2.c_v0(mu_gbar と同格に較正値から取得)
    cv0_jpn = only(filter(p -> p.name == :c_v0, aug_jpn))
    @test cv0_jpn.init == params_jpn.l2.c_v0
    @test cv0_jpn.link == :identity
    @test cv0_jpn.init_sd == 0.1
    @test cv0_jpn.rw_sd == 0.05
end

println("OK: build_m8_augmented_params 国分岐テスト green")
