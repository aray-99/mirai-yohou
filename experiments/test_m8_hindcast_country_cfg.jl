# COUNTRY_CFG の TOML 駆動ローダ(#0077, Issue #22)の回帰テスト。
#
# JPN/THA について、TOML(experiments/data/countries/*.toml)から構築した
# COUNTRY_CFG が旧ハードコードリテラルと完全一致することを確認する。
# 旧値はここにリテラルとして残し比較対象とする(Issue #22 の完了条件)。
# [hindcast] 欠落国のロードが明確なエラーになることも合わせて確認する。
#
# experiments/M8_hindcast.jl はスクリプト単体(パッケージテスト環境に無い
# Dates/DataFrames/JSON3/TOML 等に依存)のため、test/runtests.jl には
# 組み込まず experiments 側の簡易テストとして分離する(#0049 の
# test_m8_augmentation.jl と同じ流儀)。
#
# 実行: julia --project=experiments experiments/test_m8_hindcast_country_cfg.jl

using Test
using Dates

include(joinpath(@__DIR__, "M8_hindcast.jl"))

@testset "COUNTRY_CFG TOML 駆動ローダ(#0077, Issue #22)" begin
    @testset "JPN/THA が旧ハードコードリテラルと完全一致" begin
        jpn = COUNTRY_CFG["JPN"]
        @test jpn.regime == :stable
        @test jpn.calib == (5.0, 26.0)
        @test jpn.verif == (26.0, 35.0)
        @test jpn.exclude_admin1 == String[]
        @test jpn.acled_from == date_to_t(Date(2018, 1, 1))

        tha = COUNTRY_CFG["THA"]
        @test tha.regime == :volatile
        @test tha.calib == (20.0, 28.0)
        @test tha.verif == (28.0, 35.0)
        @test tha.exclude_admin1 == DEEP_SOUTH_THA
        @test tha.acled_from == date_to_t(Date(2010, 1, 1))

        @test Set(keys(COUNTRY_CFG)) == Set(["JPN", "THA"])
    end

    @testset "[hindcast] 欠落国は明確なエラー(Issue #20 参照)" begin
        cfg_no_hindcast = Dict(
            "regime" => "stable",
            "acled" => Dict(
                "exclude_admin1" => String[],
                "acled_from" => Date(2020, 1, 1),
            ),
        )
        err = try
            build_country_cfg("KOR", cfg_no_hindcast)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("KOR", err.msg)
        @test occursin("Issue #20", err.msg)
        @test occursin("[hindcast]", err.msg)
    end

    @testset "country_cfg: 実在する未確定国(KOR)も同じ明確なエラー" begin
        # countries/KOR.toml は実在するが [hindcast] 未確定(Issue #20 待ち)。
        # COUNTRY_CFG からは除外され(include は成功する)、参照時に
        # KeyError でなく build_country_cfg の明確なエラーになること。
        @test !haskey(COUNTRY_CFG, "KOR")
        err = try
            country_cfg("KOR")
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("Issue #20", err.msg)
    end
end

println("OK: COUNTRY_CFG TOML 駆動ローダ回帰テスト green")
