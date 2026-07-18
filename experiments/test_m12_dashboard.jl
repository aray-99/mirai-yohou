# M12-A ダッシュボード生成器骨格の単体テスト(Issue #11、docs/M12_DASHBOARD_DESIGN.md §5)。
#
# experiments/M12_dashboard.jl は experiments 環境の JSON3 に依存するため、
# test/runtests.jl(メインスイート)には組み込まず experiments 側の簡易テストと
# して分離する(test_m8_augmentation.jl / test_m10_prior.jl と同じ流儀)。
#
# 実行: julia --project=experiments experiments/test_m12_dashboard.jl

using Test

include(joinpath(@__DIR__, "M12_dashboard.jl"))

const FIXTURE_DIR = joinpath(@__DIR__, "test_fixtures")
const AAA_PATH = joinpath(FIXTURE_DIR, "M11_forecast_AAA.json")
const BBB_PATH = joinpath(FIXTURE_DIR, "M11_forecast_BBB.json")

@testset "M12 dashboard generator skeleton (#0074, Issue #11)" begin

    @testset "1. 生成成功" begin
        html = build_dashboard_html([AAA_PATH, BBB_PATH])
        @test html isa String
        @test length(html) > 0

        n_data_scripts = length(collect(eachmatch(r"<script type=\"application/json\" id=\"data-", html)))
        @test n_data_scripts == 2

        html_one = build_dashboard_html([AAA_PATH])
        n_one = length(collect(eachmatch(r"<script type=\"application/json\" id=\"data-", html_one)))
        @test n_one == 1
    end

    @testset "2. 決定性" begin
        html1 = build_dashboard_html([AAA_PATH, BBB_PATH])
        html2 = build_dashboard_html([AAA_PATH, BBB_PATH])
        @test html1 == html2
        @test codeunits(html1) == codeunits(html2)
    end

    @testset "3. 義務要件マーカー(データ駆動性の証明)" begin
        # (b) p_ex_fallback 警告バッジ: BBB(fallback発動)を含めば出現、
        # AAA のみ(fallback 非発動)なら出現しない。
        html_both = build_dashboard_html([AAA_PATH, BBB_PATH])
        @test occursin("ジャンプリスク評価不能", html_both)

        html_aaa_only = build_dashboard_html([AAA_PATH])
        @test !occursin("ジャンプリスク評価不能", html_aaa_only)

        # (a) 検証済み/外挿の境界ラベル等、チャート(SVG)側のマーカーは
        # Issue #12(実描画)のスコープ。ここではテストの置き場所のみ確保する。
        @test_skip occursin("外挿", html_both) && occursin("検証済み", html_both)
        # ↑ Issue #12 で SVG チャート境界ラベルを実装した際に有効化する。
    end

    @testset "4. 埋め込み忠実性" begin
        html = build_dashboard_html([AAA_PATH, BBB_PATH])
        raw_aaa = read(AAA_PATH, String)
        raw_bbb = read(BBB_PATH, String)
        @test occursin(raw_aaa, html)
        @test occursin(raw_bbb, html)
    end

end

println("OK: M12 dashboard generator skeleton テスト green")
