# M12 ダッシュボード生成器の単体テスト(Issue #11/#14/#15、docs/M12_DASHBOARD_DESIGN.md §5)。
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

@testset "M12 dashboard generator (#0074, Issue #11/#14/#15)" begin

    @testset "1. 生成成功" begin
        html = build_dashboard_html([AAA_PATH, BBB_PATH])
        @test html isa String
        @test length(html) > 0

        n_data_scripts = length(collect(eachmatch(r"<script type=\"application/json\" id=\"data-", html)))
        @test n_data_scripts == 2

        html_one = build_dashboard_html([AAA_PATH])
        n_one = length(collect(eachmatch(r"<script type=\"application/json\" id=\"data-", html_one)))
        @test n_one == 1

        # JS 無効環境の可読性(§4.3): hidden な国パネルを noscript で表示に戻す
        @test occursin("<noscript><style>.country-panel[hidden]", html)
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

        # (a) 検証済み/外挿の視覚区別(Issue #12 実装分)。
        # クライアント側 JS テンプレートに境界実装のマーカー(クラス名・ラベル文字列)が
        # 含まれること、および data-until-calendar-year 属性がフィクスチャの
        # verified_horizon.until_calendar_year 値(2028)で出力されていることを検査する。
        @test occursin("外挿", html_both)
        @test occursin("検証済み", html_both)
        @test occursin("verified-boundary", html_both)
        @test occursin("検証済み | 外挿", html_both)
        @test occursin("外挿領域(不確実性の広がり)", html_both)
        @test occursin("data-until-calendar-year=\"2028\"", html_both)
    end

    @testset "4. 埋め込み忠実性" begin
        html = build_dashboard_html([AAA_PATH, BBB_PATH])
        raw_aaa = read(AAA_PATH, String)
        raw_bbb = read(BBB_PATH, String)
        @test occursin(raw_aaa, html)
        @test occursin(raw_bbb, html)
    end

    @testset "5. JA/EN トグル(#14)" begin
        html_both = build_dashboard_html([AAA_PATH, BBB_PATH])

        # トグル UI 自体(生成器が静的に出力)
        @test occursin("lang-toggle", html_both)
        @test occursin("data-lang-option=\"ja\"", html_both)
        @test occursin("data-lang-option=\"en\"", html_both)

        # 変数ラベル対訳(§3-3 正典表、一字一句厳守)。日英どちらも同じ生成 HTML に
        # 静的に埋め込まれ、表示切替はクライアント側 CSS が担う(§2.1 の責務分割に
        # 反しない — 国切替の hidden 属性と同じ static-both-render 方式)。
        @test occursin("総人口", html_both)
        @test occursin("Total population", html_both)
        @test occursin("生産年齢人口比率(15-64歳)", html_both)
        @test occursin("Working-age population ratio (ages 15-64)", html_both)
        @test occursin("政治的分極度(V-Dem v2cacamps、高いほど分裂)", html_both)
        @test occursin("Political polarization (V-Dem v2cacamps; higher = more polarized)", html_both)
        @test occursin("Social Deborah number (dimensionless)", html_both)

        # (b) 警告バッジのデータ駆動性は EN 文字列でも成立すること
        @test occursin("Jump risk unassessable", html_both)
        html_aaa_only = build_dashboard_html([AAA_PATH])
        @test !occursin("Jump risk unassessable", html_aaa_only)
        @test !occursin("ジャンプリスク評価不能", html_aaa_only)

        # (a) 検証済み/外挿の英語表記も存在すること
        @test occursin("Verified | Extrapolated", html_both)
        @test occursin("extrapolated region (spread of uncertainty)", html_both) ||
              occursin("the extrapolated region (spread of uncertainty)", html_both)
    end

    @testset "6. 軸文字可読性(#15、案 a: 動的 viewBox)" begin
        html = build_dashboard_html([AAA_PATH, BBB_PATH])

        # 選定した案(a)の実装マーカー: ResizeObserver でコンテナ実幅から
        # viewBox を再計算する(ブラウザ描画そのものは自動テストしない、§5 末尾)。
        @test occursin("ResizeObserver", html)
        @test occursin("measureContainerWidth", html)

        # 軸文字は実 px 換算で最低 10〜11px 相当(旧: viewBox 固定 640 に対し 9px
        # 決め打ちだったため、320px 幅パネルでは実質 4px 台まで縮んでいた)。
        @test occursin("AXIS_FONT_SIZE = \"11px\"", html)
        @test !occursin("font-size: 9px", html)
    end

end

println("OK: M12 dashboard generator テスト green")
