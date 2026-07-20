# M12-A: 予報ダッシュボード生成器骨格(Issue #11、docs/M12_DASHBOARD_DESIGN.md 規範)
#
# 「Julia スクリプト1本で決定的再生成」(#0070)の生成器。HTML/CSS/JS のテンプレートは
# すべてこのファイル内の文字列リテラルとして保持する(別ファイルに分離しない)。
#
# 生成コアは純関数 build_dashboard_html(json_paths::Vector{String}) -> String。
# CLI はその薄いラッパ(ISO3 引数 → experiments/output/M11_forecast_{ISO3}.json の
# パス解決、--out 省略時 docs/dashboard/index.html)。
#
# チャート SVG の実描画は別 Issue(#12)。本スクリプトが出力するのは
# データ駆動の静的 HTML(見出し・警告バッジ・表ビュー・来歴パネル・免責注記)+
# SVG 描画用の空コンテナ + 埋め込み JSON(生バイトのまま)のみ。
#
# 実行: julia --project=experiments experiments/M12_dashboard.jl JPN THA [--out docs/dashboard/index.html]

using JSON3

# ============================================================
# 定数
# ============================================================

# 設計書 §3-3 の正典対応表(英語変数名 — 日本語ラベル / English label)。
# 日本語ラベルは一字一句厳守。英語ラベルは SPEC §3 の定義に忠実な訳とする
# (#14: 独自の言い換え禁止)。
const VARIABLE_LABELS = [
    ("P", "総人口", "Total population"),
    ("w", "生産年齢人口比率(15-64歳)", "Working-age population ratio (ages 15-64)"),
    ("y", "一人当たり産出(実質GDP/capita、合成観測)", "Output per capita (real GDP per capita, composite observation)"),
    ("g_swiid", "格差指標(可処分ジニ)", "Inequality index (disposable-income Gini)"),
    ("T_proxy", "技術フロンティア(特許出願が観測代理)", "Technology frontier (patent applications as observation proxy)"),
    ("phi", "技術普及率(ネット利用率が観測代理)", "Technology diffusion rate (internet usage as observation proxy)"),
    ("v", "情報伝播速度(モバイル契約が観測代理)", "Information propagation speed (mobile subscriptions as observation proxy)"),
    ("tau", "制度信頼(WVS政府信頼)", "Institutional trust (WVS government confidence)"),
    ("p", "政治的分極度(V-Dem v2cacamps、高いほど分裂)", "Political polarization (V-Dem v2cacamps; higher = more polarized)"),
]

# De は §3-3 の 9 系列表には含まれない診断変数(設計書 3-5「De 診断パネル」)。
# SPEC §7 の正典用語(社会的デボラ数)に忠実な訳を UI 辞書として保持する。
const DE_LABEL_JA = "体質診断(無次元)"
const DE_LABEL_EN = "Social Deborah number (dimensionless)"

const REQUIRED_TOP_KEYS = [
    :country, :regime, :N, :horizon_years, :forecast_start, :verified_horizon,
    :finite_check, :assimilation, :jump_thinning, :variables, :count_forecast,
    :years, :t, :endogenous_jumps_total, :provenance, :elapsed_sec,
]

const QUANTILE_KEYS = ["q05", "q25", "q50", "q75", "q95"]

# ============================================================
# ユーティリティ
# ============================================================

"""HTML テキストノードとしてのエスケープ(属性値にも流用可能な最小集合)。"""
function html_escape(s::AbstractString)
    s = replace(s, "&" => "&amp;")
    s = replace(s, "<" => "&lt;")
    s = replace(s, ">" => "&gt;")
    s = replace(s, "\"" => "&quot;")
    return s
end

"""
UI 文字列を JA/EN 両方 static に埋め込むヘルパー(#14: JA/EN トグル)。
生成時点で両言語のテキストノードを埋め込み、表示切替はクライアント側 CSS/JS
(`.i18n-ja` / `.i18n-en` + `html[data-lang]`)が担う。JS 無効環境では
`.i18n-en` がデフォルトで非表示(CSS)のため、日本語主体の表示が保たれる
(設計書冒頭のユーザー確認事項「UI言語 = 日本語主体」に合致)。
"""
function i18n_span(ja::AbstractString, en::AbstractString)
    return """<span class="i18n-ja">$(html_escape(ja))</span><span class="i18n-en">$(html_escape(en))</span>"""
end

"""
`i18n_span` の非エスケープ版。信頼できる静的 HTML 断片(`<code>` 等を含む
免責注記の定数文字列)を埋め込む場合にのみ使う — 外部/JSON 由来の値には
絶対に使わないこと。
"""
function i18n_span_raw(ja_html::AbstractString, en_html::AbstractString)
    return """<span class="i18n-ja">$ja_html</span><span class="i18n-en">$en_html</span>"""
end

"""数値を表示用に整形(整数はそのまま、浮動小数は有効桁を抑える)。"""
function fmt_num(x)
    x === nothing && return "—"
    if x isa Integer
        return string(x)
    elseif x isa Real
        if isnan(x) || isinf(x)
            return string(x)
        end
        return string(round(x; sigdigits=6))
    end
    return string(x)
end

fmt_opt(x) = x === nothing ? "—" : string(x)

# ============================================================
# (i) JSON 検証
# ============================================================

"""トップレベルキーの存在確認。欠落があればエラーを投げる。"""
function validate_forecast_json(doc, path::AbstractString)
    missing_keys = [k for k in REQUIRED_TOP_KEYS if !haskey(doc, k)]
    if !isempty(missing_keys)
        error("JSON 検証失敗($path): トップレベルキー欠落 $(missing_keys)")
    end
    for (name, _, _) in VARIABLE_LABELS
        if !haskey(doc[:variables], Symbol(name))
            error("JSON 検証失敗($path): variables.$(name) が欠落")
        end
    end
    if !haskey(doc[:variables], :De)
        error("JSON 検証失敗($path): variables.De が欠落")
    end
    return nothing
end

# ============================================================
# (ii)/(iii) 埋め込み JSON・マニフェスト
# ============================================================

"""
入力 JSON ファイルの生バイト列(文字列として)を検証して返す。
`<` を含む場合は `</script>` インジェクションの恐れがあるためエラーにする
(設計書 §2.2: 入力は信頼できる自作 JSON である前提で `<` 非包含のみ検証)。
"""
function read_raw_json_bytes(path::AbstractString)
    raw = read(path, String)
    if occursin("<", raw)
        error("埋め込み対象 JSON に '<' が含まれるため埋め込み不可: $path")
    end
    return raw
end

# ============================================================
# (iv) データ駆動レンダリング
# ============================================================

function render_header(countries::Vector{String}, docs::Dict{String,Any})
    first_country = countries[1]
    doc0 = docs[first_country]

    tabs = join([
        """<button type="button" class="country-tab" data-country="$(html_escape(c))">$(html_escape(c))</button>"""
        for c in countries
    ], "\n      ")

    title_i18n = i18n_span("未来予報 — 社会動態シミュレータ 将来予報", "Future Forecast — Social Dynamics Simulator")
    nav_label_ja = html_escape("国切替")
    nav_label_en = html_escape("Country switch")

    io = IOBuffer()
    print(io, """
    <header class="dashboard-header">
      <div class="header-top">
        <h1>$title_i18n</h1>
        <div class="lang-toggle" role="group" data-aria-ja="UI言語" data-aria-en="UI language" aria-label="UI言語">
          <button type="button" class="lang-btn" data-lang-option="ja" aria-pressed="true">JA</button>
          <button type="button" class="lang-btn" data-lang-option="en" aria-pressed="false">EN</button>
        </div>
      </div>
      <nav class="country-tabs" id="country-tabs" role="tablist" data-aria-ja="$nav_label_ja" data-aria-en="$nav_label_en" aria-label="$nav_label_ja">
      $tabs
      </nav>
      <div class="summary-line" id="summary-line" data-template="true">
        <span class="regime-badge" id="regime-badge"></span>
        <span id="forecast-summary"></span>
      </div>
    </header>
    """)
    return String(take!(io))
end

function panel_footnote_variables()
    text = i18n_span(
        "y=0 は予報起点断面(予報積分前のアンサンブル事後分布)であり、y=1 以降が実際の予報値である。",
        "y=0 is the forecast-origin cross-section (the ensemble posterior before forecast integration); y=1 onward are the actual forecast values.",
    )
    return """<p class="panel-footnote">$text</p>"""
end

function render_variable_table(country::String, varname::String, entry, years)
    n = length(years)
    rows = IOBuffer()
    for i in 1:n
        vals = [fmt_num(entry[Symbol(q)][i]) for q in QUANTILE_KEYS]
        print(rows, "<tr><td>$(years[i])</td>" * join(["<td>$(v)</td>" for v in vals]) * "</tr>\n")
    end
    unit = html_escape(string(entry[:unit]))
    last_obs = entry[:last_observation_year]
    last_obs_str = last_obs === nothing ? i18n_span("—(状態変数でない)", "— (not a state variable)") : string(last_obs)
    table_view_label = i18n_span("表ビュー", "Table view")
    unit_label = i18n_span("単位", "Unit")
    last_obs_label = i18n_span("最終観測年", "Last observation year")
    year_header = i18n_span("年", "Year")
    return """
        <details class="table-view">
          <summary>$table_view_label($(html_escape(varname)))</summary>
          <p class="table-meta">$unit_label: $unit / $last_obs_label: $last_obs_str</p>
          <table>
            <thead><tr><th>$year_header</th><th>q05</th><th>q25</th><th>q50</th><th>q75</th><th>q95</th></tr></thead>
            <tbody>
            $(String(take!(rows)))
            </tbody>
          </table>
        </details>
    """
end

function render_variable_panel(country::String, varname::String, label_ja::String, label_en::String, entry, years)
    unit = string(entry[:unit])
    transform = string(entry[:transform])
    table = render_variable_table(country, varname, entry, years)
    label_i18n = i18n_span(label_ja, label_en)
    transform_label = i18n_span("変換", "Transform")
    return """
      <section class="panel variable-panel" data-country="$(html_escape(country))" data-variable="$(html_escape(varname))">
        <h3>$(html_escape(varname)) — $label_i18n($(html_escape(unit)))</h3>
        <p class="panel-transform">$transform_label: $(html_escape(transform))</p>
        <div class="svg-container" data-role="chart" data-variable="$(html_escape(varname))"></div>
        $table
      </section>
    """
end

function render_count_forecast_panel(country::String, doc)
    cf = doc[:count_forecast]
    jt = doc[:jump_thinning]
    yrs_full = doc[:years]
    # count_forecast は y=1..horizon_years(配列長 horizon_years、variables より1短い)
    cf_years = yrs_full[2:end]

    badge = ""
    if get(jt, :p_ex_fallback, nothing) == "no_count_data"
        badge_text = i18n_span(
            "⚠ ジャンプリスク評価不能(同化窓内にカウントデータなし)",
            "⚠ Jump risk unassessable (no count data within assimilation window)",
        )
        badge = """<span class="warning-badge" role="status">$badge_text</span>"""
    end

    n = length(cf_years)
    rows = IOBuffer()
    for i in 1:n
        vals = [fmt_num(cf[Symbol(q)][i]) for q in QUANTILE_KEYS]
        print(rows, "<tr><td>$(cf_years[i])</td>" * join(["<td>$(v)</td>" for v in vals]) * "</tr>\n")
    end

    nu_star = fmt_num(get(cf, :nu_star, nothing))
    r_hat = get(cf, :r_hat, nothing)
    r_hat_str = r_hat === nothing ? "Poisson" : fmt_num(r_hat)

    title_i18n = i18n_span("期待騒乱イベント数 / 年", "Expected unrest event count / year")
    unit_label = i18n_span("単位", "Unit")
    table_view_label = i18n_span("表ビュー", "Table view")
    year_header = i18n_span("年", "Year")

    return """
      <section class="panel count-forecast-panel" data-country="$(html_escape(country))">
        <h3>$title_i18n $badge</h3>
        <p class="panel-note">$unit_label: $(html_escape(string(get(cf, :unit, ""))))</p>
        <p class="panel-footnote">nu_star = $nu_star / r_hat = $r_hat_str</p>
        <div class="svg-container" data-role="chart" data-variable="count_forecast"></div>
        <details class="table-view">
          <summary>$table_view_label(count_forecast)</summary>
          <table>
            <thead><tr><th>$year_header</th><th>q05</th><th>q25</th><th>q50</th><th>q75</th><th>q95</th></tr></thead>
            <tbody>
            $(String(take!(rows)))
            </tbody>
          </table>
        </details>
      </section>
    """
end

function render_de_panel(country::String, doc)
    entry = doc[:variables][:De]
    years = doc[:years]
    table = render_variable_table(country, "De", entry, years)
    title_i18n = i18n_span(DE_LABEL_JA, DE_LABEL_EN)
    footnote_i18n = i18n_span(
        "De=1: 安定/変動レジーム境界。外挿領域では個別値でなくレジーム傾向(De の帯が 1 を跨ぐか)を読む。",
        "De=1: boundary between stable and volatile regimes. In the extrapolated region, read the regime tendency (whether the De band crosses 1) rather than individual values.",
    )
    return """
      <section class="panel de-panel" data-country="$(html_escape(country))">
        <h3>De — $title_i18n</h3>
        <p class="panel-footnote">$footnote_i18n</p>
        <div class="svg-container" data-role="chart" data-variable="De"></div>
        $table
      </section>
    """
end

function render_provenance_panel(country::String, doc)
    prov = doc[:provenance]
    asm = doc[:assimilation]
    commit_full = string(get(prov, :commit, ""))
    commit_short = length(commit_full) >= 8 ? commit_full[1:8] : commit_full
    ess_range = get(asm, :ess_range, nothing)
    ess_str = ess_range === nothing ? "—" : "[$(fmt_num(ess_range[1])), $(fmt_num(ess_range[2]))]"
    window = get(asm, :window, nothing)
    window_str = window === nothing ? "—" : "[$(fmt_num(window[1])), $(fmt_num(window[2]))]"

    fetched = get(prov, :data_fetched_at, nothing)
    fetched_str = if fetched === nothing
        "—"
    else
        join(["$(k): $(v)" for (k, v) in pairs(fetched)], " / ")
    end

    provenance_title = i18n_span("来歴", "Provenance")
    return """
      <section class="panel provenance-panel" data-country="$(html_escape(country))">
        <h3>$provenance_title</h3>
        <dl class="provenance-list">
          <dt>commit</dt><dd title="$(html_escape(commit_full))">$(html_escape(commit_short))</dd>
          <dt>seed</dt><dd>$(fmt_opt(get(prov, :seed, nothing)))</dd>
          <dt>generated_at</dt><dd>$(html_escape(string(get(prov, :generated_at, ""))))</dd>
          <dt>data_fetched_at</dt><dd>$(html_escape(fetched_str))</dd>
          <dt>frozen_config</dt><dd>$(html_escape(string(get(prov, :frozen_config, ""))))</dd>
          <dt>script</dt><dd>$(html_escape(string(get(prov, :script, ""))))</dd>
          <dt>design_decision</dt><dd>$(html_escape(string(get(prov, :design_decision, ""))))</dd>
          <dt>frozen_decisions</dt><dd>$(html_escape(string(get(prov, :frozen_decisions, ""))))</dd>
          <dt>finite_check.pass</dt><dd>$(fmt_opt(get(doc[:finite_check], :pass, nothing)))</dd>
          <dt>endogenous_jumps_total</dt><dd>$(fmt_opt(get(doc, :endogenous_jumps_total, nothing)))</dd>
          <dt>assimilation.n_obs</dt><dd>$(fmt_opt(get(asm, :n_obs, nothing)))</dd>
          <dt>assimilation.window</dt><dd>$window_str</dd>
          <dt>assimilation.ess_range</dt><dd>$ess_str</dd>
          <dt>assimilation.nu_star</dt><dd>$(fmt_num(get(asm, :nu_star, nothing)))</dd>
        </dl>
      </section>
    """
end

const DISCLAIMER_TITLE_I18N = i18n_span("免責注記", "Disclaimer")

const DISCLAIMER_PARA_1_I18N = i18n_span_raw(
    """外挿領域(検証済みホライズンより先、y &gt; 6)の分位帯は、モデルの内生力学
          (Hawkes 自己励起・OU 型緩和・L3 拡大パラメータの継続ランダムウォーク)のみに
          基づく統計的投影であり、観測による制約を一切受けていない。""",
    """In the extrapolated region (beyond the verified horizon, y &gt; 6), the quantile
          bands are a purely statistical projection driven by the model's endogenous
          dynamics (Hawkes self-excitation, OU-type relaxation, ongoing random walk of
          the L3 expansion parameter), unconstrained by any observation.""",
)

const DISCLAIMER_PARA_2_I18N = i18n_span_raw(
    """特に <code>p_ex_fallback: "no_count_data"</code> が発動している国・オリジンでは、
          ジャンプ関連の力学(Γ 適用による急激な社会的応力上昇)が一切表現されないため、
          <code>count_forecast</code> は「ジャンプが起きない前提」の過小評価である
          可能性がある。""",
    """In particular, for country/origin combinations where
          <code>p_ex_fallback: "no_count_data"</code> is triggered, jump-related dynamics
          (a sudden rise in social stress via Γ) are not represented at all, so
          <code>count_forecast</code> may be an underestimate that assumes "no jumps
          occur".""",
)

const DISCLAIMER_PARA_3_I18N = i18n_span_raw(
    """検証の実証裏付けは M9/M10 の expanding-window walk-forward
          (#0052 設計・#0069 合格判定)であり、検証済み範囲は起点+6年に限られる。""",
    """The empirical basis for verification is the M9/M10 expanding-window walk-forward
          (design #0052, pass verdict #0069), and the verified range is limited to
          origin + 6 years.""",
)

const DISCLAIMER_TEXT = """
      <section class="panel disclaimer-panel">
        <h3>$DISCLAIMER_TITLE_I18N</h3>
        <p>
          $DISCLAIMER_PARA_1_I18N
        </p>
        <p>
          $DISCLAIMER_PARA_2_I18N
        </p>
        <p>
          $DISCLAIMER_PARA_3_I18N
        </p>
      </section>
"""

function render_country_section(country::String, doc)
    io = IOBuffer()
    print(io, """<div class="country-panel" id="country-$(html_escape(country))" data-country="$(html_escape(country))" hidden>\n""")

    forecast_start = doc[:forecast_start]
    verified = doc[:verified_horizon]
    print(io, """
      <div class="country-meta" data-country="$(html_escape(country))"
           data-regime="$(html_escape(string(doc[:regime])))"
           data-forecast-start-year="$(fmt_opt(get(forecast_start, :calendar_year, nothing)))"
           data-horizon-years="$(fmt_opt(get(doc, :horizon_years, nothing)))"
           data-n="$(fmt_opt(get(doc, :N, nothing)))"
           data-until-calendar-year="$(fmt_opt(get(verified, :until_calendar_year, nothing)))"
      ></div>
    """)

    print(io, """<div class="panel-grid">\n""")
    for (varname, label_ja, label_en) in VARIABLE_LABELS
        entry = doc[:variables][Symbol(varname)]
        print(io, render_variable_panel(country, varname, label_ja, label_en, entry, doc[:years]))
    end
    print(io, panel_footnote_variables())
    print(io, "</div>\n")

    print(io, render_count_forecast_panel(country, doc))
    print(io, render_de_panel(country, doc))
    print(io, render_provenance_panel(country, doc))

    print(io, "</div>\n")
    return String(take!(io))
end

# ============================================================
# CSS / JS(インライン、外部依存ゼロ)
# ============================================================

const DASHBOARD_CSS = """
/* ============================================================
   デザイントークン写経(#0081: DADS design-tokens を視覚正典として採用)
   出典: https://github.com/digital-go-jp/design-tokens (MIT License)
   commit: 298ea8c349721db463dd3feb804e5b8197061e5c
   file:   examples/tokens.css
   選抜方針: 色・タイポ・エレベーションのトークンのみ写経する。font-family は
   DADS 既定の 'Noto Sans JP'(Web フォント)を採用せず、自己完結制約(#0070)
   により既存のシステムフォントスタックを維持する。
   ============================================================ */
:root {
  /* --- DADS プリミティブ・トークン(選抜) --- */
  --dads-blue-400: #7096f8;
  --dads-blue-700: #264af4;
  --dads-gray-50: #f2f2f2;
  --dads-gray-100: #e6e6e6;
  --dads-gray-200: #cccccc;
  --dads-gray-300: #b3b3b3;
  --dads-gray-400: #999999;
  --dads-gray-536: #767676;
  --dads-gray-700: #4d4d4d;
  --dads-gray-800: #333333;
  --dads-gray-900: #1a1a1a;
  --dads-white: #ffffff;
  --dads-black: #000000;
  --dads-yellow-400: #ffc700;
  --dads-font-size-14: 0.875rem;
  --dads-font-size-16: 1rem;
  --dads-font-size-20: 1.25rem;
  --dads-line-height-150: 1.5;
  --dads-radius-4: 0.25rem;
  --dads-radius-8: 0.5rem;
  --dads-radius-full: 624.9375rem;
  --dads-elevation-1: 0 2px 8px 1px rgba(0,0,0,0.1), 0 1px 5px 0 rgba(0,0,0,0.3);

  /* --- セマンティック・トークン(ダッシュボード用途へのマッピング) ---
     系列色は DADS Blue、警告色は DADS Yellow へ写像(#0081 決定1)。
     単一色相原則(設計書 §4.2: 全パネルが同一アンサンブルの分位のため
     --color-series 一色のみを使う)は変更しない。中立色は DADS neutral
     gray スケールへ統一する(作業内容 3)。 */
  --color-page-bg: var(--dads-gray-50);
  --color-panel-bg: var(--dads-white);
  --color-text-primary: var(--dads-black);
  --color-text-secondary: var(--dads-gray-700);
  --color-axis-label: var(--dads-gray-536);
  --color-gridline: var(--dads-gray-100);
  --color-axis-line: var(--dads-gray-300);
  --color-series: var(--dads-blue-700);
  --color-warning: var(--dads-yellow-400);
}
@media (prefers-color-scheme: dark) {
  :root {
    --color-page-bg: var(--dads-black);
    --color-panel-bg: var(--dads-gray-900);
    --color-text-primary: var(--dads-white);
    --color-text-secondary: var(--dads-gray-200);
    --color-axis-label: var(--dads-gray-400);
    --color-gridline: var(--dads-gray-800);
    --color-axis-line: var(--dads-gray-700);
    --color-series: var(--dads-blue-400);
    --color-warning: var(--dads-yellow-400);
  }
}
* { box-sizing: border-box; }
body {
  margin: 0;
  font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
  font-size: var(--dads-font-size-16);
  line-height: var(--dads-line-height-150);
  background: var(--color-page-bg);
  color: var(--color-text-primary);
}
.dashboard-header {
  padding: 1rem 1.25rem;
  border-bottom: 1px solid var(--color-gridline);
}
.dashboard-header h1 { font-size: var(--dads-font-size-20); line-height: 1.3; margin: 0; }
.header-top {
  display: flex;
  align-items: baseline;
  justify-content: space-between;
  gap: 0.75rem;
  flex-wrap: wrap;
  margin-bottom: 0.5rem;
}
/* #14: JA/EN トグル。デフォルト(JS 無効・data-lang 未設定)は日本語主体表示。 */
.i18n-en { display: none; }
html[data-lang="en"] .i18n-ja { display: none; }
html[data-lang="en"] .i18n-en { display: inline; }
.lang-toggle { display: inline-flex; gap: 0.25rem; flex: none; }
.lang-btn {
  font: inherit;
  font-size: var(--dads-font-size-14);
  padding: 0.2rem 0.6rem;
  border: 1px solid var(--color-axis-line);
  border-radius: var(--dads-radius-4);
  background: var(--color-panel-bg);
  color: var(--color-text-primary);
  cursor: pointer;
}
.lang-btn[aria-pressed="true"] { border-color: var(--color-series); font-weight: 700; }
.country-tabs { display: flex; gap: 0.5rem; margin-bottom: 0.5rem; }
.country-tab {
  font: inherit;
  padding: 0.35rem 0.9rem;
  border: 1px solid var(--color-axis-line);
  border-radius: var(--dads-radius-full);
  background: var(--color-panel-bg);
  color: var(--color-text-primary);
  cursor: pointer;
}
.country-tab[aria-selected="true"] { border-color: var(--color-series); font-weight: 700; }
.summary-line { color: var(--color-text-secondary); font-size: var(--dads-font-size-14); }
.regime-badge {
  display: inline-block;
  padding: 0.1rem 0.6rem;
  border: 1px solid var(--color-axis-line);
  border-radius: var(--dads-radius-4);
  margin-right: 0.5rem;
  color: var(--color-text-primary);
}
main { padding: 1rem 1.25rem; }
.panel-grid {
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(320px, 1fr));
  gap: 1rem;
}
.panel {
  background: var(--color-panel-bg);
  border: 1px solid var(--color-gridline);
  border-radius: var(--dads-radius-8);
  box-shadow: var(--dads-elevation-1);
  padding: 0.75rem 1rem;
  margin-bottom: 1rem;
}
.panel h3 { margin: 0 0 0.35rem 0; font-size: var(--dads-font-size-16); }
.panel-transform, .panel-note, .panel-footnote, .table-meta {
  color: var(--color-text-secondary);
  font-size: var(--dads-font-size-14);
  margin: 0.15rem 0;
}
.svg-container { min-height: 2rem; position: relative; }
.chart-svg { width: 100%; height: auto; display: block; }
.chart-overlay { cursor: crosshair; }
.chart-tooltip {
  position: absolute;
  pointer-events: none;
  background: var(--color-panel-bg);
  border: 1px solid var(--color-axis-line);
  border-radius: var(--dads-radius-4);
  padding: 0.3rem 0.5rem;
  font-size: 0.72rem;
  line-height: 1.3;
  color: var(--color-text-primary);
  box-shadow: var(--dads-elevation-1);
  z-index: 5;
  white-space: nowrap;
}
.warning-badge {
  display: inline-block;
  background: var(--color-warning);
  color: var(--dads-black);
  border-radius: var(--dads-radius-4);
  padding: 0.1rem 0.5rem;
  font-size: 0.85rem;
  margin-left: 0.5rem;
}
table { border-collapse: collapse; width: 100%; font-size: var(--dads-font-size-14); margin-top: 0.4rem; }
th, td { border: 1px solid var(--color-gridline); padding: 0.2rem 0.4rem; text-align: right; }
th:first-child, td:first-child { text-align: left; }
details.table-view summary { cursor: pointer; color: var(--color-text-secondary); font-size: var(--dads-font-size-14); }
.provenance-list { display: grid; grid-template-columns: max-content 1fr; gap: 0.15rem 0.75rem; font-size: 0.82rem; }
.provenance-list dt { color: var(--color-text-secondary); }
.provenance-list dd { margin: 0; word-break: break-word; }
.disclaimer-panel p { font-size: var(--dads-font-size-14); line-height: var(--dads-line-height-150); }
footer.legend-bar { padding: 0.5rem 1.25rem; color: var(--color-text-secondary); font-size: var(--dads-font-size-14); }
.legend-item { display: flex; align-items: center; gap: 0.4rem; margin: 0.2rem 0; }
.legend-swatch { flex: none; vertical-align: middle; }
"""

const DASHBOARD_JS = """
(function () {
  "use strict";

  var VARIABLE_ORDER = ["P", "w", "y", "g_swiid", "T_proxy", "phi", "v", "tau", "p", "De"];

  // #15(案 a): viewBox をコンテナ実幅から動的決定する。CHART_WIDTH_DEFAULT/
  // CHART_HEIGHT_DEFAULT はコンテナ幅が測定不能な場合(レイアウト前・非表示)の
  // フォールバックのみに使う。アスペクト比はこの既定値を保つ。
  var CHART_WIDTH_DEFAULT = 640;
  var CHART_HEIGHT_DEFAULT = 260;
  var CHART_ASPECT = CHART_HEIGHT_DEFAULT / CHART_WIDTH_DEFAULT;
  var MIN_CHART_WIDTH = 160;
  var AXIS_FONT_SIZE = "11px";
  var MARGIN = { top: 20, right: 16, bottom: 32, left: 60 };

  // #14: チャート JS が直接生成するテキスト(生成器がデータ駆動で埋め込めない
  // ラベル)の対訳辞書。生成 HTML 自体はどちらの言語でも同一(決定性を維持)で、
  // 表示切替はここから読むランタイム値のみに依存する。
  var STRINGS = {
    ja: {
      ariaChart: "予報ファンチャート",
      verified: "検証済み",
      extrapolated: "外挿",
      boundaryLabel: "検証済み | 外挿",
      deThresholdLabel: "De=1: 安定/変動レジーム境界",
      regimePrefix: "regime: ",
      tooltipYear: function (year, status) { return String(year) + "年(" + status + ")"; },
      summaryLine: function (year, horizon, n) {
        return "予報起点 " + year + " / ホライズン " + horizon + " 年 / N=" + n;
      }
    },
    en: {
      ariaChart: "Forecast fan chart",
      verified: "Verified",
      extrapolated: "Extrapolated",
      boundaryLabel: "Verified | Extrapolated",
      deThresholdLabel: "De=1: stable/volatile regime boundary",
      regimePrefix: "regime: ",
      tooltipYear: function (year, status) { return String(year) + " (" + status + ")"; },
      summaryLine: function (year, horizon, n) {
        return "Forecast start " + year + " / Horizon " + horizon + " years / N=" + n;
      }
    }
  };
  var LANG_STORAGE_KEY = "mirai-yohou-dashboard-lang";

  function currentLang() {
    var l = document.documentElement.getAttribute("data-lang");
    return l === "en" ? "en" : "ja";
  }

  function t() {
    return STRINGS[currentLang()];
  }

  var renderedCountries = {};

  function readManifest() {
    var el = document.getElementById("manifest");
    if (!el) return { countries: [] };
    return JSON.parse(el.textContent);
  }

  function readCountryData(country) {
    var el = document.getElementById("data-" + country);
    if (!el) return null;
    return JSON.parse(el.textContent);
  }

  function svgEl(name, attrs) {
    var el = document.createElementNS("http://www.w3.org/2000/svg", name);
    if (attrs) {
      for (var k in attrs) {
        if (Object.prototype.hasOwnProperty.call(attrs, k)) {
          el.setAttribute(k, attrs[k]);
        }
      }
    }
    return el;
  }

  function isFractionLikeUnit(unit) {
    if (!unit) return false;
    return unit.indexOf("fraction") !== -1 || unit.indexOf("dimensionless") !== -1;
  }

  function formatAxisValue(v, unit) {
    if (isFractionLikeUnit(unit)) return v.toFixed(2);
    var abs = Math.abs(v);
    if (abs >= 1e9) return (v / 1e9).toFixed(1) + "G";
    if (abs >= 1e6) return (v / 1e6).toFixed(1) + "M";
    if (abs >= 1e3) return (v / 1e3).toFixed(1) + "k";
    if (abs >= 10) return v.toFixed(0);
    return v.toFixed(2);
  }

  function formatTooltipValue(v, unit) {
    if (isFractionLikeUnit(unit)) return v.toFixed(3);
    return formatAxisValue(v, unit);
  }

  function linspace(a, b, n) {
    var out = [];
    if (n === 1) { out.push(a); return out; }
    var step = (b - a) / (n - 1);
    for (var i = 0; i < n; i++) out.push(a + step * i);
    return out;
  }

  function computeYDomain(quantiles, extraValues) {
    var min = Infinity, max = -Infinity;
    ["q05", "q25", "q50", "q75", "q95"].forEach(function (key) {
      (quantiles[key] || []).forEach(function (v) {
        if (v < min) min = v;
        if (v > max) max = v;
      });
    });
    (extraValues || []).forEach(function (v) {
      if (v < min) min = v;
      if (v > max) max = v;
    });
    if (!isFinite(min) || !isFinite(max)) { min = 0; max = 1; }
    if (min === max) { min -= 1; max += 1; }
    var pad = (max - min) * 0.05;
    return [min - pad, max + pad];
  }

  function makeScaleX(years, width) {
    var y0 = years[0], y1 = years[years.length - 1];
    var x0 = MARGIN.left, x1 = width - MARGIN.right;
    var span = (y1 - y0) || 1;
    return function (year) {
      return x0 + ((year - y0) / span) * (x1 - x0);
    };
  }

  function makeScaleY(domain, height) {
    var y0 = MARGIN.top, y1 = height - MARGIN.bottom;
    var span = (domain[1] - domain[0]) || 1;
    return function (v) {
      return y1 - ((v - domain[0]) / span) * (y1 - y0);
    };
  }

  function measureContainerWidth(container) {
    var w = container.getBoundingClientRect().width;
    if (!w || w <= 0) return CHART_WIDTH_DEFAULT;
    return Math.max(MIN_CHART_WIDTH, Math.round(w));
  }

  // #15(案 a): ResizeObserver でコンテナ幅の変化を検知し、しきい値(2px)を
  // 超えたときのみ再描画する(サブピクセルの発火・再描画ループを避ける)。
  function attachResponsiveRedraw(container, redraw) {
    if (container._responsiveAttached) return;
    container._responsiveAttached = true;
    if (typeof ResizeObserver === "undefined") return;
    var lastWidth = null;
    var pending = false;
    var ro = new ResizeObserver(function () {
      if (pending) return;
      pending = true;
      window.requestAnimationFrame(function () {
        pending = false;
        var w = container.getBoundingClientRect().width;
        if (w > 0 && (lastWidth === null || Math.abs(w - lastWidth) > 2)) {
          lastWidth = w;
          redraw();
        }
      });
    });
    ro.observe(container);
  }

  function buildPolygonPoints(years, upper, lower, scaleX, scaleY) {
    var pts = [];
    var i;
    for (i = 0; i < years.length; i++) {
      pts.push(scaleX(years[i]) + "," + scaleY(upper[i]));
    }
    for (i = years.length - 1; i >= 0; i--) {
      pts.push(scaleX(years[i]) + "," + scaleY(lower[i]));
    }
    return pts.join(" ");
  }

  function buildLinePoints(years, values, scaleX, scaleY) {
    var pts = [];
    for (var i = 0; i < years.length; i++) {
      pts.push(scaleX(years[i]) + "," + scaleY(values[i]));
    }
    return pts.join(" ");
  }

  function xTicks(years) {
    var y0 = years[0], y1 = years[years.length - 1];
    var start = Math.ceil(y0 / 5) * 5;
    var ticks = [];
    for (var y = start; y <= y1; y += 5) ticks.push(y);
    if (ticks.length === 0) { ticks.push(y0); ticks.push(y1); }
    return ticks;
  }

  function clearContainer(container) {
    while (container.firstChild) container.removeChild(container.firstChild);
  }

  function renderFanChart(container, opts) {
    clearContainer(container);
    var years = opts.years;
    var quantiles = opts.quantiles;
    var unit = opts.unit || "";
    var untilYear = opts.untilYear;
    var threshold = opts.threshold;

    // #15(案 a): viewBox をコンテナの実測幅から都度決定する。1 viewBox 単位が
    // 実際の 1px に一致するため、font-size を実 px として指定でき、320px 幅の
    // グリッド縮小パネルでも軸文字が AXIS_FONT_SIZE のまま読める。
    var width = measureContainerWidth(container);
    var height = Math.round(width * CHART_ASPECT);

    var extraForDomain = threshold ? [threshold.value] : [];
    var domain = computeYDomain(quantiles, extraForDomain);
    var scaleX = makeScaleX(years, width);
    var scaleY = makeScaleY(domain, height);

    var svg = svgEl("svg", {
      viewBox: "0 0 " + width + " " + height,
      width: "100%",
      height: "auto",
      preserveAspectRatio: "xMidYMid meet",
      role: "img",
      "aria-label": t().ariaChart,
      "class": "chart-svg"
    });

    var plotX0 = MARGIN.left, plotX1 = width - MARGIN.right;
    var plotY0 = MARGIN.top, plotY1 = height - MARGIN.bottom;

    var hasSplit = untilYear !== null && untilYear !== undefined &&
      untilYear > years[0] && untilYear < years[years.length - 1];
    var boundaryX = hasSplit ? scaleX(untilYear) : null;

    if (hasSplit) {
      var wash = svgEl("rect", {
        x: boundaryX, y: plotY0, width: (plotX1 - boundaryX), height: (plotY1 - plotY0),
        style: "fill: var(--color-gridline); fill-opacity: 0.5;",
        "class": "extrapolation-wash"
      });
      svg.appendChild(wash);
    }

    // しきい値線(De=1 等)のピクセル位置・ラベル位置を先に求めておく。近接する
    // 通常の目盛はラベルを間引いて重なりを避ける(線自体が既にその値を示すため
    // 目盛ラベルは冗長。ラベルは上端(境界ラベルの行)に来る場合だけ線の下に
    // 逃がす — #15: 実 px フォント化で文字が大きくなった分、重なりが顕在化
    // しやすいため)。
    var thresholdY = threshold ? scaleY(threshold.value) : null;
    var thresholdLabelY = threshold
      ? ((thresholdY < (plotY0 + 16)) ? (thresholdY + 14) : (thresholdY - 4))
      : null;
    var tickRowHeight = (plotY1 - plotY0) / 4;

    var yTickValues = linspace(domain[0], domain[1], 5);
    yTickValues.forEach(function (v) {
      var yPix = scaleY(v);
      var line = svgEl("line", {
        x1: plotX0, x2: plotX1, y1: yPix, y2: yPix,
        style: "stroke: var(--color-gridline); stroke-width: 1;",
        "class": "chart-gridline"
      });
      svg.appendChild(line);
      var collidesThreshold = thresholdY !== null && (
        Math.abs(yPix - thresholdY) < tickRowHeight * 0.6 ||
        Math.abs(yPix - thresholdLabelY) < tickRowHeight * 0.6
      );
      if (collidesThreshold) {
        return;
      }
      var label = svgEl("text", {
        x: plotX0 - 6, y: yPix, "text-anchor": "end", "dominant-baseline": "middle",
        style: "fill: var(--color-axis-label); font-size: " + AXIS_FONT_SIZE + ";",
        "class": "chart-axis-label"
      });
      label.textContent = formatAxisValue(v, unit);
      svg.appendChild(label);
    });

    var xt = xTicks(years);
    xt.forEach(function (yr) {
      var xPix = scaleX(yr);
      var tick = svgEl("line", {
        x1: xPix, x2: xPix, y1: plotY1, y2: (plotY1 + 4),
        style: "stroke: var(--color-axis-line); stroke-width: 1;"
      });
      svg.appendChild(tick);
      var label = svgEl("text", {
        x: xPix, y: (plotY1 + 14), "text-anchor": "middle",
        style: "fill: var(--color-axis-label); font-size: " + AXIS_FONT_SIZE + ";",
        "class": "chart-axis-label"
      });
      label.textContent = String(yr);
      svg.appendChild(label);
    });

    var uid = "clip-" + Math.random().toString(36).slice(2, 10);
    var defs = svgEl("defs");
    var clipVerified = svgEl("clipPath", { id: (uid + "-verified") });
    clipVerified.appendChild(svgEl("rect", {
      x: plotX0, y: plotY0,
      width: (hasSplit ? (boundaryX - plotX0) : (plotX1 - plotX0)),
      height: (plotY1 - plotY0)
    }));
    var clipExtrapolated = svgEl("clipPath", { id: (uid + "-extrapolated") });
    clipExtrapolated.appendChild(svgEl("rect", {
      x: (hasSplit ? boundaryX : plotX1), y: plotY0,
      width: (hasSplit ? (plotX1 - boundaryX) : 0),
      height: (plotY1 - plotY0)
    }));
    defs.appendChild(clipVerified);
    defs.appendChild(clipExtrapolated);
    svg.appendChild(defs);

    var band90Points = buildPolygonPoints(years, quantiles.q95, quantiles.q05, scaleX, scaleY);
    var band50Points = buildPolygonPoints(years, quantiles.q75, quantiles.q25, scaleX, scaleY);
    var linePoints = buildLinePoints(years, quantiles.q50, scaleX, scaleY);

    function buildDataGroup(clipId, opacityFactor) {
      var g = svgEl("g", { "clip-path": ("url(#" + clipId + ")") });
      g.appendChild(svgEl("polygon", {
        points: band90Points,
        style: ("fill: var(--color-series); fill-opacity: " + (0.16 * opacityFactor) + "; stroke: none;"),
        "class": "chart-band-90"
      }));
      g.appendChild(svgEl("polygon", {
        points: band50Points,
        style: ("fill: var(--color-series); fill-opacity: " + (0.32 * opacityFactor) + "; stroke: none;"),
        "class": "chart-band-50"
      }));
      g.appendChild(svgEl("polyline", {
        points: linePoints,
        style: ("fill: none; stroke: var(--color-series); stroke-width: 2; stroke-linecap: round; stroke-linejoin: round; stroke-opacity: " + opacityFactor + ";"),
        "class": "chart-median-line"
      }));
      return g;
    }

    svg.appendChild(buildDataGroup((uid + "-verified"), 1));
    if (hasSplit) {
      svg.appendChild(buildDataGroup((uid + "-extrapolated"), 0.55));
    }

    if (threshold) {
      var tline = svgEl("line", {
        x1: plotX0, x2: plotX1, y1: thresholdY, y2: thresholdY,
        style: "stroke: var(--color-axis-line); stroke-width: 1.5; stroke-linecap: round; stroke-dasharray: 5,4;",
        "class": "threshold-line"
      });
      svg.appendChild(tline);
      var tlabel = svgEl("text", {
        x: plotX1, y: thresholdLabelY, "text-anchor": "end",
        style: "fill: var(--color-text-secondary); font-size: " + AXIS_FONT_SIZE + ";",
        "class": "threshold-label"
      });
      // 言語は opts に焼き込まず毎回 t() で読む(ResizeObserver の古いクロージャ経由の
      // 再描画でも現在の言語設定を反映するため — opts.threshold.label に固定文字列を
      // 持たせると、リサイズ再描画時に初回描画時点の言語のまま固まってしまう)。
      tlabel.textContent = t().deThresholdLabel;
      svg.appendChild(tlabel);
    }

    if (hasSplit) {
      var boundaryLine = svgEl("line", {
        x1: boundaryX, x2: boundaryX, y1: plotY0, y2: plotY1,
        style: "stroke: var(--color-axis-line); stroke-width: 1; stroke-linecap: round; stroke-dasharray: 3,3;",
        "class": "verified-boundary"
      });
      svg.appendChild(boundaryLine);
      var boundaryLabel = svgEl("text", {
        x: boundaryX, y: (plotY0 - 6), "text-anchor": "middle",
        style: "fill: var(--color-text-secondary); font-size: " + AXIS_FONT_SIZE + ";",
        "class": "verified-boundary-label"
      });
      boundaryLabel.textContent = t().boundaryLabel;
      svg.appendChild(boundaryLabel);
    }

    var crosshairH = svgEl("line", {
      x1: plotX0, x2: plotX1, y1: 0, y2: 0,
      style: "stroke: var(--color-axis-line); stroke-width: 1;",
      visibility: "hidden", "class": "crosshair-h"
    });
    var crosshairV = svgEl("line", {
      x1: 0, x2: 0, y1: plotY0, y2: plotY1,
      style: "stroke: var(--color-axis-line); stroke-width: 1;",
      visibility: "hidden", "class": "crosshair-v"
    });
    svg.appendChild(crosshairH);
    svg.appendChild(crosshairV);

    var overlay = svgEl("rect", {
      x: plotX0, y: plotY0, width: (plotX1 - plotX0), height: (plotY1 - plotY0),
      style: "fill: transparent;", "class": "chart-overlay"
    });
    svg.appendChild(overlay);

    container.appendChild(svg);

    var tooltip = document.createElement("div");
    tooltip.className = "chart-tooltip";
    tooltip.style.display = "none";
    container.appendChild(tooltip);

    function nearestIndex(year) {
      var best = 0, bestDist = Infinity;
      for (var i = 0; i < years.length; i++) {
        var d = Math.abs(years[i] - year);
        if (d < bestDist) { bestDist = d; best = i; }
      }
      return best;
    }

    function handleMove(clientX, clientY) {
      var rect = svg.getBoundingClientRect();
      if (rect.width === 0) return;
      var scale = width / rect.width;
      var px = (clientX - rect.left) * scale;
      var year = years[0] + ((px - plotX0) / (plotX1 - plotX0)) * (years[years.length - 1] - years[0]);
      var idx = nearestIndex(year);
      var xPix = scaleX(years[idx]);
      var yPix = scaleY(quantiles.q50[idx]);

      crosshairV.setAttribute("x1", xPix);
      crosshairV.setAttribute("x2", xPix);
      crosshairV.setAttribute("visibility", "visible");
      crosshairH.setAttribute("y1", yPix);
      crosshairH.setAttribute("y2", yPix);
      crosshairH.setAttribute("visibility", "visible");

      var status = (untilYear !== null && untilYear !== undefined && years[idx] <= untilYear)
        ? t().verified : t().extrapolated;

      var lines = [
        t().tooltipYear(years[idx], status),
        ("q05: " + formatTooltipValue(quantiles.q05[idx], unit)),
        ("q25: " + formatTooltipValue(quantiles.q25[idx], unit)),
        ("q50: " + formatTooltipValue(quantiles.q50[idx], unit)),
        ("q75: " + formatTooltipValue(quantiles.q75[idx], unit)),
        ("q95: " + formatTooltipValue(quantiles.q95[idx], unit))
      ];
      tooltip.innerHTML = lines.map(function (line) {
        return ("<div>" + line.replace(/&/g, "&amp;").replace(/</g, "&lt;") + "</div>");
      }).join("");
      tooltip.style.display = "block";

      var containerRect = container.getBoundingClientRect();
      var relX = clientX - containerRect.left;
      var relY = clientY - containerRect.top;
      var tooltipLeft = relX + 12;
      if (containerRect.width > 0 && (tooltipLeft + 150) > containerRect.width) {
        tooltipLeft = relX - 150;
      }
      tooltip.style.left = tooltipLeft + "px";
      tooltip.style.top = Math.max(0, (relY - 20)) + "px";
    }

    function handleLeave() {
      crosshairV.setAttribute("visibility", "hidden");
      crosshairH.setAttribute("visibility", "hidden");
      tooltip.style.display = "none";
    }

    overlay.addEventListener("mousemove", function (ev) {
      handleMove(ev.clientX, ev.clientY);
    });
    overlay.addEventListener("mouseleave", handleLeave);
    overlay.addEventListener("touchstart", function (ev) {
      if (ev.touches && ev.touches.length > 0) {
        handleMove(ev.touches[0].clientX, ev.touches[0].clientY);
      }
      ev.preventDefault();
    }, { passive: false });
    overlay.addEventListener("touchmove", function (ev) {
      if (ev.touches && ev.touches.length > 0) {
        handleMove(ev.touches[0].clientX, ev.touches[0].clientY);
      }
      ev.preventDefault();
    }, { passive: false });
    overlay.addEventListener("touchend", handleLeave);

    attachResponsiveRedraw(container, function () {
      renderFanChart(container, opts);
    });
  }

  function renderCountryCharts(country, data) {
    var panel = document.querySelector('.country-panel[data-country="' + country + '"]');
    if (!panel || !data) return;
    var years = data.years;
    var untilYear = data.verified_horizon ? data.verified_horizon.until_calendar_year : null;

    VARIABLE_ORDER.forEach(function (name) {
      var entry = data.variables && data.variables[name];
      if (!entry) return;
      var container = panel.querySelector('.svg-container[data-variable="' + name + '"]');
      if (!container) return;
      var opts = {
        years: years,
        quantiles: { q05: entry.q05, q25: entry.q25, q50: entry.q50, q75: entry.q75, q95: entry.q95 },
        unit: entry.unit,
        untilYear: untilYear
      };
      if (name === "De") {
        opts.threshold = { value: 1 };
      }
      renderFanChart(container, opts);
    });

    var cfContainer = panel.querySelector('.svg-container[data-variable="count_forecast"]');
    if (cfContainer && data.count_forecast) {
      var cfYears = years.slice(1);
      var cf = data.count_forecast;
      renderFanChart(cfContainer, {
        years: cfYears,
        quantiles: { q05: cf.q05, q25: cf.q25, q50: cf.q50, q75: cf.q75, q95: cf.q95 },
        unit: cf.unit,
        untilYear: untilYear
      });
    }
  }

  function updateSummaryLine(country) {
    var meta = document.querySelector(
      '.country-meta[data-country="' + country + '"]'
    );
    var badge = document.getElementById("regime-badge");
    var summary = document.getElementById("forecast-summary");
    if (meta && badge && summary) {
      badge.textContent = t().regimePrefix + meta.getAttribute("data-regime");
      summary.textContent = t().summaryLine(
        meta.getAttribute("data-forecast-start-year"),
        meta.getAttribute("data-horizon-years"),
        meta.getAttribute("data-n")
      );
    }
  }

  function selectCountry(country) {
    var panels = document.querySelectorAll(".country-panel");
    panels.forEach(function (p) {
      p.hidden = p.getAttribute("data-country") !== country;
    });
    var tabs = document.querySelectorAll(".country-tab");
    tabs.forEach(function (tabEl) {
      var isSel = tabEl.getAttribute("data-country") === country;
      tabEl.setAttribute("aria-selected", isSel ? "true" : "false");
    });
    updateSummaryLine(country);
    if (window.location.hash !== "#" + country) {
      window.location.hash = country;
    }
    if (!renderedCountries[country]) {
      var data = readCountryData(country);
      if (data) {
        renderCountryCharts(country, data);
        renderedCountries[country] = true;
      }
    }
  }

  // #14: JA/EN トグル。html[data-lang] 属性の付替えは CSS(.i18n-ja/.i18n-en)が
  // 表示切替を担う。チャート SVG・サマリ行はテキストを直接生成するため、ここで
  // 現在表示中の国のみ再描画してテキストを追従させる。
  function applyLang(lang) {
    document.documentElement.setAttribute("data-lang", lang);
    document.documentElement.lang = lang;

    document.querySelectorAll(".lang-btn").forEach(function (btn) {
      var isSel = btn.getAttribute("data-lang-option") === lang;
      btn.setAttribute("aria-pressed", isSel ? "true" : "false");
    });

    document.querySelectorAll("[data-aria-ja][data-aria-en]").forEach(function (el) {
      el.setAttribute("aria-label", lang === "en" ? el.getAttribute("data-aria-en") : el.getAttribute("data-aria-ja"));
    });

    try {
      window.localStorage.setItem(LANG_STORAGE_KEY, lang);
    } catch (e) { /* localStorage 不可(private mode 等)— 保存は諦めて継続 */ }

    // 描画キャッシュを破棄する: 非表示国のチャートは旧言語の SVG ラベル
    // (境界・しきい値)を焼き込んだままなので、次に国切替されたとき
    // selectCountry の遅延描画パスで現在の言語で再描画させる。表示中の国だけ
    // ここで即時再描画する(ResizeObserver は幅不変だと再発火しないため、
    // キャッシュ破棄なしでは切替後も旧言語が残る)。
    renderedCountries = {};
    var activePanel = document.querySelector(".country-panel:not([hidden])");
    if (activePanel) {
      var country = activePanel.getAttribute("data-country");
      updateSummaryLine(country);
      var data = readCountryData(country);
      if (data) {
        renderCountryCharts(country, data);
        renderedCountries[country] = true;
      }
    }
  }

  function readStoredLang() {
    try {
      var stored = window.localStorage.getItem(LANG_STORAGE_KEY);
      if (stored === "ja" || stored === "en") return stored;
    } catch (e) { /* localStorage 不可 — 既定の ja にフォールバック */ }
    return "ja";
  }

  function init() {
    var manifest = readManifest();
    var countries = manifest.countries || [];
    if (countries.length === 0) return;

    var tabs = document.querySelectorAll(".country-tab");
    tabs.forEach(function (tab) {
      tab.addEventListener("click", function () {
        selectCountry(tab.getAttribute("data-country"));
      });
    });

    document.querySelectorAll(".lang-btn").forEach(function (btn) {
      btn.addEventListener("click", function () {
        applyLang(btn.getAttribute("data-lang-option"));
      });
    });

    var initial = window.location.hash
      ? window.location.hash.replace("#", "")
      : countries[0];
    if (countries.indexOf(initial) === -1) initial = countries[0];
    selectCountry(initial);

    applyLang(readStoredLang());
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", init);
  } else {
    init();
  }
})();
"""

# ============================================================
# 生成コア(純関数): 入力 JSON パス配列 → HTML 文字列
# ============================================================

"""
    build_dashboard_html(json_paths::Vector{String}) -> String

`json_paths` の各 M11 予報 JSON からデータ駆動の静的ダッシュボード HTML を生成する。
決定的(同一入力 → バイト同一の HTML)。自己完結(外部リクエストなし)。
チャート SVG の実描画は行わない(空コンテナのみ。Issue #12 のスコープ)。
"""
function build_dashboard_html(json_paths::Vector{String})
    isempty(json_paths) && error("json_paths が空です")

    countries = String[]
    docs = Dict{String,Any}()
    raw_bytes = Dict{String,String}()

    for path in json_paths
        isfile(path) || error("入力ファイルが見つかりません: $path")
        raw = read_raw_json_bytes(path)
        doc = JSON3.read(raw)
        country = String(doc[:country])
        validate_forecast_json(doc, path)
        push!(countries, country)
        docs[country] = doc
        raw_bytes[country] = raw
    end

    manifest_json = JSON3.write(Dict("countries" => countries))

    data_scripts = join([
        """<script type="application/json" id="data-$(html_escape(c))">$(raw_bytes[c])</script>"""
        for c in countries
    ], "\n    ")

    country_sections = join([render_country_section(c, docs[c]) for c in countries], "\n")

    header = render_header(countries, docs)

    legend_1_i18n = i18n_span(
        "90% 区間(q05–q95)/ 50% 区間(q25–q75)/ 中央値(q50)",
        "90% interval (q05–q95) / 50% interval (q25–q75) / median (q50)",
    )
    legend_2_i18n = i18n_span(
        "検証済み予報(起点+6年、#0052/#0069) と 外挿領域(不確実性の広がり)の区別(不透明度・背景ウォッシュ・境界ラベルで表現)",
        "Distinction between verified forecast (origin + 6 years, #0052/#0069) and the extrapolated region (spread of uncertainty), expressed via opacity, background wash, and boundary label",
    )

    html = """<!doctype html>
<html lang="ja">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>未来予報 — 社会動態シミュレータ 将来予報</title>
<style>
$DASHBOARD_CSS
</style>
<noscript><style>.country-panel[hidden] { display: block; }</style></noscript>
</head>
<body>
$header
<footer class="legend-bar">
  <div class="legend-item">
    <svg class="legend-swatch" width="48" height="16" viewBox="0 0 48 16" aria-hidden="true">
      <rect x="0" y="1" width="48" height="14" style="fill: var(--color-series); fill-opacity: 0.16;"></rect>
      <rect x="0" y="4" width="48" height="8" style="fill: var(--color-series); fill-opacity: 0.32;"></rect>
      <line x1="0" y1="8" x2="48" y2="8" style="stroke: var(--color-series); stroke-width: 2;"></line>
    </svg>
    <span>$legend_1_i18n</span>
  </div>
  <div class="legend-item">
    <svg class="legend-swatch" width="48" height="16" viewBox="0 0 48 16" aria-hidden="true">
      <rect x="0" y="1" width="24" height="14" style="fill: var(--color-series); fill-opacity: 0.32;"></rect>
      <rect x="24" y="1" width="24" height="14" style="fill: var(--color-gridline); fill-opacity: 0.5;"></rect>
      <rect x="24" y="1" width="24" height="14" style="fill: var(--color-series); fill-opacity: 0.176;"></rect>
      <line x1="24" y1="0" x2="24" y2="16" style="stroke: var(--color-axis-line); stroke-width: 1; stroke-dasharray: 3,3;"></line>
    </svg>
    <span>$legend_2_i18n</span>
  </div>
</footer>
<main>
$country_sections
$(DISCLAIMER_TEXT)
</main>
<script type="application/json" id="manifest">$manifest_json</script>
$data_scripts
<script>
$DASHBOARD_JS
</script>
</body>
</html>
"""
    return html
end

# ============================================================
# CLI ラッパ
# ============================================================

function default_output_path()
    return joinpath("docs", "dashboard", "index.html")
end

function resolve_input_path(iso3::AbstractString)
    return joinpath("experiments", "output", "M11_forecast_$(iso3).json")
end

function print_regeneration_help(iso3::AbstractString, path::AbstractString)
    println(stderr, "入力ファイルが見つかりません: $path")
    println(stderr, "以下のコマンドで再生成してください(FORECAST_JSON.md L8-19 参照):")
    println(stderr, "  julia --project=experiments -t 8 experiments/M11_forecast.jl $(iso3)")
end

function main(args::Vector{String})
    iso3_list = String[]
    out_path = default_output_path()

    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--out"
            i == length(args) && error("--out には値が必要です")
            out_path = args[i+1]
            i += 2
        else
            push!(iso3_list, a)
            i += 1
        end
    end

    if isempty(iso3_list)
        iso3_list = ["JPN", "THA", "KOR", "TUR", "EGY"]
    end

    json_paths = String[]
    for iso3 in iso3_list
        path = resolve_input_path(iso3)
        if !isfile(path)
            print_regeneration_help(iso3, path)
            exit(1)
        end
        push!(json_paths, path)
    end

    html = build_dashboard_html(json_paths)

    out_dir = dirname(out_path)
    if !isempty(out_dir) && !isdir(out_dir)
        mkpath(out_dir)
    end
    write(out_path, html)
    println("生成完了: $out_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
