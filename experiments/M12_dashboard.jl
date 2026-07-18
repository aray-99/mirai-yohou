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

# 設計書 §3-3 の正典対応表(英語変数名 — 日本語ラベル)。一字一句厳守。
const VARIABLE_LABELS = [
    ("P", "総人口"),
    ("w", "生産年齢人口比率(15-64歳)"),
    ("y", "一人当たり産出(実質GDP/capita、合成観測)"),
    ("g_swiid", "格差指標(可処分ジニ)"),
    ("T_proxy", "技術フロンティア(特許出願が観測代理)"),
    ("phi", "技術普及率(ネット利用率が観測代理)"),
    ("v", "情報伝播速度(モバイル契約が観測代理)"),
    ("tau", "制度信頼(WVS政府信頼)"),
    ("p", "政治的分極度(V-Dem v2cacamps、高いほど分裂)"),
]

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
    for (name, _) in VARIABLE_LABELS
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

    io = IOBuffer()
    print(io, """
    <header class="dashboard-header">
      <h1>未来予報 — 社会動態シミュレータ 将来予報</h1>
      <nav class="country-tabs" id="country-tabs" role="tablist" aria-label="国切替">
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
    return """<p class="panel-footnote">y=0 は予報起点断面(予報積分前のアンサンブル事後分布)であり、y=1 以降が実際の予報値である。</p>"""
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
    last_obs_str = last_obs === nothing ? "—(状態変数でない)" : string(last_obs)
    return """
        <details class="table-view">
          <summary>表ビュー($(html_escape(varname)))</summary>
          <p class="table-meta">単位: $unit / 最終観測年: $last_obs_str</p>
          <table>
            <thead><tr><th>年</th><th>q05</th><th>q25</th><th>q50</th><th>q75</th><th>q95</th></tr></thead>
            <tbody>
            $(String(take!(rows)))
            </tbody>
          </table>
        </details>
    """
end

function render_variable_panel(country::String, varname::String, label::String, entry, years)
    unit = string(entry[:unit])
    transform = string(entry[:transform])
    table = render_variable_table(country, varname, entry, years)
    return """
      <section class="panel variable-panel" data-country="$(html_escape(country))" data-variable="$(html_escape(varname))">
        <h3>$(html_escape(varname)) — $(html_escape(label))($(html_escape(unit)))</h3>
        <p class="panel-transform">変換: $(html_escape(transform))</p>
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
        badge = """<span class="warning-badge" role="status">⚠ ジャンプリスク評価不能(同化窓内にカウントデータなし)</span>"""
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

    return """
      <section class="panel count-forecast-panel" data-country="$(html_escape(country))">
        <h3>期待騒乱イベント数 / 年 $badge</h3>
        <p class="panel-note">単位: $(html_escape(string(get(cf, :unit, ""))))</p>
        <p class="panel-footnote">nu_star = $nu_star / r_hat = $r_hat_str</p>
        <div class="svg-container" data-role="chart" data-variable="count_forecast"></div>
        <details class="table-view">
          <summary>表ビュー(count_forecast)</summary>
          <table>
            <thead><tr><th>年</th><th>q05</th><th>q25</th><th>q50</th><th>q75</th><th>q95</th></tr></thead>
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
    return """
      <section class="panel de-panel" data-country="$(html_escape(country))">
        <h3>De — 体質診断(無次元)</h3>
        <p class="panel-footnote">De=1: 安定/変動レジーム境界。外挿領域では個別値でなくレジーム傾向(De の帯が 1 を跨ぐか)を読む。</p>
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

    return """
      <section class="panel provenance-panel" data-country="$(html_escape(country))">
        <h3>来歴</h3>
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

const DISCLAIMER_TEXT = """
      <section class="panel disclaimer-panel">
        <h3>免責注記</h3>
        <p>
          外挿領域(検証済みホライズンより先、y &gt; 6)の分位帯は、モデルの内生力学
          (Hawkes 自己励起・OU 型緩和・L3 拡大パラメータの継続ランダムウォーク)のみに
          基づく統計的投影であり、観測による制約を一切受けていない。
        </p>
        <p>
          特に <code>p_ex_fallback: "no_count_data"</code> が発動している国・オリジンでは、
          ジャンプ関連の力学(Γ 適用による急激な社会的応力上昇)が一切表現されないため、
          <code>count_forecast</code> は「ジャンプが起きない前提」の過小評価である
          可能性がある。
        </p>
        <p>
          検証の実証裏付けは M9/M10 の expanding-window walk-forward
          (#0052 設計・#0069 合格判定)であり、検証済み範囲は起点+6年に限られる。
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
    for (varname, label) in VARIABLE_LABELS
        entry = doc[:variables][Symbol(varname)]
        print(io, render_variable_panel(country, varname, label, entry, doc[:years]))
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
:root {
  --color-page-bg: #f9f9f7;
  --color-panel-bg: #fcfcfb;
  --color-text-primary: #0b0b0b;
  --color-text-secondary: #52514e;
  --color-axis-label: #898781;
  --color-gridline: #e1e0d9;
  --color-axis-line: #c3c2b7;
  --color-series: #2a78d6;
  --color-warning: #fab219;
}
@media (prefers-color-scheme: dark) {
  :root {
    --color-page-bg: #0d0d0d;
    --color-panel-bg: #1a1a19;
    --color-text-primary: #ffffff;
    --color-text-secondary: #c3c2b7;
    --color-axis-label: #898781;
    --color-gridline: #2c2c2a;
    --color-axis-line: #383835;
    --color-series: #3987e5;
    --color-warning: #fab219;
  }
}
* { box-sizing: border-box; }
body {
  margin: 0;
  font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
  background: var(--color-page-bg);
  color: var(--color-text-primary);
}
.dashboard-header {
  padding: 1rem 1.25rem;
  border-bottom: 1px solid var(--color-gridline);
}
.dashboard-header h1 { font-size: 1.25rem; margin: 0 0 0.5rem 0; }
.country-tabs { display: flex; gap: 0.5rem; margin-bottom: 0.5rem; }
.country-tab {
  font: inherit;
  padding: 0.35rem 0.9rem;
  border: 1px solid var(--color-axis-line);
  border-radius: 999px;
  background: var(--color-panel-bg);
  color: var(--color-text-primary);
  cursor: pointer;
}
.country-tab[aria-selected="true"] { border-color: var(--color-series); font-weight: 600; }
.summary-line { color: var(--color-text-secondary); font-size: 0.9rem; }
.regime-badge {
  display: inline-block;
  padding: 0.1rem 0.6rem;
  border: 1px solid var(--color-axis-line);
  border-radius: 4px;
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
  border-radius: 6px;
  padding: 0.75rem 1rem;
  margin-bottom: 1rem;
}
.panel h3 { margin: 0 0 0.35rem 0; font-size: 1rem; }
.panel-transform, .panel-note, .panel-footnote, .table-meta {
  color: var(--color-text-secondary);
  font-size: 0.8rem;
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
  border-radius: 4px;
  padding: 0.3rem 0.5rem;
  font-size: 0.72rem;
  line-height: 1.3;
  color: var(--color-text-primary);
  box-shadow: 0 1px 4px rgba(0, 0, 0, 0.18);
  z-index: 5;
  white-space: nowrap;
}
.warning-badge {
  display: inline-block;
  background: var(--color-warning);
  color: #0b0b0b;
  border-radius: 4px;
  padding: 0.1rem 0.5rem;
  font-size: 0.85rem;
  margin-left: 0.5rem;
}
table { border-collapse: collapse; width: 100%; font-size: 0.8rem; margin-top: 0.4rem; }
th, td { border: 1px solid var(--color-gridline); padding: 0.2rem 0.4rem; text-align: right; }
th:first-child, td:first-child { text-align: left; }
details.table-view summary { cursor: pointer; color: var(--color-text-secondary); font-size: 0.85rem; }
.provenance-list { display: grid; grid-template-columns: max-content 1fr; gap: 0.15rem 0.75rem; font-size: 0.82rem; }
.provenance-list dt { color: var(--color-text-secondary); }
.provenance-list dd { margin: 0; word-break: break-word; }
.disclaimer-panel p { font-size: 0.85rem; line-height: 1.5; }
footer.legend-bar { padding: 0.5rem 1.25rem; color: var(--color-text-secondary); font-size: 0.8rem; }
.legend-item { display: flex; align-items: center; gap: 0.4rem; margin: 0.2rem 0; }
.legend-swatch { flex: none; vertical-align: middle; }
"""

const DASHBOARD_JS = """
(function () {
  "use strict";

  var VARIABLE_ORDER = ["P", "w", "y", "g_swiid", "T_proxy", "phi", "v", "tau", "p", "De"];

  var CHART_WIDTH = 640;
  var CHART_HEIGHT = 260;
  var MARGIN = { top: 20, right: 16, bottom: 32, left: 60 };

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

  function makeScaleX(years) {
    var y0 = years[0], y1 = years[years.length - 1];
    var x0 = MARGIN.left, x1 = CHART_WIDTH - MARGIN.right;
    var span = (y1 - y0) || 1;
    return function (year) {
      return x0 + ((year - y0) / span) * (x1 - x0);
    };
  }

  function makeScaleY(domain) {
    var y0 = MARGIN.top, y1 = CHART_HEIGHT - MARGIN.bottom;
    var span = (domain[1] - domain[0]) || 1;
    return function (v) {
      return y1 - ((v - domain[0]) / span) * (y1 - y0);
    };
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

    var extraForDomain = threshold ? [threshold.value] : [];
    var domain = computeYDomain(quantiles, extraForDomain);
    var scaleX = makeScaleX(years);
    var scaleY = makeScaleY(domain);

    var svg = svgEl("svg", {
      viewBox: "0 0 " + CHART_WIDTH + " " + CHART_HEIGHT,
      width: "100%",
      height: "auto",
      preserveAspectRatio: "xMidYMid meet",
      role: "img",
      "aria-label": "予報ファンチャート",
      "class": "chart-svg"
    });

    var plotX0 = MARGIN.left, plotX1 = CHART_WIDTH - MARGIN.right;
    var plotY0 = MARGIN.top, plotY1 = CHART_HEIGHT - MARGIN.bottom;

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

    var yTickValues = linspace(domain[0], domain[1], 5);
    yTickValues.forEach(function (v) {
      var yPix = scaleY(v);
      var line = svgEl("line", {
        x1: plotX0, x2: plotX1, y1: yPix, y2: yPix,
        style: "stroke: var(--color-gridline); stroke-width: 1;",
        "class": "chart-gridline"
      });
      svg.appendChild(line);
      var label = svgEl("text", {
        x: plotX0 - 6, y: yPix, "text-anchor": "end", "dominant-baseline": "middle",
        style: "fill: var(--color-axis-label); font-size: 9px;",
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
        style: "fill: var(--color-axis-label); font-size: 9px;",
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
        style: ("fill: none; stroke: var(--color-series); stroke-width: 2; stroke-opacity: " + opacityFactor + ";"),
        "class": "chart-median-line"
      }));
      return g;
    }

    svg.appendChild(buildDataGroup((uid + "-verified"), 1));
    if (hasSplit) {
      svg.appendChild(buildDataGroup((uid + "-extrapolated"), 0.55));
    }

    if (threshold) {
      var ty = scaleY(threshold.value);
      var tline = svgEl("line", {
        x1: plotX0, x2: plotX1, y1: ty, y2: ty,
        style: "stroke: var(--color-axis-line); stroke-width: 1.5; stroke-dasharray: 5,4;",
        "class": "threshold-line"
      });
      svg.appendChild(tline);
      var tlabel = svgEl("text", {
        x: plotX1, y: (ty - 4), "text-anchor": "end",
        style: "fill: var(--color-text-secondary); font-size: 9px;",
        "class": "threshold-label"
      });
      tlabel.textContent = threshold.label;
      svg.appendChild(tlabel);
    }

    if (hasSplit) {
      var boundaryLine = svgEl("line", {
        x1: boundaryX, x2: boundaryX, y1: plotY0, y2: plotY1,
        style: "stroke: var(--color-axis-line); stroke-width: 1; stroke-dasharray: 3,3;",
        "class": "verified-boundary"
      });
      svg.appendChild(boundaryLine);
      var boundaryLabel = svgEl("text", {
        x: boundaryX, y: (plotY0 - 6), "text-anchor": "middle",
        style: "fill: var(--color-text-secondary); font-size: 9px;",
        "class": "verified-boundary-label"
      });
      boundaryLabel.textContent = "検証済み | 外挿";
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
      var scale = CHART_WIDTH / rect.width;
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
        ? "検証済み" : "外挿";

      var lines = [
        (String(years[idx]) + "年(" + status + ")"),
        ("q05: " + formatTooltipValue(quantiles.q05[idx], unit)),
        ("q25: " + formatTooltipValue(quantiles.q25[idx], unit)),
        ("q50: " + formatTooltipValue(quantiles.q50[idx], unit)),
        ("q75: " + formatTooltipValue(quantiles.q75[idx], unit)),
        ("q95: " + formatTooltipValue(quantiles.q95[idx], unit))
      ];
      tooltip.innerHTML = lines.map(function (t) {
        return ("<div>" + t.replace(/&/g, "&amp;").replace(/</g, "&lt;") + "</div>");
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
        opts.threshold = { value: 1, label: "De=1: 安定/変動レジーム境界" };
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

  function selectCountry(country) {
    var panels = document.querySelectorAll(".country-panel");
    panels.forEach(function (p) {
      p.hidden = p.getAttribute("data-country") !== country;
    });
    var tabs = document.querySelectorAll(".country-tab");
    tabs.forEach(function (t) {
      var isSel = t.getAttribute("data-country") === country;
      t.setAttribute("aria-selected", isSel ? "true" : "false");
    });
    var meta = document.querySelector(
      '.country-meta[data-country="' + country + '"]'
    );
    var badge = document.getElementById("regime-badge");
    var summary = document.getElementById("forecast-summary");
    if (meta && badge && summary) {
      badge.textContent = "regime: " + meta.getAttribute("data-regime");
      summary.textContent =
        "予報起点 " + meta.getAttribute("data-forecast-start-year") +
        " / ホライズン " + meta.getAttribute("data-horizon-years") + " 年" +
        " / N=" + meta.getAttribute("data-n");
    }
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

    var initial = window.location.hash
      ? window.location.hash.replace("#", "")
      : countries[0];
    if (countries.indexOf(initial) === -1) initial = countries[0];
    selectCountry(initial);
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
    <span>90% 区間(q05–q95)/ 50% 区間(q25–q75)/ 中央値(q50)</span>
  </div>
  <div class="legend-item">
    <svg class="legend-swatch" width="48" height="16" viewBox="0 0 48 16" aria-hidden="true">
      <rect x="0" y="1" width="24" height="14" style="fill: var(--color-series); fill-opacity: 0.32;"></rect>
      <rect x="24" y="1" width="24" height="14" style="fill: var(--color-gridline); fill-opacity: 0.5;"></rect>
      <rect x="24" y="1" width="24" height="14" style="fill: var(--color-series); fill-opacity: 0.176;"></rect>
      <line x1="24" y1="0" x2="24" y2="16" style="stroke: var(--color-axis-line); stroke-width: 1; stroke-dasharray: 3,3;"></line>
    </svg>
    <span>検証済み予報(起点+6年、#0052/#0069) と 外挿領域(不確実性の広がり)の区別(不透明度・背景ウォッシュ・境界ラベルで表現)</span>
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
        iso3_list = ["JPN", "THA"]
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
