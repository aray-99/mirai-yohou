# M12 予報円ダッシュボード 設計仕様(凍結: DECISIONS #0073)

M11 将来予報 JSON(契約 = docs/FORECAST_JSON.md)を可視化する自己完結静的 HTML
ダッシュボードの実装仕様。実装 Issue はすべて本書を規範とする。上位の枠組みは
docs/PHASE3_DESIGN.md §3(#0070 凍結)。

- UI 言語: **日本語主体 + 英語変数名**(2026-07-18 ユーザー確認)
- 公開方式: **main の /docs 配下**(`docs/dashboard/index.html`、Pages ソース = main /docs)
  (2026-07-18 ユーザー確認)

## 1. 目的と義務要件

M11 の予報 JSON(JPN/THA、horizon 30 年)をブラウザだけで閲覧できる形にする。
科学的誠実さの線引き(#0070)を UI として実装することが本体である。

**義務要件(凍結。受入判定で個別に確認する):**

- (a) **検証済み/外挿の視覚区別**(#0070、FORECAST_JSON.md L169-178):
  `verified_horizon.until_calendar_year` を境界として、検証済み区間(起点+6年)と
  外挿領域を全パネルで視覚的に区別する。外挿領域を「予報」と表記しない
  (「外挿領域(不確実性の広がり)」と表記する)。
- (b) **p_ex_fallback 警告表示**(#0071、FORECAST_JSON.md L75-84):
  `jump_thinning.p_ex_fallback == "no_count_data"` のとき、count_forecast パネルに
  警告バッジ「⚠ ジャンプリスク評価不能(同化窓内にカウントデータなし)」を表示し、
  通常の p_ex > 0 ケースと区別する。現行 JPN/THA は null(バッジ非表示)が正常系だが、
  M13 横展開国で発動しうるため**実装とテストは必須**。

**凍結済み完了条件**(#0070。M12 クローズ判定 = fable-required):

1. JPN/THA の予報円が国切替つきで表示される
2. 来歴(コミット SHA・シード・データ取得日)が UI で確認できる
3. GitHub Pages URL で閲覧できる

## 2. アーキテクチャ

### 2.1 生成器

- **`experiments/M12_dashboard.jl`** 1本(#0070「Julia スクリプト1本で決定的再生成」)。
  依存は experiments 環境の既存パッケージ(JSON3 等)のみ。新規パッケージ追加不可。
- CLI: `julia --project=experiments experiments/M12_dashboard.jl JPN THA [--out docs/dashboard/index.html]`
  - 引数の各 ISO3 について `experiments/output/M11_forecast_{ISO3}.json` を読む
    (無ければ再生成コマンドを表示してエラー終了。FORECAST_JSON.md L8-19 参照)。
  - `--out` 省略時は `docs/dashboard/index.html`。出力ディレクトリは無ければ作成。
- HTML/CSS/JS テンプレートは M12_dashboard.jl 内の文字列リテラルとして保持
  (別ファイルに分けない — 「スクリプト1本」の凍結文言に従う)。生成器がやることは
  (i) JSON 検証(トップレベルキーの存在確認)、(ii) 国ごとの JSON をそのまま
  `<script type="application/json" id="data-{ISO3}">` として埋め込み、(iii) 国リスト等の
  マニフェストを `<script type="application/json" id="manifest">` に埋め込み。
  **描画ロジックは持たない**(描画はクライアント側 JS)。

### 2.2 決定性の規約

- **同一入力 JSON → バイト同一の HTML**。生成器内で `now()`・乱数・環境依存値
  (ホスト名、絶対パス)を出力に含めることを禁止。日時・シード・コミット SHA は
  すべて入力 JSON の `provenance` から転記する。
- 埋め込み JSON は入力ファイルのバイト列をそのまま挿入する(再シリアライズしない —
  キー順・浮動小数表現の揺れを避ける)。`</script>` 系のエスケープ問題は入力が
  信頼できる自作 JSON であるため、`<` を含まないことの検証のみ行う(含めばエラー)。
- 生成 HTML は git コミットする(ダッシュボード自体が来歴内蔵 — PHASE3_PROPOSAL §4.2)。
  `experiments/output/` は従来通り gitignore。

### 2.3 自己完結性

- 外部リクエストゼロ: CDN・Web フォント・外部画像・fetch 禁止。CSS/JS はすべて
  インライン。フォントはシステムサンズ
  (`system-ui, -apple-system, "Segoe UI", sans-serif`)。
- Claude Artifact の CSP 制約下でもプレビュー可能であること(開発時の確認手段)。

## 3. 画面構成

単一ページ。上から:

1. **ヘッダ**: タイトル「未来予報 — 社会動態シミュレータ 将来予報」、国切替タブ
   (JPN / THA — manifest から生成、M13 で国が増える前提のループ実装)、
   レジームバッジ(`regime`: stable / volatile — 中立色 + テキスト、状態色は使わない)、
   予報起点(`forecast_start.calendar_year`)・ホライズン(`horizon_years` 年)・
   アンサンブル数 `N` の一行サマリ。
2. **凡例バー**(全パネル共通、1箇所のみ): 90% 区間(q05–q95)/ 50% 区間(q25–q75)/
   中央値(q50)の帯サンプル + 「検証済み予報(起点+6年、#0052/#0069)」「外挿領域
   (不確実性の広がり)」の区別サンプル。
3. **9 実系列パネル**(グリッド、`variables` の P, w, y, g_swiid, T_proxy, phi, v, tau, p):
   パネル題は「英語変数名 — 日本語ラベル(単位)」。日本語ラベルは
   P=人口、w=労働分配率相当、y=一人当たり GDP、g_swiid=所得格差(可処分ジニ)、
   T_proxy=技術(特許出願)、phi=政治参加、v=情報流通、tau=社会的信頼、p=政治的自由。
   単位は JSON の `unit` を併記。y=0 が起点断面(予報積分前のアンサンブル事後)である
   旨をパネル群の脚注に1回記載。
4. **count_forecast パネル**: 「期待騒乱イベント数 / 年」。x は予報年 y=1..horizon
   (配列長 horizon — variables と1つずれる点に注意、FORECAST_JSON.md L140-157)。
   `nu_star`・`r_hat`(null なら「Poisson」)を脚注表示。**義務要件 (b) のバッジは
   このパネルのヘッダに置く**(status 色 + ⚠ アイコン + テキスト。色のみで意味を
   伝えない)。
5. **De 診断パネル**: 「De — 体質診断(無次元)」。De=1 のしきい線(破線 + ラベル
   「De=1: 安定/変動レジーム境界」)。外挿領域の主提示物として、パネル脚注に
   「外挿領域では個別値でなくレジーム傾向(De の帯が 1 を跨ぐか)を読む」旨を記載。
6. **来歴パネル**(`provenance` + トップレベルから): commit(短縮 SHA + 全文 title 属性)、
   seed、generated_at、data_fetched_at、frozen_config、script、design_decision、
   frozen_decisions、finite_check.pass、endogenous_jumps_total、
   assimilation の要約(n_obs、window、ess_range、nu_star)。国切替に追従する。
7. **免責注記**(FORECAST_JSON.md L180-189 を日本語で転記): 外挿領域(y>6)の分位帯は
   モデル内生力学のみに基づき観測制約を受けないこと、p_ex フォールバック発動国の
   count_forecast は過小評価の可能性があること、検証の実証裏付けは M9/M10
   walk-forward(#0052/#0069)であること。

## 4. チャート仕様(ファンチャート = 予報円)

### 4.1 幾何

- SVG をクライアント側 vanilla JS で構築(viewBox 固定、width 100% でレスポンシブ)。
- x 軸: 西暦年(`years`)。線形。目盛は 5 年刻み。
- y 軸: 線形、範囲は [min(q05), max(q95)] に 5% パディング。目盛 4〜5 本、
  大きい数値は SI 略記(1.2M 等)、fraction は小数 2 桁。
- 帯: q05–q95 ポリゴン(不透明度 0.16)、q25–q75 ポリゴン(0.32)、q50 折れ線
  (太さ 2px)。全パネル同一色相(下記 series 色)— パネル間で色を変えない
  (色は系列の同一性を運ぶ。全パネルが同一アンサンブルの分位なので単一色相が正)。
- **検証済み/外挿の区別(義務要件 (a) の実装)**: `until_calendar_year` に境界縦線
  (破線、baseline 色)+ 直上ラベル「検証済み | 外挿」。外挿領域は
  (i) 背景に極薄ウォッシュ(gridline 色、不透明度 0.5)、(ii) 帯・中央値線の不透明度を
  0.55 倍。境界は clip-path で帯を 2 分割して実装する(色相は変えない —
  区別は明度・背景・ラベルの3チャネルで冗長化)。
- グリッド線はヘアライン(gridline 色)、軸・ラベルは muted 色。データより前に出ない。

### 4.2 色(dataviz 参照パレット。ライト/ダーク両対応)

CSS カスタムプロパティで定義し、`@media (prefers-color-scheme: dark)` で切替:

| ロール | Light | Dark |
|---|---|---|
| ページ地 | `#f9f9f7` | `#0d0d0d` |
| チャート面 | `#fcfcfb` | `#1a1a19` |
| 主文字 | `#0b0b0b` | `#ffffff` |
| 副文字 | `#52514e` | `#c3c2b7` |
| 軸・目盛ラベル | `#898781` | `#898781` |
| グリッド線 | `#e1e0d9` | `#2c2c2a` |
| 軸線・境界縦線 | `#c3c2b7` | `#383835` |
| series(帯・中央値) | `#2a78d6` | `#3987e5` |
| 警告バッジ (b) | `#fab219` | `#fab219` |

- series は単一色相(青)のみ使用 — 検証済みパレットの categorical slot 1。
  複数色相を使わないため CVD 検証は自明に通る(区別は色でなく不透明度 + 位置)。
- 警告バッジは status "warning"。ライト面ではコントラスト不足のため
  **⚠ アイコン + 黒文字ラベルを必ず併記**(色のみ禁止)。
- 文字は常に文字色トークン(series 色の文字を書かない)。

### 4.3 インタラクション

- 各パネルにホバー十字線 + ツールチップ: 年、q05/q25/q50/q75/q95 の値、
  検証済み/外挿の別。タッチ環境ではタップで同等表示。
- 国切替: タブクリックで全パネルを対象国の埋め込み JSON から再描画。
  URL ハッシュ(`#THA`)に反映し、ロード時に復元(デフォルトは manifest 先頭)。
- 各パネルに `<details>` の**表ビュー**(年 × 5 分位の table)— アクセシビリティ
  要件。JS 無効環境でも表ビューと来歴・免責テキストは読めること
  (SVG は JS 必須で可)。

## 5. テスト(Julia テストスイートに追加)

フィクスチャ: 実出力を縮約した小 JSON(horizon 5 年・2 か国相当、
`p_ex_fallback: "no_count_data"` の国を 1 つ含む)を `test/fixtures/` にコミット。

1. **生成成功**: フィクスチャ入力で HTML が生成され、サイズ > 0、
   `<script type="application/json"` が国数分存在。
2. **決定性**: 同一入力で 2 回生成しバイト一致。
3. **義務要件マーカー**: 生成 HTML に (a) 境界ラベル文字列「外挿」と
   until_calendar_year の値、(b) fallback 国向け警告バッジ文字列
   「ジャンプリスク評価不能」が含まれる(非 fallback 国のみの生成では (b) が
   **含まれない**ことも検査 — バッジはデータ駆動であること)。
4. **埋め込み忠実性**: 埋め込み JSON 断片が入力ファイルのバイト列と一致。

(ブラウザ描画の自動テストは行わない — 目視確認は受入判定(Issue C)で実施。)

## 6. スコープ外(M12 では行わない)

- 観測実績値の重ね描き(予報 JSON に観測系列が含まれないため。必要なら M13 以降で
  スキーマ拡張とセットで検討)
- 国別比較ビュー・アニメーション・サーバサイド機能
- スキーマ(FORECAST_JSON.md)の変更

## 7. 受入判定(Issue C、fable-required)

1. 完了条件 1〜3(§1)を実 URL で確認
2. 義務要件 (a)(b) を DOM / 表示で確認((b) はフィクスチャ HTML で確認)
3. 決定性: `M11_forecast_{JPN,THA}.json` から 2 回生成してバイト一致
4. フルテスト green(`julia --project -e 'using Pkg; Pkg.test()'`)
5. Artifact プレビューで表示崩れ・ラベル衝突の目視確認(ライト/ダーク両方)
