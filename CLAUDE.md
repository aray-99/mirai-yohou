# CLAUDE.md — 未来予報(社会動態シミュレータ)

Julia/SciML 製の社会動態シミュレータ。13次元 SDE + Hawkes ジャンプ + EnKF/EnKS 同化。
原典仕様は docs/SPEC.md、全決定事項は docs/DECISIONS.md(連番エントリが一次資料)、
Phase 2 設計は docs/PHASE2_DESIGN.md。

## エージェント運用方針(トークン最適化)

**メインセッション(Fable)はオーケストレーターに徹する。**

- 実装・コード探索・テスト実行・ログ調査・パラメータスイープ等の作業タスクは
  **Sonnet のサブエージェント**(Agent ツール、`model: "sonnet"`)に委譲する。
  読み取り専用の探索は Explore、実装や実行を伴うものは general-purpose を使う。
- サブエージェントへの指示は「目的・対象ファイル・完了条件・報告してほしい形式」を
  明示した自己完結プロンプトにする(サブエージェントは本ファイルを読めるが
  会話履歴は見えない)。
- **Fable(メイン)に吸い上げるもの**: 科学的な意思決定(基準・プロトコル・
  数式変更)、DECISIONS エントリの起草、トラブルの診断方針、ユーザーへの
  判断依頼。数値結果の解釈と判定は Fable が行う。
- docs/SPEC.md・DECISIONS の全文読み込みをメインで行わない。必要箇所の抜粋を
  サブエージェントに取らせる。

## 友人委譲(GitHub Actions + @claude、#0076)

- 友人(GitHub: Rodan-n)への実装委譲は **Issue 起票 + `@claude` メンション**で行う。
  `.github/workflows/claude.yml` が claude-code-action@v1 を起動し、友人のクォータ
  (secrets.CLAUDE_CODE_OAUTH_TOKEN)で実装 → PR 作成まで自動実行される。
  メッセージアプリでの人間中継は廃止。
- 起票する Issue 本文は**自己完結プロンプト**にする(目的・対象ファイル・完了条件・
  凍結ファイル不可侵・テスト方法を明記。Action 側は会話履歴を見えない)。
- 起動は aray-99 / Rodan-n の発話のみ(第三者による友人クォータ消費の防止)。
  無関係な Issue コメントに `@claude` を書かないこと。
- Actions ランナーには Julia 環境がないため、**フルテストはオーナーがレビュー時に
  実行**する(確立済み手順: checkout → diff 精査 → テスト自走 → --no-ff マージ)。
- PR の自動監視ループは回さない(ユーザーが GitHub 通知で検知し、指示を受けて
  Fable がレビューする)。

## セッション引き継ぎ(常時運用)

- **セッション開始時は必ず HANDOVER_PROMPT.md(リポジトリ直下、未追跡)を読んで
  そこから再開する**。ユーザーからの明示指示を待たない。
- **状況が変わるたびに HANDOVER_PROMPT.md を即座に更新する**(マイルストーン
  進捗、PR レビュー状態、実行中のラン、未完のタスク)。長時間ラン起動前・
  ターン終了前が更新チェックポイント。いつ /clear されても次セッションが
  そのまま引き継げる状態を常に保つこと。

## 意思決定権限(#0039、2026-07-10 ユーザー委任)

- **手法・是正方針の選択は Fable が自律決定してよい**(個別のユーザー承認不要)。
  ただし決定と理由は必ず DECISIONS エントリとして**事前に**記録する。
- 例外(引き続き明示のユーザー承認が必要): 凍結済み合格基準・しきい値の変更。
  結果を見た後に基準を動かさない原則は不変。
- マイルストーンの合格認定・main マージ・タグ付与は結果とともに事後報告。

## 開発運用ルール(SPEC §0.5 の要点。詳細は SPEC 参照)

- Conventional Commits。scope は coordinates/parameters/drift/diffusion/diagnostics/
  jumps/integrator/observation/enkf/weights/experiments/docs/ci から。
  **件名は英語**(`type(scope): summary`)、本文は日本語可。DECISIONS 参照は
  件名でなく本文末尾の `Refs: docs/DECISIONS.md#XXXX` 行に書く(件名の
  `(#XXXX)` は GitHub の issue/PR 参照と誤読されるため禁止)。
- feature/<milestone>-<topic> ブランチ → develop へ `--no-ff` マージ。マージ body に
  マイルストーン番号・通過テスト名・DECISIONS ID を必ず記載。main へは
  マイルストーン完了時のみ(vM* タグ)。
- **テスト red でコミットしない**(`julia --project -e 'using Pkg; Pkg.test()'`、
  約3分、299本)。
- 合格基準・凍結済み設定値の変更は**ユーザーの明示承認 + DECISIONS エントリ必須**。
  基準は結果を見た後に動かさない(不合格なら診断を記録して報告)。
- 実験成果物には来歴(コミット SHA・シード・生成日時)のサイドカーを付す。

## よく使うコマンド

```bash
# ACLED 認証(.bashrc は非対話シェルで早期 return するため直接読む)
eval "$(grep '^export ACLED_' ~/.bashrc)"

# データ取得・抽出(全てキャッシュ再現)
julia --project=experiments experiments/data/fetch_data.jl        # WB + ACLED
julia --project=experiments experiments/data/extract_vdem_p.jl    # p(V-Dem)
julia --project=experiments experiments/data/extract_swiid_g.jl   # g(SWIID)

# M8 ヒンドキャスト(--smoke = 較正窓のみ、--calib = 凍結済み較正値を使用)
julia --project=experiments -t 8 experiments/M8_hindcast.jl JPN THA --smoke --calib
julia --project=experiments -t 8 experiments/M8_hindcast.jl JPN THA --calib   # 検証ラン

# EKI 較正(結果は experiments/output/M8_calib_*.json → 凍結は M8_frozen_config.toml)
julia --project=experiments -t 8 experiments/M8_calibrate.jl THA JPN --J 24 --iters 4 --N 100
```

## 技術的な罠(解決済み。再導入しないこと)

- **σ_s の log 座標剛性**: 載荷 L<0 で drift の L·e^{−ξ} が発散 → EM に tame+floor
  ガード導入済み(#0032、src/integrator.jl)。双子実験の領域では不発火。
- **ACLED カウントの過分散**: 素のポアソン重みでは ESS が 1 に崩壊 → 1/ν
  テンパリング(#0033)。ν は EKI で較正せずポアソン最尤 ν* = ΣN/ΣΛ で
  プロファイル化。
- **phaseb_agreement テスト**: stable レジームの分散比はジャンプ未経験メンバー
  条件付き(#0029)。重い裾 × FP 環境差で全体分散比は CI で不安定。
- **VariableRateJump はスレッド並列で発火しない**(#0021)。Phase B は逐次実行。
- タイの ACLED は深南部4県(Pattani/Narathiwat/Yala/Songkhla)の慢性反乱が
  死者の85% — 国政レベル分析では除外必須(#0026/#0030)。
- 外生入力(netgrowth, wbar)は §8.4 の双子実験定数でなく較正窓データから
  フィットする(experiments/M8_hindcast.jl の fit_exogenous)。

## データ来歴の要点

- ACLED: OAuth password grant(ENV: ACLED_USERNAME/PASSWORD)。Research access は
  イベント単位データに直近12ヶ月エンバーゴ。JPN カバレッジ 2018〜、THA 2010〜。
  認証フローの検証記録: https://github.com/aray-99/acled-client
- 観測系列は experiments/data/raw/ に CSV + .meta.json(来歴)で全てコミット済み。
