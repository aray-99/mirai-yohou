# 次セッション引き継ぎ(2026-07-20 更新)

> 本ファイルはセッション引き継ぎ用の**状態スナップショット**(git 追跡、#0092 で
> 転換)。一次資料は docs/DECISIONS.md と GitHub Issues — 齟齬があればそちらが正。
> @claude ランナーへの指示ではない(タスクは Issue 本文で自己完結させる)。
> 更新はセッション内で随時、コミットは節目(ターン終了・長時間ラン起動前・
> バッチ完了)ごと。

## 進行中: Issue #33 バッチ 2(EGY)— **合格(#0091)**、M11 予報生成中(2026-07-20)

- **EKI 較正凍結 #0090**(ν*=17.6)→ **walk-forward 合格 #0091(是正 0 回)**:
  被覆 0.974(n=38)/RMSE 7/8/logL −406.8 > −613.7/テスト green。develop
  6f8ff5e push 済み
- **実行中: M11_forecast.jl EGY**(setsid、ログ = experiments/output/
  EGY_m11_20260720.log、マーカー = 同 .done、ETA ≈30 分)→ 成果物 =
  experiments/output/M11_forecast_EGY.json
- **完了後の残作業(#0086 手順 5)**: ① M11 JSON 検分・コミット ②
  M12_dashboard.jl:1385 の iso3_list に "EGY" 追加 + EGY 用色スロット(#0085
  色表に追記、dataviz でパレット再検証。5 カ国目 — 既存 = JPN 青/THA 青緑/
  KOR 紫/TUR 橙)③ テストフィクスチャ(KOR/TUR 時は CCC/DDD 追加が前例)
  ④ ダッシュボード再生成・コミット(オーナー、#0074 決定 5)⑤ フルテスト →
  バッチクローズ DECISIONS(#0085 相当)→ main 反映判断(vM* タグ)は
  ユーザー相談
- **#27 完了・クローズ済み**: @claude(Fable 5 + skill、#0089 是正後)が
  MODEL_GUIDE.md を起草 → PR #36 → オーナーレビュー(数値突合 14 項目、
  是正 2 点: パラメータ数 6+ν 訂正・SPEC §8.1 出典分離)→ フルテスト green →
  develop fc39be5 マージ。ワークフローは修正済みで今後の @claude 委譲は
  Fable 5 + tag モードで動く(claude[bot] 自己コメントの重複発火は if 条件で
  skip、無害)

## (旧)EKI 較正メモ(2026-07-20 起動)

- 完了済み: #0086(バッチ凍結)→ 手順 1 データ整備 + 手順 2 診断 #0087(develop
  0ec96cd)→ #35 tau 収載・クローズ(44dfac2、**EGY は Wave 7 欠測で 2001/2013 の
  2 点 — 続行判断は #0088 に記録**)→ #0088 実行プロトコル凍結(9b43152: calib_t
  [21,29]/verif [29,35]/オリジン 29:33 = 5 個/EKI J24/iters4/N100 シード 20260710/
  theta_sig 含む)→ 付随コード変更(9d140ea、テスト 2510 green)
- **実行中: EGY EKI 較正**(experiments/M8_calibrate.jl EGY --J 24 --iters 4 --N 100)。
  setsid 切り離し済み。ログ = experiments/output/EGY_calib_20260720.log、完了マーカー =
  同 .done(exit=0 を確認すること)。ETA 30-60 分
- **較正完了後の手順**: ① 結果 JSON(experiments/output/M8_calib_*.json)をオーナー
  検分 → M8_frozen_config.toml に [EGY] 凍結(#0050 流儀、8771269 が前例)→ コミット
  ② M10 walk-forward(5 オリジン ≈1.3-1.5h、setsid)→ ③ #0052 判定(基準不変更・
  是正上限 2 回)→ ④ 合格ならダッシュボード追加(M12_dashboard.jl iso3_list に EGY、
  M11 予報 JSON 生成)+ #0085 色表に EGY スロット追記 + パレット再検証
- develop は push 済み 44dfac2 まで + 未 push = 9b43152/9d140ea(較正緑後にまとめて push)

## 現在地(要点)

- **Phase 3(M11-M13)完了**: M13 合格クローズ = DECISIONS #0085、**vM13 タグ、main 4cba108**。
  公開ダッシュボード = 4 カ国(JPN/THA/KOR/TUR)稼働、Pages 配信バイト一致確認済み:
  <https://aray-99.github.io/mirai-yohou/dashboard/>
- **Phase 4 はユーザー指示で一旦ステイ**(#28 Hawkes 相互励起 / #29 高頻度観測(金融=センサー
  結論記録済み)/ #30 wbar 定数近似の限界 — 全てバックログ維持、着手しない)
- DECISIONS 最新 = #0085(M13 クローズ)。#0083 に国別色・国旗の追記、#0084 = #0070 改定
  (フラッグシップ実装基盤、ユーザー承認済み)

## フラッグシップ可視化(次の主戦線)

- **プロトタイプ v6 完成・ユーザー好評**: `~/tmp_proto/mirai_flagship_proto.html`
  (テンプレート = scratchpad/proto_template.html + world_slim.json + hist_obs.json +
  M11_forecast JSON ×4 を python で `__PLACEHOLDER__` 差し込み)。機能: 3D グローブ
  (自転方向・停止)/2D 地図/タイムライン再生/スパン 1-30 年/観測重畳(w・phi は
  %→分数整合済み。ACLED 生カウントはモデル観測と定義不一致のためヒーロー非重畳)/
  歴史イベント注釈/時間スケール群化/国別ランキング/拡大モーダル + 国比較
  (国別固定色 + 簡略 SVG 国旗)/デュアルテーマ(ライト = DADS 写像)
- **本実装 = Issue #32(エピック)**: #0084 の不変原則(実行時外部リクエストゼロ・
  決定的ビルド・Pages・義務要件 a/b・JA/EN)下で vendored + ビルド可。正典は段階的置換
  (M12 生成器はパリティ達成後に凍結保守へ)。**Pages リリース予定はユーザー確認済み**
- デザイン正典 = #0083(+追記): デュアルテーマトークン、国別色
  JPN 青 #7096f8/#264af4・THA 青緑 #2fbf9a/#197a3e・KOR 紫 #c792ea/#6f42c1・
  TUR 橙 #f2a33c/#b25000(dataviz 検証済み、ダークの明度帯逸脱は国旗+直接ラベルで担保)

## 進行中の委譲(要監視 — ユーザーの GitHub 通知で検知、自動監視ループ禁止)

1. **#27 モデル解説+V&V 文書(docs/MODEL_GUIDE.md)を @claude に起草委譲済み**
   (2026-07-19、issuecomment-5015075805)。完了通知が来たら: PR 起票(base=develop)→
   科学的正確性レビュー(数値は DECISIONS 出典 ID 付き引用の検証)→ マージ。
   Julia テスト不要(docs のみ)だがフルテストは流儀どおり実行
2. **#35 [cowork] WVS tau 転記(EGY + IDN/PHL/MYS 先行)** — ユーザーが Claude Cowork で
   実施する運用。成果はコメント貼付 → コミットはオーナー(来歴主義)
3. workflow prompt に常時遵守事項を注入済み(5b0a15e、main 反映 4cba108)—
   PR #25 の件名規約違反の再発防止

## 次の意思決定(ユーザー相談待ち)

- **#32(フラッグシップ実装)の着手タイミング**(#33 EGY と並行可能 — フロントエンド vs
  科学パイプラインで独立)。@claude 余剰枠は「実装専任」タスクに割当てる方針
  (ユーザー指示。意思決定・重い計算を挟むものは不可)
- #33 は着手済み(上記「進行中」参照)。EGY の較正着手(手順 3 以降)のみ #35 待ち

## バックログ Issue 一覧

- #26 レポート機能(スコープ別、フラッグシップ実装後)
- #27 V&V 文書(@claude 起草中)
- #28 Phase 4: Hawkes 相互励起(ステイ)
- #29 高頻度観測: 金融=センサー・LLM 指標(ステイ)
- #30 wbar/netgrowth 定数近似 → UN WPP 時変化(ステイ。V&V 文書に限界として記載)
- #32 フラッグシップ実装エピック
- #33 国拡大常設トラック(EGY → IDN/PHL/MYS。IND/BRA/MEX は枠組み不適合で見送り)
- #34 Plurality アイデア(ユーザー好感触。優先順位づけ時に再訪)
- #35 [cowork] WVS tau 転記

## 運用ルール(変更なし)

- セッション開始時は本ファイルから再開。状況が変わるたびに即更新
- Fable はオーケストレーター。探索・実装はサブエージェント/@claude 委譲
- 凍結基準・しきい値の変更はユーザー明示承認 + DECISIONS 必須。決定は事前記録(#0039)
- 長時間ラン(>30分)は setsid + 完了マーカー + Monitor、起動前 ETA 報告
- 友人 PR レビュー手順: worktree で checkout → diff 精査(凍結ファイル即差し戻し)→
  スタンドアロン + フルテスト自走(2510 本 ≈4-6 分)→ 決定性確認 → --no-ff マージ →
  公開物再生成はオーナー(#0074 決定 5)→ main 反映(Pages バイト一致確認)
- コミット規約: 件名英語・(#XXXX) 禁止・Refs 行。gh 書き込みは gh api が確実
- プロトタイプ共有: SendUserFile が使えない文脈があるため ~/tmp_proto/ 経由も併用
