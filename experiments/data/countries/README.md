# 国別データ設定(Issue #3)

`<ISO3>.toml` を1本置くだけで、データスクリプト群(`fetch_data.jl` /
`extract_vdem_p.jl` / `extract_swiid_g.jl` / `prepare_events.jl` /
`build_observations.jl`)の対象国に追加される(引数省略時は本ディレクトリの
全 ISO3 が対象)。

- `name_en`: ACLED / World Bank 検索用の英語国名。
- `swiid_name_en`(省略可): SWIID の国名表記が `name_en` と異なる場合の上書き
  (例: KOR は ACLED/WB が "South Korea"、SWIID が "Korea"。#0077)。
- `regime`: `build_params` に渡すレジーム(`:stable` / `:volatile` 等)。
- `[acled] year_from`: ACLED 取得リクエストの下限年。実際に返るデータ範囲が
  これより後ろから始まる場合がある(例: JPN は 2010 をリクエストしても
  ACLED 側に 2018 年より前のデータが無い)ため、実カバレッジは
  `acled_from` を参照する。
- `[acled] exclude_admin1`: 国政レベル分析で除外する admin1 のリスト
  (慢性的な地域紛争など。空配列も可)。
- `[acled] acled_from`: ACLED の実カバレッジ開始日(TOML 日付リテラル)。
  `experiments/data/raw/<ISO3>_events.csv.meta.json` の `date_range` 開始日と
  一致させる(#0022)。M8 ヒンドキャストの観測窓カバレッジ判定に使う。
- `[acled] calib_from` / `calib_to`: ジャンプカタログの較正期間(TOML 日付リテラル)。
- `[hindcast] calib_t` / `verif_t`: M8 ヒンドキャストの較正/検証ウィンドウ
  (1990-01-01 起点の経過年、Float の2要素配列)。セクション自体が省略可能な
  国もあり、その場合はウィンドウ未確定として `M8_hindcast.jl` のロードが
  明確なエラーで停止する(#0077。KOR/TUR は Issue #20 でオーナーが確定するまで
  未設定)。

`regime` と `exclude_admin1` は国固有の科学的判断であり(#0026/#0030 での
検討に相当)、新規国追加時はオーナーによる fable-required なレビューが必要。
`[hindcast]` のウィンドウ設定も同様に科学的判断であり、オーナーレビューが必要。
