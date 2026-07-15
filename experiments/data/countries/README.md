# 国別データ設定(Issue #3)

`<ISO3>.toml` を1本置くだけで、データスクリプト群(`fetch_data.jl` /
`extract_vdem_p.jl` / `extract_swiid_g.jl` / `prepare_events.jl` /
`build_observations.jl`)の対象国に追加される(引数省略時は本ディレクトリの
全 ISO3 が対象)。

- `name_en`: ACLED / SWIID 検索用の英語国名。
- `regime`: `build_params` に渡すレジーム(`:stable` / `:volatile` 等)。
- `[acled] year_from`: ACLED 取得開始年。
- `[acled] exclude_admin1`: 国政レベル分析で除外する admin1 のリスト
  (慢性的な地域紛争など。空配列も可)。
- `[acled] calib_from` / `calib_to`: ジャンプカタログの較正期間(TOML 日付リテラル)。

`regime` と `exclude_admin1` は国固有の科学的判断であり(#0026/#0030 での
検討に相当)、新規国追加時はオーナーによる fable-required なレビューが必要。
