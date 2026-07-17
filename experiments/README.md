# experiments/

Mirai予報のデータ取得・観測整形・較正・walk-forward検証・将来予報を行う
実験スクリプト群。

- `data/`: データ取得・整形スクリプト(`fetch_data.jl`、`extract_vdem_p.jl`、
  `extract_swiid_g.jl`、`build_observations.jl`、`prepare_events.jl`)と
  国別設定(`data/countries/<ISO3>.toml`)。
- `M8_*` / `M9_*` / `M10_*` / `M11_*`: 較正(M8)・walk-forward検証
  (M9/M10)・将来予報(M11)の各実験スクリプト。
- `output/`: 各スクリプトの生成物(JSON + 来歴サイドカー)。git 管理外
  (`.gitignore` 参照)。

## forecast_pipeline.jl の使い方(Issue #5)

データ取得 → 観測整形 → 較正 → walk-forward 自己検証 → 将来予報を
単一コマンドで通すドライバ。各段は既存スクリプトをサブプロセスとして
呼び出すだけで、ロジックは複製しない。

```
julia --project=experiments -t 8 experiments/forecast_pipeline.jl ISO3 \
    [--from STAGE] [--to STAGE] [--smoke] [--recalibrate] [--force]
```

### 引数

| 引数 | 説明 |
| --- | --- |
| `ISO3` | 対象国コード(1個。`experiments/data/countries/<ISO3>.toml` が必要) |
| `--from` | 開始段(既定 `fetch`) |
| `--to` | 終了段(既定 `forecast`) |
| `--smoke` | walkforward/forecast をスモークモード(短時間)で実行 |
| `--recalibrate` | calibrate 段で凍結値を使わず再較正(M8_calibrate.jl)を実行 |
| `--force` | fetch 段でキャッシュ済み CSV を再取得 |

### 段

| 段名 | 呼び出すスクリプト | 成果物 | ETA |
| --- | --- | --- | --- |
| `fetch` | `data/fetch_data.jl` → `data/extract_vdem_p.jl` → `data/extract_swiid_g.jl` | `data/raw/<ISO3>_*.csv`(+ `.meta.json`) | (キャッシュ有無で変動) |
| `obs` | `data/build_observations.jl` → `data/prepare_events.jl` | 検証サマリのみ(観測はダウンストリームでメモリ内構築)。`<ISO3>_events.csv` | 数分未満 |
| `calibrate` | (既定)凍結値を再利用してスキップ。`--recalibrate` 時は `M8_calibrate.jl` | `output/M8_calib_<ISO3>.json` | 約 30〜60 分/国(再較正時) |
| `walkforward` | `M10_walkforward.jl ISO3 --mu-gbar-sd 0.3 [--smoke]` | `output/M10_walkforward_<ISO3>[_smoke].json` | smoke: 数分 / 通常: 約 45 分〜2 時間/国(環境依存) |
| `forecast` | `M11_forecast.jl ISO3 [--smoke]` | `output/M11_forecast_<ISO3>[_smoke].json`(+ `.meta.json`) | smoke: 約 1 分 / 通常: 約 7 分/国 |

### 実行例

全段通し(新規国、fetch から forecast まで):

```
julia --project=experiments -t 8 experiments/forecast_pipeline.jl KEN
```

キャッシュ済みの JPN について観測整形〜予報だけ実行:

```
julia --project=experiments -t 8 experiments/forecast_pipeline.jl JPN --from obs --to forecast
```

スモーク動作確認(walkforward + forecast を短時間で):

```
julia --project=experiments -t 8 experiments/forecast_pipeline.jl JPN --from walkforward --to forecast --smoke
```

### 前提と注意

- (i) fetch の ACLED 取得には `ENV["ACLED_USERNAME"]` / `ENV["ACLED_PASSWORD"]`
  が必要。キャッシュ済み CSV(`data/raw/`)があれば再取得はスキップされる。
- (ii) JPN/THA は `M8_frozen_config.toml` の凍結較正値を再利用するのが既定
  (calibrate 段)。
- (iii) 新規国は `data/countries/<ISO3>.toml` の作成に加え、現行の
  M8〜M11 スクリプトが参照する `COUNTRY_CFG`(`experiments/M8_hindcast.jl`)が
  JPN/THA 前提のため、横展開(M13)ではこの externalize が別途必要
  (Issue #1 棚卸しの既知事項)。
- (iv) 較正値の凍結(`M8_frozen_config.toml` への記録)はオーナーの
  DECISIONS 手続きが必要(#0050 流儀)。

予報 JSON のスキーマは `docs/FORECAST_JSON.md` を参照。
