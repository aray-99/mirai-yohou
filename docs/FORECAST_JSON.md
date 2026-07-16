# M11 将来予報 JSON 契約(M12 実装者向け)

`experiments/M11_forecast.jl` が生成する `experiments/output/M11_forecast_{ISO3}.json`
(+ `.meta.json` 来歴サイドカー)のスキーマを定義する。M12 ダッシュボードは
本ファイルの構造を前提に実装してよい。output/ は git 管理外(.gitignore、
§0.5.7)なので、常にスクリプト再実行で再生成するのが正である。

## 生成コマンド

```
julia --project=experiments -t 8 experiments/M11_forecast.jl JPN THA
julia --project=experiments -t 8 experiments/M11_forecast.jl JPN THA --horizon 30 --N 100
julia --project=experiments -t 8 experiments/M11_forecast.jl THA --smoke   # 動作確認(N=40, 5年)
```

引数: 国コード(省略時 `JPN THA` の両方)、`--horizon`(既定 30.0 年)、
`--N`(既定 100)、`--smoke`(動作確認モード。N=40・EKI J=12/iters=2/N_eki=40、
`--horizon` 未指定時は 5 年に短縮、出力ファイル名に `_smoke` サフィックス)。
シードは `20260711`(#0063)固定で変更不可。

## ファイルパス規約

- `experiments/output/M11_forecast_{ISO3}.json`(スモークは `_{ISO3}_smoke.json`)
- 同名 `.meta.json`: 来歴(`provenance`)のみを複製したサイドカー。JSON 本体
  を読まずに来歴だけ確認したい用途向け(`M11_horizon_smoke.jl` と同じ流儀)。
- `output/` は `.gitignore` 対象(大容量アンサンブル生成物、§0.5.7)。
  リポジトリに存在しない場合は上記コマンドで再生成すること。

## トップレベルキー一覧

| キー | 型 | 意味 |
|---|---|---|
| `country` | string | 国 ISO3 コード(`"JPN"` / `"THA"`) |
| `regime` | string | モデルレジーム(`"stable"` / `"volatile"`) |
| `N` | int | アンサンブルメンバー数 |
| `horizon_years` | int | 予報ホライズン(年、整数) |
| `forecast_start` | object | `{t, calendar_year}`。予報起点(データ末尾 t_k) |
| `verified_horizon` | object | `{years, until_t, until_calendar_year, note}`。検証済み範囲(6年、#0052/#0069)。これ以遠は外挿領域 — M12 は視覚的に区別する義務がある |
| `finite_check` | object | `{pass, nonfinite_values}`。予報アンサンブル全状態の NaN/Inf 検査 |
| `assimilation` | object | 同化・EKI 較正の要約(下記) |
| `jump_thinning` | object | Γ ジャンプ超過確率シンニング p_ex の要約(下記) |
| `variables` | object | 9系列 + De の年次分位点(下記) |
| `count_forecast` | object | 年次期待騒乱イベント数の年次分位点(下記) |
| `years` | array[int] | 年次グリッドの西暦年(`t` と対応、長さ `horizon_years + 1`) |
| `t` | array[float] | 年次グリッドの保持時刻(`forecast_start.t + 0, 1, ..., horizon_years`) |
| `endogenous_jumps_total` | int | 全メンバー・全予報期間の内生ジャンプ発生数合計 |
| `provenance` | object | 来歴(下記) |
| `elapsed_sec` | float | 当該国の実行所要時間(秒) |

### `assimilation`

| キー | 意味 |
|---|---|
| `theta_hat` | EKI 最終オリジン再較正後の較正パラメータ中心値(`{param名 => 値}`) |
| `nu_star` | カウント尤度スケール ν*(NegBin/ポアソン共通) |
| `r_hat` | NegBin 分散パラメータ r̂(カウント窓なしなら `null` = 実質ポアソン) |
| `include_theta_sig` | theta_sig の L3 拡大適用可否(#0052/#0054 データ規則) |
| `sigma_n_total_counts` | theta_sig 判定に使った同化窓内フィルタ後週次カウント合計 ΣN |
| `mu_gbar_prior` | `{center, sd, rule}`。mu_gbar アンカリング prior(#0062、sd=0.3 #0063) |
| `nresample` | 同化ラン中の系統リサンプリング回数 |
| `ess_range` | `[min, max]`。有効サンプルサイズの範囲 |
| `n_obs` | 較正・同化に使った観測レコード数 |
| `window` | `[win_start, t_k]`。較正・同化窓 |

### `jump_thinning`

| キー | 意味 |
|---|---|
| `p_ex` | Γ ジャンプ採択確率のシンニング係数(#0068 式、将来予報モードは #0071 規則) |
| `n_forced_weeks` | 較正・同化窓内の強制ジャンプ週数 |
| `n_count_weeks` | 較正・同化窓内のカウントデータ週数(ACLED カバレッジ内) |
| `p_ex_fallback` | `"no_count_data"`(0/0 フォールバック発動、#0071) または `null` |
| `design_decision` | 適用規則の参照(`"#0068 式 / #0071 将来予報フォールバック規則"`) |

**p_ex_fallback の意味**(#0071): 同化窓内にカウントデータ週が1つも無い国
(例: ACLED カバレッジ開始が同化窓より後)は p_ex の分母がゼロになるため、
将来予報モードでは `p_ex = 0.0`(予報 Γ ジャンプを一切発生させない)に
フォールバックする。#0068(M10 walk-forward)のフォールバック値 1.0 とは
**異なる**値である点に注意(1年先予報では実害が小さいという #0068 の前提が
30年ホライズンでは成立しないため、#0071 で将来予報モード専用に凍結された
規則)。`p_ex_fallback: "no_count_data"` が立っている場合、当該国・オリジン
はジャンプリスク評価が実質不能であることを意味する。M12 ダッシュボードは
この状態を(例えば count_forecast パネルに警告バッジを出すなどして)通常の
p_ex > 0 のケースと区別して表示すること。

### `variables`

キーは `P, w, y, g_swiid, T_proxy, phi, v, tau, p, De` の10系列(9実系列 + 1診断量)。
各エントリの構造:

| キー | 型 | 意味 |
|---|---|---|
| `unit` | string | 自然単位(下表) |
| `transform` | string | 保持座標 → 自然単位の変換の一言説明 |
| `last_observation_year` | int または null | 当該系列の最終観測年(西暦)。De は null(状態変数でないため) |
| `q05`, `q25`, `q50`, `q75`, `q95` | array[float] | 年次分位点(長さ `horizon_years + 1`、`years`/`t` と同じ添字。メンバー横断) |

#### 座標変換の規約(src/coordinates.jl の規約、保持座標→自然単位)

モデル内部は保持座標(log または logit)で状態を持つ。各系列は状態インデックス
`xi` から次のとおり自然単位へ変換する(`base_*` は raw CSV の1990年値、
無ければ最初の観測年値。`experiments/data/build_observations.jl` の
`_load_series`/`_baseline` と同一基準):

| 変数 | 状態変数 | 保持座標 | 自然単位への変換 | 単位 |
|---|---|---|---|---|
| P | xi_P | log | `base_P * exp(xi_P)` | persons |
| w | xi_w | logit | `sigmoid(xi_w)` | fraction (0-1) |
| y | (合成観測) | — | `base_y * exp(h_logy(xi))` | USD (2015, per capita) |
| g_swiid | xi_g | logit | `sigmoid(xi_g)` | fraction (0-1) |
| T_proxy | xi_T | log | `base_T * exp(xi_T)` | patent applications (resident) |
| phi | xi_phi | logit | `sigmoid(xi_phi)` | fraction (0-1) |
| v | xi_v | log | `base_v * exp(xi_v)` | subscriptions per 100 people |
| tau | xi_tau | logit | `sigmoid(xi_tau)` | fraction (0-1) |
| p | xi_p | logit | `sigmoid(xi_p)` | fraction (0-1) |
| De | (診断量) | — | `deborah_number(l2, sigmoid(xi_g), sigmoid(xi_p))` | dimensionless |

**y の合成観測**: 一人当たり産出 y は状態に直接含まれず、TFP・資本・人的
資本・労働参加率から合成する観測演算子 `h_logy`(SPEC §9.1、
`build_observations.jl` の `_h_logy`)で毎時刻代数的に計算する:

```
A(T, phi) = A0 * exp(theta_T * xi_T + theta_phi * sigmoid(xi_phi))
y = A * exp(alpha * xi_k + (1 - alpha) * (xi_h + log(sigmoid(xi_w))))
```

(自然単位版。実装は log 空間で合成し `base_y * exp(...)` で戻す。)

**De の定義式**(SPEC §7、社会的デボラ数): `l2` は共有パラメータ(較正後
`params.l2`。De の構成パラメータ eta_g/g_c/eta_p/delta_sig は L3 状態拡大
(#0046)の対象外の固定値):

```
De = ( eta_g * max(gbar - g_c, 0) + eta_p * pbar ) / delta_sig
```

De < 1 は応力が忘却で緩和される安定レジーム、De > 1 は単調蓄積する変動
レジームを意味する診断量(状態でも観測でもない)。

### `count_forecast`

騒乱リスクの年次予報(**期待イベント数**の分位点であり、実現カウントの
分位点ではない点に注意)。

| キー | 型 | 意味 |
|---|---|---|
| `unit` | string | `"expected political disorder events per year, nu* × ∫intensity"` |
| `nu_star` | float | カウント尤度スケール(`assimilation.nu_star` と同値) |
| `r_hat` | float または null | NegBin 分散パラメータ(`assimilation.r_hat` と同値) |
| `note` | string | 「NegBin 分散は r_hat 参照」の注記 |
| `q05`〜`q95` | array[float] | 予報年 y=1..horizon_years の年次期待イベント数分位点(長さ `horizon_years`。y=0 の起点断面は含まない) |

各メンバー・各予報年について `nu_star * ∫ intensity(state, params) dt`
(年内 dt=0.01 刻みの積分、intensity は Hawkes 到着レート `min(lam_bar, lam_b + lam_e)`)
をメンバー横断で分位点化したものである。個々のメンバー値そのものが
「その年に起きるイベント数」ではなく、その年の期待到着数(ポアソン/NegBin
の平均パラメータ)である — 分散(過分散)は `r_hat` 側で表現される。

## 年次グリッドの定義

`years[i]` / `t[i]`(i = 1..horizon_years+1)は予報起点からの経過年
y = 0, 1, ..., horizon_years に対応する。**y=0 は予報起点の断面**(初期条件
= 同化窓末尾 t_k のアンサンブル事後分布そのもの。予報積分は未実施)であり、
y=1 以降が実際の予報値である。`variables` 内の各分位配列は `years`/`t` と
同じ添字を共有する(重複を避けるため配列自体は `variables` 側に持たず、
トップレベルの `years`/`t` を参照する設計)。`count_forecast` は y=0 を
含まないため長さが1つ短い(y=1..horizon_years の `horizon_years` 個)。

## `verified_horizon` の意味

検証済みホライズンは**起点 + 6 年**(Issue #4 完了条件、PHASE3_DESIGN §2/§3。
実証的裏付けは M9/M10 の expanding-window walk-forward — #0052 設計・#0069
合格判定 — による)。`verified_horizon.until_t` /
`until_calendar_year` はこの検証済み範囲の終端を表す。**これより先
(y > 6)は統計的検証を経ていない外挿領域**であり、M12 ダッシュボードは
検証済み区間と外挿領域を視覚的に区別して表示する義務がある(#0070)。
検証済み範囲を「予報」、外挿領域を無条件に同列表示することは、モデルの
実証的裏付けの範囲について利用者に誤解を与えるため避けること。

## 免責と誠実さの注記

`count_forecast` および `variables` の分位帯は、外挿領域(y > 6)では
モデルの内生力学(Hawkes 自己励起・OU 型緩和・L3 拡大パラメータの継続 RW)
のみに基づく統計的投影であり、観測による制約を受けていない。特に p_ex の
0/0 フォールバック(`p_ex_fallback: "no_count_data"`)が発動している国・
オリジンでは、ジャンプ関連の力学(Γ 適用による急激な社会的応力上昇)が
一切表現されないため、`count_forecast` は「ジャンプが起きない前提」の
過小評価である可能性がある。この限界は `jump_thinning.p_ex_fallback` を
参照して利用者に明示すること。
