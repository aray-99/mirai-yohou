# 【ユーザー作業】WVS から tau(制度信頼)系列を作成する手順

M8 の観測系列のうち唯一残っているユーザー側作業です。
**完了したらこのファイルは削除して構いません**(手順の恒久版は
`experiments/data/fetch_data.jl` の `load_tau` docstring にあります)。

## 作るもの

以下の2ファイル(CSV、ヘッダ行つき、値は 0〜1):

```
experiments/data/raw/JPN_tau.csv
experiments/data/raw/THA_tau.csv
```

形式(1行目はヘッダ固定):

```csv
year,value
1981,0.28
1990,0.31
...
```

## 値の定義

- 設問: **E069_11 "Confidence: The Government"**(政府への信頼)
- 値 = **(A great deal + Quite a lot の回答数)÷ 有効回答数**
  (Don't know / No answer 等の無効回答は分母に含めない)
- year = その波のフィールドワーク実施年

## 対象の波

| 国 | 年(フィールドワーク) |
|----|------|
| 日本 | 1981, 1990, 1995, 2000, 2005, 2010, 2019 |
| タイ | 2007, 2013, 2018 |

(取得できない波があれば飛ばして構いません。行が欠けても
ローダはそのまま動きます)

## 取得経路(どちらでも可)

### 経路A: Online Analysis(登録不要、推奨)

1. <https://www.worldvaluessurvey.org> → **Data and Documentation**
   → **Online Analysis**
2. 波(Wave)と国を選択
3. 設問リストから **E069_11**(Confidence: The Government)を選び度数表を表示
4. 「A great deal」と「Quite a lot」の割合(有効回答ベース)を合算して転記

### 経路B: データダウンロード(登録あり)

1. 同サイトで登録し **WVS Time-series (1981–2022)** をダウンロード
2. 変数 E069_11 を国・波でクロス集計して上記定義で計算

## 完了後

ファイルを置いたら Claude に「tau を置いた」と伝えてください。
ローダ(`load_tau`)での読み込み確認と M8 への組み込みはこちらで行います。
