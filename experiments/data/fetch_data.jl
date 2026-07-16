# Phase 2 データ基盤(M7、DECISIONS #0023/#0024): 観測系列の取得・整形・来歴管理
#
# - World Bank API から自動取得(P, w, y, g, T代理, phi, v)。年次。
# - tau(WVS)・p(V-Dem)はローカルファイルのローダ(手動ダウンロード。
#   各 load_* の docstring 参照)。
# - ACLED は OAuth(password grant)で取得(ENV["ACLED_USERNAME"] /
#   ["ACLED_PASSWORD"])。認証仕様と権限(Research access)の経緯は
#   DECISIONS #0026 と https://github.com/aray-99/acled-client 参照。
# - 出力: experiments/data/raw/<ISO3>_<file>.csv(year,value)+ 来歴サイドカー
#   JSON(ソース URL・指標 ID・取得日時・API の lastupdated)。キャッシュ済み
#   ファイルは再取得しない(force=true で上書き)。
#
# 実行: julia --project=experiments experiments/data/fetch_data.jl [ISO3...] [--force]

using Downloads
using JSON3
using Dates

include(joinpath(@__DIR__, "country_config.jl"))

const RAW_DIR = joinpath(@__DIR__, "raw")

# 変数 → World Bank 指標 ID・出力ファイル名(#0023。SWIID は手動ローダ、WB Gini は
# 自動側)。:P だけ file="pop"(<ISO3>_pol.csv[分極度]との大小文字衝突回避、Issue #3)。
const WB_INDICATORS = [
    (:P, "SP.POP.TOTL", "pop"),           # 総人口
    (:w, "SP.POP.1564.TO.ZS", "w"),       # 生産年齢人口比率(%)
    (:y, "NY.GDP.PCAP.KD", "y"),          # 実質 GDP per capita(2015 USD、年次。四半期化は M8)
    (:g, "SI.POV.GINI", "g"),             # ジニ係数(疎。SWIID はローカルローダ側)
    (:T_proxy, "IP.PAT.RESD", "T_proxy"), # 特許出願(居住者)
    (:phi, "IT.NET.USER.ZS", "phi"),      # インターネット利用率(%)
    (:v, "IT.CEL.SETS.P2", "v"),          # モバイル契約 /100人
]

function fetch_wb(country::String, indicator::String; force::Bool = false)
    mkpath(RAW_DIR)
    file = first(file for (_, id, file) in WB_INDICATORS if id == indicator)
    csvpath = joinpath(RAW_DIR, "$(country)_$(file).csv")
    if isfile(csvpath) && !force
        println("cached: $csvpath")
        return csvpath
    end
    url = "https://api.worldbank.org/v2/country/$country/indicator/$indicator" *
          "?format=json&per_page=200&date=1960:2030"
    buf = IOBuffer()
    Downloads.download(url, buf)
    doc = JSON3.read(String(take!(buf)))
    meta, rows = doc[1], doc[2]
    open(csvpath, "w") do io
        println(io, "year,value")
        for r in sort(collect(rows); by = r -> r.date)
            r.value === nothing && continue
            println(io, r.date, ",", r.value)
        end
    end
    write(csvpath * ".meta.json", JSON3.write(Dict(
        "source" => "World Bank API v2", "url" => url,
        "indicator" => indicator, "country" => country,
        "api_lastupdated" => get(meta, "lastupdated", ""),
        "fetched_at" => Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))))
    println("fetched: $csvpath")
    return csvpath
end

"CSV(year,value)を (years, values) に読む共通ローダ"
function load_series(csvpath::String)
    years = Int[]; values = Float64[]
    for (i, line) in enumerate(eachline(csvpath))
        i == 1 && continue
        y, v = split(line, ",")
        push!(years, parse(Int, y)); push!(values, parse(Float64, v))
    end
    return years, values
end

"""
tau(制度信頼)ローダ: World Values Survey の政府信頼設問を国・波別に集計し、
手動で experiments/data/raw/<ISO3>_tau.csv(year,value。value は 0〜1)に置く。

作成手順(ユーザー側作業):
- 設問: E069_11 "Confidence: The Government"。
  value = (A great deal + Quite a lot の回答割合) / 有効回答(0〜1)。
  year = 各波のフィールドワーク年(日本: 1981/1990/1995/2000/2005/2010/2019、
  タイ: 2007/2013/2018)。
- 経路A(登録不要): https://www.worldvaluessurvey.org → Data and Documentation
  → Online Analysis で国・波ごとに設問の度数表を表示し、割合を転記。
- 経路B(登録): 同サイトから Time-series (1981-2022) データを登録ダウンロードし、
  E069_11 を国・波で集計。
"""
function load_tau(country::String)
    p = joinpath(RAW_DIR, "$(country)_tau.csv")
    isfile(p) || error("$(p) がありません。WVS から国別の制度信頼系列を作成して置いてください(値は 0〜1)")
    return load_series(p)
end

"""
p(分極度)ローダ: V-Dem の政治分極指数(v2cacamps 等を 0〜1 に正規化)を
手動で experiments/data/raw/<ISO3>_pol.csv に置く。
V-Dem は https://v-dem.net からデータセット(CSV)をダウンロード。
"""
function load_p(country::String)
    p = joinpath(RAW_DIR, "$(country)_pol.csv")
    isfile(p) || error("$(p) がありません。V-Dem の分極指数系列を作成して置いてください(値は 0〜1)")
    return load_series(p)
end

# ---- ACLED(イベントデータ。DECISIONS #0026)----
#
# OAuth password grant(access_token 24時間)→ /api/acled/read を Bearer で
# ページング取得。Research access ティアは直近12ヶ月のイベント単位データが
# エンバーゴ(#0026)。認証フローの検証記録:
# https://github.com/aray-99/acled-client (TROUBLESHOOTING.md)

const ACLED_TOKEN_URL = "https://acleddata.com/oauth/token"
const ACLED_READ_URL = "https://acleddata.com/api/acled/read"
const ACLED_FIELDS = ["event_id_cnty", "event_date", "year", "disorder_type",
                      "event_type", "sub_event_type", "admin1", "fatalities"]

urlencode(s::AbstractString) =
    join(c in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-._~" ?
         string(c) : join("%" * uppercase(string(b, base = 16, pad = 2))
                          for b in codeunits(string(c)))
         for c in s)

formbody(pairs) =
    join(("$(urlencode(string(k)))=$(urlencode(string(v)))" for (k, v) in pairs), "&")

function acled_token()
    user = get(ENV, "ACLED_USERNAME", "")
    pass = get(ENV, "ACLED_PASSWORD", "")
    (isempty(user) || isempty(pass)) &&
        error("ACLED の認証情報が未設定です。myACLED アカウント(Research access 承認済み)の " *
              "ENV[\"ACLED_USERNAME\"] と ENV[\"ACLED_PASSWORD\"] を設定してください")
    body = formbody(["username" => user, "password" => pass,
                     "grant_type" => "password", "client_id" => "acled",
                     "scope" => "authenticated"])
    out = IOBuffer()
    resp = Downloads.request(ACLED_TOKEN_URL; method = "POST",
        input = IOBuffer(body), output = out,
        headers = ["Content-Type" => "application/x-www-form-urlencoded"])
    resp.status == 200 || error("ACLED トークン取得に失敗 (HTTP $(resp.status)): $(String(take!(out)))")
    return JSON3.read(String(take!(out))).access_token
end

"""
ACLED イベント取得: 国別のイベント単位データ(event_date, event_type, fatalities)を
取得し `<ISO3>_events.csv` にキャッシュする。年範囲は ACLED のカバレッジ開始〜現在
(エンバーゴ分はサーバ側で欠ける)。403 の場合はアカウント権限の問題
(https://github.com/aray-99/acled-client/blob/main/TROUBLESHOOTING.md)。
"""
function fetch_acled(country::String; force::Bool = false,
                     year_from::Union{Int, Nothing} = nothing)
    mkpath(RAW_DIR)
    csvpath = joinpath(RAW_DIR, "$(country)_events.csv")
    if isfile(csvpath) && !force
        println("cached: $csvpath")
        return csvpath
    end
    cfg = load_country_config(country)
    year_from = something(year_from, cfg["acled"]["year_from"])
    token = acled_token()
    cname = cfg["name_en"]
    rows = Any[]
    page, limit = 1, 5000
    while true
        query = join(["country=$(urlencode(cname))",
                      "year=$(year_from)|$(Dates.year(now()))", "year_where=BETWEEN",
                      "fields=$(join(ACLED_FIELDS, '|'))",
                      "limit=$limit", "page=$page", "_format=json"], "&")
        out = IOBuffer()
        resp = Downloads.request("$ACLED_READ_URL?$query"; method = "GET", output = out,
            headers = ["Authorization" => "Bearer $token"])
        resp.status == 200 || error("ACLED read に失敗 (HTTP $(resp.status)): " *
                                    "$(String(take!(out)))(403 なら権限未付与。#0026 参照)")
        doc = JSON3.read(String(take!(out)))
        batch = get(doc, "data", [])
        append!(rows, batch)
        length(batch) < limit && break
        page += 1
    end
    sort!(rows; by = r -> string(r.event_date))
    open(csvpath, "w") do io
        println(io, join(ACLED_FIELDS, ","))
        for r in rows   # 選択フィールドはいずれもカンマを含まない語彙(#0026)
            println(io, join((string(get(r, Symbol(f), "")) for f in ACLED_FIELDS), ","))
        end
    end
    write(csvpath * ".meta.json", JSON3.write(Dict(
        "source" => "ACLED API (OAuth, Research access)", "url" => ACLED_READ_URL,
        "country" => cname, "year_from" => year_from, "n_events" => length(rows),
        "date_range" => isempty(rows) ? "" :
            "$(rows[1].event_date)–$(rows[end].event_date)",
        "embargo_note" => "Research access は直近12ヶ月のイベント単位データを含まない(#0026)",
        "fetched_at" => Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))))
    println("fetched: $csvpath ($(length(rows)) events)")
    return csvpath
end

function main()
    force = "--force" in ARGS
    countries = country_args(ARGS)
    for c in countries, (var, id, file) in WB_INDICATORS
        try
            fetch_wb(c, id; force)
        catch e
            println("WARN: $c/$id failed: ", sprint(showerror, e))
        end
    end
    println("\n-- 取得サマリ --")
    for c in countries, (var, _, file) in WB_INDICATORS
        p = joinpath(RAW_DIR, "$(c)_$(file).csv")
        if isfile(p)
            ys, vs = load_series(p)
            println(rpad("$(c)_$(file)", 14), length(ys), " 点  ",
                    isempty(ys) ? "" : "$(ys[1])–$(ys[end])")
        end
    end
    if haskey(ENV, "ACLED_USERNAME") && haskey(ENV, "ACLED_PASSWORD")
        for c in countries
            try
                fetch_acled(c; force)
            catch e
                println("WARN: ACLED $c failed: ", sprint(showerror, e))
            end
        end
    else
        println("\nACLED はスキップ(ENV[\"ACLED_USERNAME\"]/[\"ACLED_PASSWORD\"] 未設定。#0026)")
    end
    println("\ntau(WVS)・p(V-Dem)は手動配置が必要(各 docstring 参照)")
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
