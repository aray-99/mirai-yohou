# Phase 2 データ基盤(M7、DECISIONS #0023/#0024): 観測系列の取得・整形・来歴管理
#
# - World Bank API から自動取得(P, w, y, g, T代理, phi, v)。年次。
# - tau(WVS)・p(V-Dem)はローカルファイルのローダ(手動ダウンロード。
#   各 load_* の docstring 参照)。
# - ACLED はアクセスキーが必要(ENV["ACLED_ACCESS_KEY"] / ["ACLED_EMAIL"])。
#   キー未設定ではスタブがガイダンスを出して停止する。
# - 出力: experiments/data/raw/<ISO3>_<var>.csv(year,value)+ 来歴サイドカー
#   JSON(ソース URL・指標 ID・取得日時・API の lastupdated)。キャッシュ済み
#   ファイルは再取得しない(force=true で上書き)。
#
# 実行: julia --project=experiments experiments/data/fetch_data.jl [--force]

using Downloads
using JSON3
using Dates

const RAW_DIR = joinpath(@__DIR__, "raw")
const COUNTRIES = ["JPN", "THA"]                      # 日本・タイ(#0023)

# 変数 → World Bank 指標 ID(#0023。SWIID は手動ローダ、WB Gini は自動側)
const WB_INDICATORS = [
    (:P, "SP.POP.TOTL"),          # 総人口
    (:w, "SP.POP.1564.TO.ZS"),    # 生産年齢人口比率(%)
    (:y, "NY.GDP.PCAP.KD"),       # 実質 GDP per capita(2015 USD、年次。四半期化は M8)
    (:g, "SI.POV.GINI"),          # ジニ係数(疎。SWIID はローカルローダ側)
    (:T_proxy, "IP.PAT.RESD"),    # 特許出願(居住者)
    (:phi, "IT.NET.USER.ZS"),     # インターネット利用率(%)
    (:v, "IT.CEL.SETS.P2"),       # モバイル契約 /100人
]

function fetch_wb(country::String, indicator::String; force::Bool = false)
    mkpath(RAW_DIR)
    var = first(name for (name, id) in WB_INDICATORS if id == indicator)
    csvpath = joinpath(RAW_DIR, "$(country)_$(var).csv")
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
tau(制度信頼)ローダ: World Values Survey の信頼設問(国・波別集計)を
手動で experiments/data/raw/<ISO3>_tau.csv(year,value。value は 0〜1)に
置く。WVS は https://www.worldvaluessurvey.org からの登録ダウンロード。
"""
function load_tau(country::String)
    p = joinpath(RAW_DIR, "$(country)_tau.csv")
    isfile(p) || error("$(p) がありません。WVS から国別の制度信頼系列を作成して置いてください(値は 0〜1)")
    return load_series(p)
end

"""
p(分極度)ローダ: V-Dem の政治分極指数(v2cacamps 等を 0〜1 に正規化)を
手動で experiments/data/raw/<ISO3>_p.csv に置く。
V-Dem は https://v-dem.net からデータセット(CSV)をダウンロード。
"""
function load_p(country::String)
    p = joinpath(RAW_DIR, "$(country)_p.csv")
    isfile(p) || error("$(p) がありません。V-Dem の分極指数系列を作成して置いてください(値は 0〜1)")
    return load_series(p)
end

"""
ACLED イベント取得スタブ: ENV["ACLED_ACCESS_KEY"] と ENV["ACLED_EMAIL"] が必要
(https://acleddata.com で無償登録)。キー設定後に実装を有効化する(M8)。
"""
function fetch_acled(country::String; force::Bool = false)
    haskey(ENV, "ACLED_ACCESS_KEY") && haskey(ENV, "ACLED_EMAIL") ||
        error("ACLED のアクセスキーが未設定です。https://acleddata.com で登録し、" *
              "ENV[\"ACLED_ACCESS_KEY\"] と ENV[\"ACLED_EMAIL\"] を設定してください")
    error("ACLED フェッチャは M8 で実装(キー確認後)")
end

function main()
    force = "--force" in ARGS
    for c in COUNTRIES, (var, id) in WB_INDICATORS
        try
            fetch_wb(c, id; force)
        catch e
            println("WARN: $c/$id failed: ", sprint(showerror, e))
        end
    end
    println("\n-- 取得サマリ --")
    for c in COUNTRIES, (var, _) in WB_INDICATORS
        p = joinpath(RAW_DIR, "$(c)_$(var).csv")
        if isfile(p)
            ys, vs = load_series(p)
            println(rpad("$(c)_$(var)", 14), length(ys), " 点  ",
                    isempty(ys) ? "" : "$(ys[1])–$(ys[end])")
        end
    end
    println("\ntau(WVS)・p(V-Dem)・ACLED は手動/キー設定が必要(各 docstring 参照)")
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
