# p(分極度)の抽出(DECISIONS #0027): V-Dem の政治分極指数 v2cacamps を
# 公式 R パッケージ vdemdata(GitHub)の vdem.RData から抽出し、
# raw/<ISO3>_p.csv(year,value。value = v2cacamps_osp / 4 ∈ [0,1])を生成する。
#
# - ソース: https://github.com/vdeminstitute/vdemdata (data/vdem.RData)。
#   来歴(パッケージ版・データファイルのコミット SHA)はサイドカー JSON に記録。
# - vdem.RData(約34MB)は Git 管理外のスクラッチにキャッシュし、コミットしない。
# - 測定不確実性(v2cacamps_osp_sd / 4)は meta.json に系列で保存(M8 の R 設計用)。
#
# 実行: julia --project=experiments experiments/data/extract_vdem_p.jl [--force]

using Downloads
using JSON3
using Dates
using RData
using DataFrames

const RAW_DIR = joinpath(@__DIR__, "raw")
const SCRATCH = joinpath(@__DIR__, "scratch")          # .gitignore 対象
const VDEM_URL = "https://github.com/vdeminstitute/vdemdata/raw/master/data/vdem.RData"
const COUNTRIES = ["JPN", "THA"]
const YEAR_FROM = 1990

function vdem_dataframe()
    mkpath(SCRATCH)
    rdata = joinpath(SCRATCH, "vdem.RData")
    if !isfile(rdata)
        println("downloading vdem.RData (~34MB)...")
        Downloads.download(VDEM_URL, rdata)
    end
    return load(rdata)["vdem"]
end

"vdemdata リポジトリの来歴(データファイルの最新コミット SHA)を取得。失敗時は空文字。"
function vdem_provenance()
    try
        buf = IOBuffer()
        Downloads.download("https://api.github.com/repos/vdeminstitute/vdemdata" *
                           "/commits?path=data/vdem.RData&per_page=1", buf)
        c = JSON3.read(String(take!(buf)))[1]
        return string(c.sha), string(c.commit.committer.date)
    catch
        return "", ""
    end
end

function main()
    force = "--force" in ARGS
    targets = [joinpath(RAW_DIR, "$(c)_p.csv") for c in COUNTRIES]
    if all(isfile, targets) && !force
        foreach(p -> println("cached: $p"), targets)
        return
    end
    df = vdem_dataframe()
    sha, sha_date = vdem_provenance()
    mkpath(RAW_DIR)
    for iso in COUNTRIES
        sub = dropmissing(df[(df.country_text_id .=== iso) .& (df.year .>= YEAR_FROM),
                             [:year, :v2cacamps_osp, :v2cacamps_osp_sd]])
        sort!(sub, :year)
        csvpath = joinpath(RAW_DIR, "$(iso)_p.csv")
        open(csvpath, "w") do io
            println(io, "year,value")
            for r in eachrow(sub)
                println(io, Int(r.year), ",", r.v2cacamps_osp / 4)
            end
        end
        write(csvpath * ".meta.json", JSON3.write(Dict(
            "source" => "V-Dem v2cacamps (political polarization) via vdemdata R package",
            "url" => VDEM_URL, "country" => iso,
            "normalization" => "v2cacamps_osp / 4 (original scale 0-4 -> 0-1)",
            "vdem_data_commit" => sha, "vdem_data_commit_date" => sha_date,
            "measurement_sd" => Dict(string(Int(r.year)) => r.v2cacamps_osp_sd / 4
                                     for r in eachrow(sub)),
            "fetched_at" => Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))))
        println("extracted: $csvpath ($(nrow(sub)) 点  $(Int(sub.year[1]))–$(Int(sub.year[end])))")
    end
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
