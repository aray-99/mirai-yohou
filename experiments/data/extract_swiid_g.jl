# g(ジニ)の抽出(DECISIONS #0023/#0030): SWIID の gini_disp(可処分所得
# ベース)を Harvard Dataverse から取得し、raw/<ISO3>_g_swiid.csv
# (year,value。value = gini_disp / 100 ∈ [0,1])を生成する。
#
# - ソース: doi:10.7910/DVN/LM4OWF(SWIID 9.92、datafile id は下記定数)。
# - SWIID はモデルベース補間系列のため観測ノイズ R は膨らませて使う(#0023)。
#   各年の標準誤差(gini_disp_se / 100)を meta.json に保存し、M8 の R 設計で
#   R = (se × r_inflation)^2 として利用する。
# - World Bank Gini(<ISO3>_g.csv、生・疎)は感度分析用に併存。
#
# 実行: julia --project=experiments experiments/data/extract_swiid_g.jl [--force]

using Downloads
using JSON3
using Dates

const RAW_DIR = joinpath(@__DIR__, "raw")
const SCRATCH = joinpath(@__DIR__, "scratch")          # .gitignore 対象
const SWIID_DOI = "doi:10.7910/DVN/LM4OWF"
const SWIID_VERSION = "9.92"
const SWIID_DATAFILE_ID = 13657070                     # swiid9_92.zip
const SWIID_URL = "https://dataverse.harvard.edu/api/access/datafile/$(SWIID_DATAFILE_ID)"
const COUNTRIES = [("JPN", "Japan"), ("THA", "Thailand")]

function swiid_summary_path()
    mkpath(SCRATCH)
    csv = joinpath(SCRATCH, "swiid$(replace(SWIID_VERSION, "." => "_"))_summary.csv")
    if !isfile(csv)
        zip = joinpath(SCRATCH, "swiid.zip")
        if !isfile(zip)
            println("downloading SWIID $(SWIID_VERSION) (~25MB)...")
            Downloads.download(SWIID_URL, zip)
        end
        run(`unzip -o -q -j $zip -d $SCRATCH`)
        isfile(csv) || error("展開後に $(basename(csv)) が見つかりません")
    end
    return csv
end

function main()
    force = "--force" in ARGS
    targets = [joinpath(RAW_DIR, "$(iso)_g_swiid.csv") for (iso, _) in COUNTRIES]
    if all(isfile, targets) && !force
        foreach(p -> println("cached: $p"), targets)
        return
    end
    src = swiid_summary_path()
    mkpath(RAW_DIR)
    for (iso, cname) in COUNTRIES
        rows = Tuple{Int, Float64, Float64}[]   # (year, gini_disp, se)
        for (i, line) in enumerate(eachline(src))
            i == 1 && continue
            f = split(line, ",")
            f[1] == cname || continue
            (isempty(f[3]) || isempty(f[4])) && continue
            push!(rows, (parse(Int, f[2]), parse(Float64, f[3]), parse(Float64, f[4])))
        end
        sort!(rows; by = first)
        csvpath = joinpath(RAW_DIR, "$(iso)_g_swiid.csv")
        open(csvpath, "w") do io
            println(io, "year,value")
            for (y, g, _) in rows
                println(io, y, ",", g / 100)
            end
        end
        write(csvpath * ".meta.json", JSON3.write(Dict(
            "source" => "SWIID $(SWIID_VERSION) gini_disp (Solt, Harvard Dataverse)",
            "doi" => SWIID_DOI, "datafile_id" => SWIID_DATAFILE_ID,
            "url" => SWIID_URL, "country" => cname,
            "normalization" => "gini_disp / 100 -> 0-1",
            "note" => "SWIID はモデルベース補間系列。R は se を基に膨らませて使用(#0023/#0030)",
            "measurement_se" => Dict(string(y) => se / 100 for (y, _, se) in rows),
            "fetched_at" => Dates.format(now(), dateformat"yyyy-mm-ddTHH:MM:SS"))))
        println("extracted: $csvpath ($(length(rows)) 点  $(rows[1][1])–$(rows[end][1]))")
    end
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
