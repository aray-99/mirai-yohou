# ACLED イベントの前処理ユーティリティ(M8。設計は PHASE2_DESIGN §9.1)
#
# raw/<ISO3>_events.csv(fetch_data.jl の ACLED スキーマ)から
#   (1) 週次カウント系列(ポアソン重み用)
#   (2) 強制ジャンプ候補カタログ(週次規模の分位閾値)
# を生成する。フィルタ・閾値はすべて引数であり、値の凍結は M8 の
# DECISIONS エントリで行う(§9.1 は未承認のドラフト)。
#
# 実行(サマリ表示): julia --project=experiments experiments/data/prepare_events.jl [ISO3...]

using Dates

include(joinpath(@__DIR__, "country_config.jl"))

const EVENTS_RAW_DIR = joinpath(@__DIR__, "raw")

struct AcledEvent
    id::String
    date::Date
    disorder_type::String
    event_type::String
    sub_event_type::String
    admin1::String
    fatalities::Int
end

function load_events(country::String)
    path = joinpath(EVENTS_RAW_DIR, "$(country)_events.csv")
    isfile(path) || error("$(path) がありません。fetch_data.jl を先に実行してください")
    events = AcledEvent[]
    for (i, line) in enumerate(eachline(path))
        i == 1 && continue
        f = split(line, ",")
        push!(events, AcledEvent(f[1], Date(f[2]), f[4], f[5], f[6], f[7],
                                 parse(Int, f[8])))
    end
    return events
end

"""
PHASE2_DESIGN §9.1 のドラフト既定: 深南部4県(タイの慢性反乱)。
一次ソースは countries/THA.toml(#0026/#0030)。M8_hindcast.jl の
COUNTRY_CFG がこの定数を参照するため、TOML 由来の値のまま定数として残置する。
"""
const DEEP_SOUTH_THA = Vector{String}(load_country_config("THA")["acled"]["exclude_admin1"])

"""
政治的騒乱イベントの抽出(§9.1 ドラフト既定)。
disorder_type が Political violence(併記含む)のもの + Riots。
exclude_admin1 の県は除外(タイの深南部フィルタ)。
"""
function political_events(events::Vector{AcledEvent};
                          exclude_admin1::Vector{String} = String[])
    return filter(events) do e
        e.admin1 in exclude_admin1 && return false
        occursin("Political violence", e.disorder_type) || e.event_type == "Riots"
    end
end

"""
週次カウント系列: 開始日 t0(週境界)から終了日まで、7日窓のイベント件数と
死者数を返す。返り値 (week_start::Vector{Date}, counts, fatalities)。
"""
function weekly_counts(events::Vector{AcledEvent};
                       t0::Date = firstdayofweek(minimum(e.date for e in events)),
                       t1::Date = maximum(e.date for e in events))
    nweeks = Int(cld(Dates.value(t1 - t0) + 1, 7))
    counts = zeros(Int, nweeks)
    fatal = zeros(Int, nweeks)
    for e in events
        e.date < t0 && continue
        w = Int(fld(Dates.value(e.date - t0), 7)) + 1
        w > nweeks && continue
        counts[w] += 1
        fatal[w] += e.fatalities
    end
    return [t0 + Day(7 * (w - 1)) for w in 1:nweeks], counts, fatal
end

"""
強制ジャンプ候補カタログ: 週次規模 s_w = fatalities_w + count_scale * counts_w が
較正期間 [calib_from, calib_to] 内の分位 q を超える週を返す。
閾値は較正期間のみから計算し、返り値には全期間の該当週を含む
(検証期間への凍結適用、§9.2)。min_size は閾値の下限フロア
(イベントが疎な国で分位が退化し単発週を拾うのを防ぐ。安定国では
カタログが空になるのが正常)。
返り値: (weeks::Vector{Date}, magnitudes::Vector{Float64}, threshold)
"""
function jump_catalog(week_start::Vector{Date}, counts::Vector{Int},
                      fatal::Vector{Int};
                      calib_from::Date, calib_to::Date,
                      q::Float64 = 0.98, count_scale::Float64 = 0.1,
                      min_size::Float64 = 5.0)
    s = fatal .+ count_scale .* counts
    incal = [calib_from <= w <= calib_to for w in week_start]
    any(incal) || error("較正期間内に週がありません")
    scal = sort(s[incal])
    thr = max(scal[clamp(ceil(Int, q * length(scal)), 1, length(scal))], min_size)
    sel = s .>= thr
    return week_start[sel], s[sel], thr
end

function summarize(country::String; exclude_admin1 = String[],
                   calib_from::Date, calib_to::Date)
    ev = political_events(load_events(country); exclude_admin1)
    isempty(ev) && (println("$country: 政治的騒乱イベントなし"); return)
    ws, c, f = weekly_counts(ev)
    println("$country: $(length(ev)) 件 ($(ws[1])〜)、週数 $(length(ws))、",
            "週平均 $(round(sum(c)/length(c), digits=2)) 件")
    jw, jm, thr = jump_catalog(ws, c, f; calib_from, calib_to)
    println("  ジャンプ候補週(q=0.98、較正 $(calib_from)〜$(calib_to)、閾値 $(round(thr, digits=1))):")
    for (w, m) in zip(jw, jm)
        println("    ", w, "  規模 ", round(m, digits = 1))
    end
end

function main()
    for c in country_args(ARGS)
        cfg = load_country_config(c)
        summarize(c; exclude_admin1 = Vector{String}(cfg["acled"]["exclude_admin1"]),
                 calib_from = cfg["acled"]["calib_from"], calib_to = cfg["acled"]["calib_to"])
    end
end

abspath(PROGRAM_FILE) == (@__FILE__) && main()
