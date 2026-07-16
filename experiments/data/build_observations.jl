# 実データ観測ビルダ(M8、DECISIONS #0030): raw/ の系列を保持座標の
# ObservationRecord 列に変換する。
#
# - 時間座標: t = (year − 1990) + 0.5(年央近似、§14-6。t=0 が 1990-01-01)。
# - スケール正規化(#0030-2): P, y, T代理, v は 1990 年値(欠測なら最初の
#   観測年の値)で除して log。w, g, phi, tau, p は logit。
# - 観測ノイズ sd は保持座標(log/logit)単位。点別の測定不確実性が
#   ある系列(g_swiid の se、p の V-Dem osp_sd)はデルタ法で logit 座標へ
#   変換して使う。数値は暫定であり、検証ラン前に DECISIONS で凍結する
#   (#0030-7)。
#
# 実行(サマリ): julia --project=experiments experiments/data/build_observations.jl [ISO3...]

using JSON3
using MiraiYohou
using MiraiYohou: ObservationSpec, ObservationRecord,
                  IX_P, IX_W, IX_K, IX_G, IX_T, IX_PHI, IX_V, IX_TAU, IX_PP,
                  logit

include(joinpath(@__DIR__, "country_config.jl"))

const OBS_RAW_DIR = joinpath(@__DIR__, "raw")
const BASE_YEAR = 1990.0                      # #0030-2

"CSV(year,value)ローダ(fetch_data.jl と同形式)"
function _load_series(csvpath::String)
    years = Int[]; values = Float64[]
    for (i, line) in enumerate(eachline(csvpath))
        i == 1 && continue
        y, v = split(line, ",")
        push!(years, parse(Int, y)); push!(values, parse(Float64, v))
    end
    return years, values
end

_series_path(country, file) = joinpath(OBS_RAW_DIR, "$(country)_$(file).csv")

"meta.json のサイドカーから点別不確実性(year→sd、自然座標)を読む"
function _load_point_sd(csvpath::String, key::String)
    meta = JSON3.read(read(csvpath * ".meta.json", String))
    haskey(meta, key) || return Dict{Int, Float64}()
    return Dict(parse(Int, string(k)) => Float64(v) for (k, v) in pairs(meta[key]))
end

"1990 年値(なければ最初の観測値)で正規化する基準値"
function _baseline(years, values)
    i = findfirst(==(Int(BASE_YEAR)), years)
    return values[something(i, 1)]
end

"logit 座標の sd へのデルタ法変換: sd_logit = sd_nat / (x(1−x))"
_logit_sd(sd_nat, x) = sd_nat / max(x * (1 - x), 1e-3)

"""
実データ観測系列の定義(#0030-1)。各要素は
(変数名, ファイル名, 変換, 状態インデックスまたは :log_y, 暫定 sd, 点別sdキー)。
ObservationSpec 名は var のまま(Symbol(s.var))、ファイルパスは file を使う
(P は raw/<ISO3>_pop.csv、p は raw/<ISO3>_pol.csv。Issue #3 の大小文字衝突回避)。
sd は保持座標単位の暫定値(検証前に DECISIONS で凍結、#0030-7)。
"""
const REAL_SERIES = [
    (var = "P",       file = "pop",     kind = :log_norm,  target = IX_P,    sd = 0.005, sdkey = ""),
    (var = "w",       file = "w",       kind = :logit_pct, target = IX_W,    sd = 0.01,  sdkey = ""),
    (var = "y",       file = "y",       kind = :log_y,     target = :log_y,  sd = 0.02,  sdkey = ""),
    (var = "g_swiid", file = "g_swiid", kind = :logit,     target = IX_G,    sd = 0.05,  sdkey = "measurement_se", sdinflate = 2.0),
    (var = "T_proxy", file = "T_proxy", kind = :log_norm,  target = IX_T,    sd = 0.15,  sdkey = ""),
    (var = "phi",     file = "phi",     kind = :logit_pct, target = IX_PHI,  sd = 0.05,  sdkey = ""),
    (var = "v",       file = "v",       kind = :log_norm,  target = IX_V,    sd = 0.05,  sdkey = ""),
    (var = "tau",     file = "tau",     kind = :logit,     target = IX_TAU,  sd = 0.15,  sdkey = ""),
    (var = "p",       file = "pol",     kind = :logit,     target = IX_PP,   sd = 0.10,  sdkey = "measurement_sd", sdinflate = 1.0),
]

"log y の観測演算子(標準観測 §9.1 と同一の合成)"
function _h_logy(params)
    l1, l2 = params.l1, params.l2
    return xi -> begin
        phi = 1 / (1 + exp(-xi[IX_PHI]))
        logA = log(l1.A0) + l2.theta_T * xi[IX_T] + l2.theta_phi * phi
        logA + l1.alpha * xi[IX_K] +
            (1 - l1.alpha) * (xi[IX_H] + log(1 / (1 + exp(-xi[IX_W]))))
    end
end

"""
    build_observations(country, params; t1) -> Vector{ObservationRecord}

raw/ の実データから保持座標の観測レコード列(時刻順)を構築する。
`t1` より後(および t < 0、つまり 1990 年より前)の観測は捨てる。
tau など任意系列はファイルがなければスキップ(警告表示)。
"""
function build_observations(country::String, params; t1::Float64)
    records = ObservationRecord[]
    for s in REAL_SERIES
        path = _series_path(country, s.file)
        if !isfile(path)
            println("WARN: $(basename(path)) なし — $(s.var) をスキップ")
            continue
        end
        years, values = _load_series(path)
        point_sd = isempty(s.sdkey) ? Dict{Int, Float64}() :
                   _load_point_sd(path, s.sdkey)
        base = s.kind in (:log_norm, :log_y) ? _baseline(years, values) : 1.0
        for (y, v) in zip(years, values)
            t = (y - BASE_YEAR) + 0.5
            (0 <= t <= t1) || continue
            if s.kind === :log_norm
                z = log(v / base)
            elseif s.kind === :logit_pct
                z = logit(clamp(v / 100, 1e-4, 1 - 1e-4))
            elseif s.kind === :logit
                z = logit(clamp(v, 1e-4, 1 - 1e-4))
            else  # :log_y
                z = log(v / base)
            end
            sd = s.sd
            if haskey(point_sd, y)
                # 点別不確実性をデルタ法で logit 座標へ。sdinflate は補間系列の
                # R 膨張(SWIID: 2.0、#0030-1)。V-Dem は測定モデルの事後 sd
                # そのままなので 1.0。
                x = s.kind === :logit_pct ? v / 100 : v
                sd = max(sd, get(s, :sdinflate, 1.0) * _logit_sd(point_sd[y], x))
            end
            h = s.target === :log_y ? _h_logy(params) :
                (ix -> (xi -> xi[ix]))(s.target)
            # target_ix(DECISIONS #0043): 恒等写像の実データ観測のみ状態行を記録
            # (:log_y のような合成観測は 0 のまま = スプレッド床の対象外)。
            target_ix = s.target === :log_y ? 0 : s.target
            spec = ObservationSpec(Symbol(s.var), 1.0, sd, h, target_ix)
            push!(records, ObservationRecord(t, spec, z))
        end
    end
    sort!(records; by = r -> r.t)
    return records
end

function _summarize(country)
    params = build_params(Symbol(load_country_config(country)["regime"]))  # regime は国別 TOML 由来(Issue #3)
    recs = build_observations(country, params; t1 = 35.0)
    println("== $country: $(length(recs)) 観測 ==")
    for s in REAL_SERIES
        rs = [r for r in recs if r.spec.name === Symbol(s.var)]
        isempty(rs) && continue
        vals = [r.value for r in rs]
        sds = unique(round.([r.spec.sd for r in rs], digits = 3))
        println(rpad(s.var, 9), length(rs), " 点  t=[",
                round(rs[1].t, digits = 1), ",", round(rs[end].t, digits = 1),
                "]  z=[", round(minimum(vals), digits = 2), ",",
                round(maximum(vals), digits = 2), "]  sd=",
                length(sds) <= 3 ? sds : "$(minimum(sds))..$(maximum(sds))")
        @assert all(abs.(vals) .< 10) "$(s.var): |z| >= 10(§9.4 警告水準)"
    end
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    for c in country_args(ARGS)
        _summarize(c)
    end
end
