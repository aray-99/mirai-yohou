# M11 パイプラインドライバ(Issue #5): データ取得 → 観測整形 → 較正 →
# walk-forward 自己検証 → 将来予報を単一コマンドで通す。
#
# 各段は既存スクリプトをサブプロセスとして呼び出すだけで、ロジックは複製しない
# (MiraiYohou を using しない。TOML/JSON3/Dates のみ)。
#
# 実行: julia --project=experiments -t 8 experiments/forecast_pipeline.jl ISO3
#       [--from STAGE] [--to STAGE] [--smoke] [--recalibrate] [--force]
# STAGE = fetch / obs / calibrate / walkforward / forecast(この順序)。

using TOML
using JSON3
using Dates

include(joinpath(@__DIR__, "data", "country_config.jl"))

const STAGES = ["fetch", "obs", "calibrate", "walkforward", "forecast"]
const REPO_DIR = dirname(@__DIR__)
const EXP_DIR = @__DIR__
const RAW_DIR = joinpath(EXP_DIR, "data", "raw")
const OUT_DIR = joinpath(EXP_DIR, "output")

"experiments/ 配下のスクリプトをサブプロセスとして実行する(--project は
experiments ディレクトリ絶対パス、スレッド数は親の Threads.nthreads() を引き継ぐ)。
失敗時は run が例外を送出するので、段名を添えて rethrow する。"
function julia_script(script, args...)
    julia_bin = Base.julia_cmd().exec[1]
    script_path = joinpath(EXP_DIR, script)
    # 引数は明示的な argv 構築で渡す(cmd リテラル補間の splat は引数を落とす
    # ことがあるため使わない)。実行コマンドをログに残し、引数の脱落を可視化する。
    argv = String[julia_bin, "--project=$(EXP_DIR)", "-t", string(Threads.nthreads()),
                  script_path, String.(collect(args))...]
    println("  実行: ", join([basename(script_path); String.(collect(args))], " "))
    run(Cmd(Cmd(argv); dir = REPO_DIR))
end

"path が存在すれば成果物として表示し、来歴サイドカー(path.meta.json)や
JSON 本体の provenance キーがあれば主要項目を1行で表示する。"
function print_provenance(path)
    if !isfile(path)
        return
    end
    println("  成果物: $path")
    keys_of_interest = ["source", "fetched_at", "placed_at", "commit", "seed", "generated_at"]
    meta_path = path * ".meta.json"
    if isfile(meta_path)
        meta = JSON3.read(read(meta_path, String))
        parts = String[]
        for k in keys_of_interest
            if haskey(meta, Symbol(k)) || haskey(meta, k)
                v = haskey(meta, Symbol(k)) ? meta[Symbol(k)] : meta[k]
                push!(parts, "$k=$v")
            end
        end
        if !isempty(parts)
            println("    来歴: ", join(parts, " "))
        end
    else
        # M10 walkforward JSON は .meta.json サイドカーが無く、本体に
        # provenance キーを持つ(commit / seed / generated_at)。
        try
            doc = JSON3.read(read(path, String))
            if haskey(doc, :provenance)
                prov = doc[:provenance]
                parts = String[]
                for k in ["commit", "seed", "generated_at"]
                    if haskey(prov, Symbol(k))
                        push!(parts, "$k=$(prov[Symbol(k)])")
                    end
                end
                if !isempty(parts)
                    println("    来歴: ", join(parts, " "))
                end
            end
        catch
            # 本体が JSON でない、または provenance が無い場合は何もしない
        end
    end
end

"段の開始を告知する(eta が nothing でなければ所要時間の目安を表示)"
function stage_banner(name, eta)
    println("\n===== 段 $name =====")
    if eta !== nothing
        println("  予想所要時間(ETA): $eta")
    end
end

# ---- 各段の実装 ----

function run_fetch(iso3::String; force::Bool = false)
    cfg_path = joinpath(EXP_DIR, "data", "countries", "$(iso3).toml")
    isfile(cfg_path) || error("国別設定がありません。countries/README.md 参照")
    extra = force ? ("--force",) : ()
    julia_script(joinpath("data", "fetch_data.jl"), iso3, extra...)
    julia_script(joinpath("data", "extract_vdem_p.jl"), iso3, extra...)
    julia_script(joinpath("data", "extract_swiid_g.jl"), iso3, extra...)
    if isdir(RAW_DIR)
        for f in sort(readdir(RAW_DIR))
            if startswith(f, "$(iso3)_") && endswith(f, ".csv")
                print_provenance(joinpath(RAW_DIR, f))
            end
        end
    end
end

function ensure_fetched(iso3::String)
    isdir(RAW_DIR) || error("fetch 段が未実行です(--from fetch から実行してください)")
    any(f -> startswith(f, "$(iso3)_") && endswith(f, ".csv"), readdir(RAW_DIR)) ||
        error("fetch 段が未実行です(--from fetch から実行してください)")
end

function run_obs(iso3::String)
    ensure_fetched(iso3)
    julia_script(joinpath("data", "build_observations.jl"), iso3)
    julia_script(joinpath("data", "prepare_events.jl"), iso3)
    println("  観測はダウンストリームがメモリ内で構築(build_observations.jl)。本段は整形の検証サマリのみ")
    print_provenance(joinpath(RAW_DIR, "$(iso3)_events.csv"))
end

function frozen_section(iso3::String)
    frozen_path = joinpath(EXP_DIR, "M8_frozen_config.toml")
    isfile(frozen_path) || return nothing, frozen_path
    cfg = TOML.parsefile(frozen_path)
    return get(cfg, iso3, nothing), frozen_path
end

function run_calibrate(iso3::String; recalibrate::Bool = false)
    ensure_fetched(iso3)
    section, frozen_path = frozen_section(iso3)
    if section !== nothing && !recalibrate
        println("  凍結較正値を再利用(--recalibrate で再較正)")
        println("  凍結値ファイル: $frozen_path")
        for (k, v) in sort(collect(section); by = first)
            println("    $k = $v")
        end
        return
    end
    if !recalibrate
        error("凍結較正値がありません。--recalibrate で較正を実行してください")
    end
    println("  予想所要時間(ETA): 約 30〜60 分/国(M8 EKI 既定設定)")
    julia_script("M8_calibrate.jl", iso3)
    print_provenance(joinpath(OUT_DIR, "M8_calib_$(iso3).json"))
    println("  注意: 較正値の凍結(M8_frozen_config.toml への記録)はオーナーの DECISIONS 手続きが必要(#0050 流儀)。walkforward/forecast 段は凍結値を参照する")
end

function ensure_frozen(iso3::String)
    section, _ = frozen_section(iso3)
    section === nothing && error("凍結較正値がありません(calibrate 段 + オーナー凍結が必要)")
end

function run_walkforward(iso3::String; smoke::Bool = false)
    ensure_frozen(iso3)
    args = smoke ? (iso3, "--mu-gbar-sd", "0.3", "--smoke") : (iso3, "--mu-gbar-sd", "0.3")
    julia_script("M10_walkforward.jl", args...)
    suffix = smoke ? "_smoke" : ""
    print_provenance(joinpath(OUT_DIR, "M10_walkforward_$(iso3)$(suffix).json"))
end

function run_forecast(iso3::String; smoke::Bool = false)
    ensure_frozen(iso3)
    args = smoke ? (iso3, "--smoke") : (iso3,)
    julia_script("M11_forecast.jl", args...)
    suffix = smoke ? "_smoke" : ""
    print_provenance(joinpath(OUT_DIR, "M11_forecast_$(iso3)$(suffix).json"))
end

# ---- ETA(walkforward/forecast は run_* 内で表示するため、ここでは
# stage_banner に渡す値を段ごとに用意する) ----

function eta_for(stage::String, smoke::Bool)
    if stage == "fetch" || stage == "obs"
        return nothing
    elseif stage == "calibrate"
        return nothing  # calibrate は凍結値スキップか --recalibrate かで分岐するため run_calibrate 内で表示
    elseif stage == "walkforward"
        return smoke ? "数分" : "約 45 分〜2 時間/国(環境依存。M10 実測)"
    elseif stage == "forecast"
        return smoke ? "約 1 分" : "約 7 分/国(M11 実測)"
    end
    return nothing
end

function usage()
    println("""
    使い方: julia --project=experiments -t 8 experiments/forecast_pipeline.jl ISO3 [--from STAGE] [--to STAGE] [--smoke] [--recalibrate] [--force]
    STAGE = fetch / obs / calibrate / walkforward / forecast(この順序)
    """)
end

function main()
    args = ARGS
    smoke = "--smoke" in args
    recalibrate = "--recalibrate" in args
    force = "--force" in args
    from_stage = "fetch"
    to_stage = "forecast"
    for (i, a) in enumerate(args)
        if a == "--from"
            from_stage = args[i + 1]
        elseif a == "--to"
            to_stage = args[i + 1]
        end
    end

    flag_names = Set(["--from", "--to", "--smoke", "--recalibrate", "--force"])
    positional = String[]
    skip_next = false
    for (i, a) in enumerate(args)
        if skip_next
            skip_next = false
            continue
        end
        if a in ("--from", "--to")
            skip_next = true
            continue
        end
        if a in flag_names
            continue
        end
        push!(positional, a)
    end

    if length(positional) != 1
        usage()
        error("ISO3 を1個指定してください(現在 $(length(positional)) 個)")
    end
    iso3 = positional[1]

    if !(from_stage in STAGES)
        error("不明な段: $from_stage(有効: fetch/obs/calibrate/walkforward/forecast)")
    end
    if !(to_stage in STAGES)
        error("不明な段: $to_stage(有効: fetch/obs/calibrate/walkforward/forecast)")
    end
    from_idx = findfirst(==(from_stage), STAGES)
    to_idx = findfirst(==(to_stage), STAGES)
    if from_idx > to_idx
        error("--from ($from_stage) が --to ($to_stage) より後です")
    end

    timings = Dict{String, Float64}()
    total_t0 = time()

    for idx in from_idx:to_idx
        stage = STAGES[idx]
        eta = eta_for(stage, smoke)
        stage_banner(stage, eta)
        t0 = time()
        try
            if stage == "fetch"
                run_fetch(iso3; force = force)
            elseif stage == "obs"
                run_obs(iso3)
            elseif stage == "calibrate"
                run_calibrate(iso3; recalibrate = recalibrate)
            elseif stage == "walkforward"
                run_walkforward(iso3; smoke = smoke)
            elseif stage == "forecast"
                run_forecast(iso3; smoke = smoke)
            end
        catch e
            println(stderr, "段 $stage で失敗しました")
            rethrow(e)
        end
        timings[stage] = time() - t0
    end

    total = time() - total_t0
    breakdown = join(["$s $(round(timings[s], digits = 1))s" for s in STAGES[from_idx:to_idx]], " / ")
    println("\n== パイプライン完了: $iso3 $from_stage→$to_stage 合計 $(round(total, digits = 1)) 秒(段別: $breakdown)==")
end

if abspath(PROGRAM_FILE) == (@__FILE__)
    main()
end
