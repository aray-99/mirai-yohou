# 逐次・非同期同化ドライバ(SPEC §9.2/§9.3、双子実験 §13 用)
#
# - 観測は届いた時刻に、届いた観測だけで解析(グリッドに丸めた tstops)。
# - イベントは ExogenousEvents 相当: カタログ時刻で全メンバー強制ジャンプ
#   (マーク rho は各メンバー独立)。内生発火は同化中オフ(§9.3)。
# - 無イベント尤度: 週次窓(0.02年 ≈ 7.3日で近似、DECISIONS #0009)で
#   Λ_i = ∫ lam_i dt をトラッキングし、ポアソン重み+ESS<N/2 で系統再抽選。
# - 乗法的インフレーション 1.02 常時、強制ジャンプ直後の解析は 1.05(§9.3)。
# - E1b: theta_sig を log 座標で状態拡大(行 14。d(param)=0+微小ノイズ、§8.3)。
# - 駆動パラメータの L3 状態拡大の汎用化(DECISIONS #0046): augmented_params
#   (AugmentedParam のリスト)で任意個数・任意名のパラメータを状態拡大できる。
#   augmented=true(theta_sig 1個)は内部でこの記述子1個に変換される特例経路。

"""
同化の設定(既定値は SPEC §9/§13)。

`inflation_mode`(DECISIONS #0012/#0013):
- `:per_analysis` — 解析毎に偏差 × rho_inf(SPEC v1.0 の原記述。多レート観測
  では弱観測部分空間が複利膨張しフィルタ崩壊する — #0012)
- `:per_time` — 解析時に rho_inf^(Δt_前回解析から / tau_ref)。単位時間あたり
  注入率を解析頻度から切り離す
- `:rtps` — relaxation to prior spread(rtps_alpha)。観測に拘束されない
  成分への注入がゼロで、複利膨張が原理的に起きない
強制ジャンプ直後の「一時的に強める」(§9.3)は rho_inf_jump / rtps_alpha_jump。
"""
Base.@kwdef struct AssimConfig
    t0::Float64 = 0.0
    t1::Float64 = 45.0
    dt::Float64 = 0.01
    inflation_mode::Symbol = :rtps   # 採用方式(#0013)。:per_analysis は SPEC v1.0 原案
    rho_inf::Float64 = 1.02          # 乗法モードの基礎レート(§9.2)
    rho_inf_jump::Float64 = 1.05     # 強制ジャンプ直後(§9.3)
    tau_ref::Float64 = 0.25          # :per_time の正規化時定数(四半期)
    rtps_alpha::Float64 = 0.7        # :rtps の緩和係数(#0013 の診断マトリクスで選定)
    rtps_alpha_jump::Float64 = 0.8   # 強制ジャンプ直後の :rtps 係数
    event_window::Float64 = 0.02     # 週次バッチ近似(#0009)
    ess_ratio::Float64 = 0.5         # ESS < N * ratio で再抽選(§9.3)
    param_noise_sd::Float64 = 0.01   # 状態拡大パラメータの微小ノイズ(/√年)
    smoother_lag::Float64 = 0.0      # 固定ラグ EnKS のラグ(年)。0 = 平滑化オフ(#0024)
    smoother_dt::Float64 = 0.1       # 平滑化スナップショットの間隔(年)
    # 平滑化更新は period ≥ この値の観測を含む解析のみで行う。高頻度観測は
    # 過去状態への実情報が乏しく、クロス共分散の標本雑音だけが累積するため(#0025)
    smoother_min_period::Float64 = 2.0
    # 平滑化で更新する状態行(変数局所化、#0025)。実情報のない変数
    # (k 等)への雑音蓄積を防ぐ。既定は制度ブロック + 格差。
    smoother_vars::Vector{Int} = [IX_G, IX_TAU, IX_TAUA, IX_SIG, IX_PP]
    # tauA(IX_TAUA)への緩い擬似観測の倍率(DECISIONS #0036)。tau 観測と
    # 同時刻・同値の擬似観測を sd = このスカラー × tau 観測 sd で追加する。
    # 既定 0.0 = オフ(従来動作。E1・既存テストの記録結果を保護)。
    tauA_pseudo_sd_mult::Float64 = 0.0
    # 現在時刻解析の変数局所化(DECISIONS #0040-(α))。ここに列挙した状態行は、
    # 解析バッチが analysis_unmask_names のいずれの観測名も含まない場合、
    # 現在時刻の EnKF 更新(K の該当行)をマスクする。EnKS の smoother_vars/
    # smooth_rows(過去平滑化の局所化、#0025)と対になる現在時刻側の局所化。
    # 既定は空 = 従来動作(後方互換)。
    analysis_masked_vars::Vector{Int} = Int[]
    # analysis_masked_vars のマスクを解除する観測名(#0040-(α))。解析バッチに
    # この name の観測が1つでも含まれればマスクを解除する。
    analysis_unmask_names::Vector{Symbol} = Symbol[]
    # 観測座標の加法的スプレッド下限(DECISIONS #0043)。解析直前の事前
    # アンサンブルで、観測座標(target_ix ≠ 0 の恒等写像観測に限る)の sd が
    # floor = このスカラー × 観測 sd を下回る場合、対応する状態行へ独立
    # ガウス摂動を加えて sd を floor まで回復してから解析する。既定 0.0 =
    # オフ(後方互換)。ゲイン消失(#0042-1)への対処。
    obs_spread_floor_frac::Float64 = 0.0
    # 再抽選後の Liu-West 若返り縮小係数(DECISIONS #0057)。系統再抽選が
    # 発火した窓で、再抽選前の重み付きアンサンブル平均 x̄・座標別標準偏差 σ
    # (状態+拡大の全行)を用い、再抽選後の各粒子を
    # x′ = a·x + (1−a)·x̄ + √(1−a²)·σ·ε(ε〜N(0,I)) で置換する
    # (1次・2次モーメント保存)。既定 1.0 = 若返り無効(従来動作と完全同一、
    # 統計計算もスキップ)。EnKS ラグ窓スナップショットは対象外
    # (インデックス再抽選のみ、過去断面への遡及ジッタはしない)。
    rejuvenation_a::Float64 = 1.0
end

"同化ランの結果(X は状態行 × 時刻 × メンバー。拡大時は N_STATE+1 行目以降が拡大パラメータ、#0046)"
struct AssimResult
    t::Vector{Float64}
    X::Array{Float64,3}
    ranks::Dict{Symbol,Vector{Int}}      # 解析直前の順位(ランクヒストグラム用)
    ess::Vector{Float64}                 # 各週次窓の ESS
    nresample::Int
    ts_snap::Vector{Float64}             # EnKS スナップショット時刻(平滑化オフなら空)
    Xs::Array{Float64,3}                 # 平滑化アンサンブル(行 × スナップ × メンバー)
    count_observed::Vector{Int}          # 各週次窓の観測カウント(-1 = データなし、#0054)
    count_logscore::Vector{Float64}      # 各週次窓の1ステップ先予測 log スコア(NaN = データなし、#0054)
    # g 計装診断(DECISIONS #0065/#0066)。`instrument_g = false`(既定)なら
    # 常に空ベクトル(計算自体を行わない、後方互換・判定数値に無関係)。各要素は
    # Dict("t", "phase" [:pre/:post], "update_type"
    # [:count_weekly/:g_swiid_annual/:other_obs], "g_mean", "g_sd",
    # "last_g_swiid_obs", "aug")。"aug" は拡大パラメータ別の
    # Dict(name => Dict("link", "mean", "sd"))(#0066。統計は**内部リンク座標**
    # — :log リンクなら log(自然値) の平均・sd、:identity なら自然値そのまま。
    # 拡大なしのランでは "aug" キー自体を省略)。読み取り専用の診断で、
    # 乱数・状態・重みは変更しない。
    g_diag::Vector{Dict{String,Any}}
end

"""
    pathological(xi) -> Bool

メンバーが数値的に病的な領域にいるか(DECISIONS #0011)。
§9.4 の警告水準(|logit ξ| > 10)を大きく超えた |logit ξ| > 15、
または σ_s > e³ ≈ 20(降伏応力 σ_Y = 1 の20倍)、または非有限値。
"""
function pathological(xi::AbstractVector{Float64})
    all(isfinite, xi) || return true
    for i in (IX_W, IX_G, IX_PHI, IX_TAU, IX_PP)
        abs(xi[i]) > 15 && return true
    end
    return xi[IX_SIG] > 3
end

"""
    liu_west_rejuvenate!(E, xbar, sigma, a, rng) -> E

Liu-West 縮小核による若返り(DECISIONS #0057)。系統再抽選**後**の `E`
(全行 × N メンバー)の各粒子列を
`x′ = a·x + (1−a)·xbar + √(1−a²)·sigma·ε`(ε〜N(0,I))で置換する
(再抽選前の重み付き平均 `xbar`・座標別標準偏差 `sigma` を用い、1次・2次
モーメントを保存)。`sigma` が 0 または非有限の座標はジッタしない
(その行は `a·x + (1−a)·xbar` のみ、`a=1` かつ `xbar` 無関係なら不変)。
"""
function liu_west_rejuvenate!(E::AbstractMatrix{Float64}, xbar::AbstractVector{Float64},
                               sigma::AbstractVector{Float64}, a::Float64, rng::AbstractRNG)
    n, N = size(E)
    b = sqrt(max(1 - a^2, 0.0))
    for i in 1:N
        for k in 1:n
            jitter = (isfinite(sigma[k]) && sigma[k] > 0) ? b * sigma[k] * randn(rng) : 0.0
            E[k, i] = a * E[k, i] + (1 - a) * xbar[k] + jitter
        end
    end
    return E
end

"""
    weighted_mean_std(E, w) -> (xbar, sigma)

重み `w`(正規化済み、和1)によるアンサンブル `E`(全行 × N メンバー)の
座標別加重平均 `xbar` と加重標準偏差 `sigma`(不偏補正なし、母標準偏差)を
返す。DECISIONS #0057 の若返り統計(再抽選前に評価)。
"""
function weighted_mean_std(E::AbstractMatrix{Float64}, w::AbstractVector{Float64})
    n, N = size(E)
    xbar = [sum(w[i] * E[k, i] for i in 1:N) for k in 1:n]
    var = [sum(w[i] * (E[k, i] - xbar[k])^2 for i in 1:N) for k in 1:n]
    sigma = sqrt.(max.(var, 0.0))
    return xbar, sigma
end

"""
    select_masked_rows(cfg, batch) -> Vector{Int}

現在時刻解析の変数局所化(DECISIONS #0040-(α))の対象行を選ぶ。
`cfg.analysis_masked_vars` が空なら常に `Int[]`(既定・従来動作)。
非空でも、`batch` に `cfg.analysis_unmask_names` のいずれかの観測名が
1つでも含まれればマスク解除(`Int[]`)。それ以外は `cfg.analysis_masked_vars`
をそのまま返す。
"""
function select_masked_rows(cfg::AssimConfig, batch::AbstractVector{ObservationRecord})
    isempty(cfg.analysis_masked_vars) && return Int[]
    any(o.spec.name in cfg.analysis_unmask_names for o in batch) && return Int[]
    return cfg.analysis_masked_vars
end

"""
    apply_obs_spread_floor!(E, batch, floor_frac, rng) -> E

観測座標の加法的スプレッド床(DECISIONS #0043)。`floor_frac <= 0` は
no-op(既定・後方互換)。`batch` 内の各観測について、`target_ix == 0`
(合成観測、恒等写像でない h)は対象外。恒等写像観測(`target_ix != 0`)は
事前アンサンブル `E` での観測座標の sd を評価し、`floor = floor_frac *
o.spec.sd` を下回れば、対応する状態行(`E` の 1:N_STATE 内、augmented でも
theta_sig 行には触れない)へ独立ガウス摂動 N(0, floor² − sd²) を全メンバー
独立に加える。同一バッチで複数観測が同じ状態行を対象にする場合は、
必要な追加分散が最大の1回のみ適用する(二重加算を避ける)。
`rng` は呼び出し側(run_assimilation の既存ストリーム)から渡す。
"""
function apply_obs_spread_floor!(E::AbstractMatrix{Float64},
                                 batch::AbstractVector{ObservationRecord},
                                 floor_frac::Float64, rng::AbstractRNG)
    floor_frac > 0 || return E
    N = size(E, 2)
    extra_var = Dict{Int,Float64}()   # 状態行 → 追加分散(バッチ内最大)
    for o in batch
        ix = o.spec.target_ix
        ix == 0 && continue
        z = [o.spec.h(view(E, 1:N_STATE, j)) for j in 1:N]
        zbar = sum(z) / N
        sd = sqrt(sum(abs2, z .- zbar) / (N - 1))
        floor = floor_frac * o.spec.sd
        sd < floor || continue
        v = floor^2 - sd^2
        extra_var[ix] = max(get(extra_var, ix, 0.0), v)
    end
    for (ix, v) in extra_var
        sdadd = sqrt(v)
        for j in 1:N
            E[ix, j] += sdadd * randn(rng)
        end
    end
    return E
end

"theta_sig を差し替えた ModelParameters(状態拡大メンバー用。後方互換の特例経路)"
with_theta_sig(p::ModelParameters, theta::Real) =
    ModelParameters(p.regime, p.l1, p.l2, L3Params(theta_sig = float(theta)),
                    p.exo, p.x0_nat, p.x0)

"""
駆動パラメータの L3 状態拡大の記述子(DECISIONS #0046)。

`augmented::Bool` + `param_noise_sd`(theta_sig 専用のハードコード、#0010系)を
一般化し、状態拡大するパラメータの集合を記述子のリストとして表現する。

- `name`: L2Params / L3Params / ConstantExogenous のいずれかのフィールド名
  (`_aug_location` が名前だけで所属構造体を判定する)。
- `link`: `:log`(状態行 = log(自然値)。正値パラメータ用)または
  `:identity`(状態行 = 自然値。負値も可)。
- `init`: 初期値(自然単位)。初期アンサンブル構築(`augment_ensemble`)にのみ使う。
- `init_sd`: 初期アンサンブルスプレッド(リンク座標)。
- `rw_sd`: 予測ステップのランダムウォーク sd(リンク座標、/√年)。
"""
Base.@kwdef struct AugmentedParam
    name::Symbol
    link::Symbol = :identity
    init::Float64 = 0.0
    init_sd::Float64 = 0.0
    rw_sd::Float64 = 0.0
end

"リンク座標 → 自然単位(#0046)。未知の link はエラー"
function _link_from(link::Symbol, x::Real)
    link === :log && return exp(x)
    link === :identity && return x
    throw(ArgumentError("unknown AugmentedParam link :$link (expected :log or :identity)"))
end

"自然単位 → リンク座標(#0046)。未知の link はエラー"
function _link_to(link::Symbol, x::Real)
    link === :log && return log(x)
    link === :identity && return x
    throw(ArgumentError("unknown AugmentedParam link :$link (expected :log or :identity)"))
end

"""
    _aug_location(name) -> Symbol

拡大パラメータ名がどの構造体に属するか(:l2 / :l3 / :exo)を、フィールド名
だけから判定する(DECISIONS #0046)。L3Params → ConstantExogenous → L2Params
の順で探し、いずれにも無ければエラー(現行モデルは名前衝突なし)。
"""
function _aug_location(name::Symbol)
    name in fieldnames(L3Params) && return :l3
    name in fieldnames(ConstantExogenous) && return :exo
    name in fieldnames(L2Params) && return :l2
    throw(ArgumentError("augmented param :$name not found in L2Params/L3Params/ConstantExogenous"))
end

"""
    _with_field(x::T, field, value) -> T

kwdef struct `x` のフィールド `field` だけを `value` に差し替えたコピー
(他フィールドはそのままコピー)。Setfield 等の新規依存を避けるための
手書きヘルパ(DECISIONS #0046)。
"""
function _with_field(x::T, field::Symbol, value) where {T}
    vals = Dict{Symbol,Any}(f => getfield(x, f) for f in fieldnames(T))
    haskey(vals, field) || throw(ArgumentError("$T has no field :$field"))
    vals[field] = value
    return T(; vals...)
end

"""
    _inject_param(p::ModelParameters, name, value) -> ModelParameters

拡大パラメータ1個(自然単位の値)を注入した `ModelParameters` を返す
(所属構造体は `_aug_location` で判定)。
"""
function _inject_param(p::ModelParameters, name::Symbol, value::Float64)
    loc = _aug_location(name)
    if loc === :l3
        return ModelParameters(p.regime, p.l1, p.l2, _with_field(p.l3, name, value),
                               p.exo, p.x0_nat, p.x0)
    elseif loc === :exo
        return ModelParameters(p.regime, p.l1, p.l2, p.l3, _with_field(p.exo, name, value),
                               p.x0_nat, p.x0)
    else
        return ModelParameters(p.regime, p.l1, _with_field(p.l2, name, value), p.l3, p.exo,
                               p.x0_nat, p.x0)
    end
end

"""
    build_member_params(params, augmented_params, E, state_rows, i) -> ModelParameters

拡大行(`E` の `state_rows+1:state_rows+length(augmented_params)`、メンバー `i`)
をそれぞれのリンク座標から自然単位に逆変換して `params` に順次注入する
(DECISIONS #0046)。`augmented_params` が空なら `params` をそのまま返す。
"""
function build_member_params(params::ModelParameters,
                             augmented_params::Vector{AugmentedParam},
                             E::AbstractMatrix{Float64}, state_rows::Int, i::Int)
    p = params
    for (k, ap) in enumerate(augmented_params)
        value = _link_from(ap.link, E[state_rows + k, i])
        p = _inject_param(p, ap.name, value)
    end
    return p
end

"""
    augment_ensemble(E0_state, augmented_params; rng) -> Matrix

状態行列 `E0_state`(N_STATE × N)に `augmented_params` の初期アンサンブル行
(記述子順)を追加した拡大初期アンサンブル(n × N、n = N_STATE +
length(augmented_params))を返す(DECISIONS #0046)。行 k は
`_link_to(link, init) + init_sd * randn()`(メンバー独立)。
`augmented_params` が空なら `E0_state` のコピーをそのまま返す(後方互換)。
"""
function augment_ensemble(E0_state::AbstractMatrix{Float64},
                         augmented_params::Vector{AugmentedParam};
                         rng::AbstractRNG)
    isempty(augmented_params) && return Matrix{Float64}(E0_state)
    N = size(E0_state, 2)
    extra = Matrix{Float64}(undef, length(augmented_params), N)
    for (k, ap) in enumerate(augmented_params)
        c = _link_to(ap.link, ap.init)
        extra[k, :] .= c .+ ap.init_sd .* randn(rng, N)
    end
    return vcat(Matrix{Float64}(E0_state), extra)
end

"アンサンブル行1本の加重なし平均・標本標準偏差(不偏補正、N>1。読み取り専用)"
function _row_mean_sd(row::AbstractVector{Float64})
    n = length(row)
    m = sum(row) / n
    sd = n > 1 ? sqrt(sum(abs2, row .- m) / (n - 1)) : 0.0
    return m, sd
end

"""
    g_diag_entry(t, phase, update_type, E, last_g_swiid_val,
                 aug_params = AugmentedParam[]) -> Dict

g 計装診断1レコード(DECISIONS #0065/#0066)。`E`(全行 × N メンバー)の
`IX_G` 行と各拡大パラメータ行(`aug_params` の記述子順、N_STATE+1 以降)から
読み取り専用で加重なし平均・標本標準偏差(不偏補正、N>1)を計算するだけで、
`E` そのものは一切変更しない。拡大パラメータの統計は**内部リンク座標**
(`:log` なら log(自然値))で記録し、各エントリに "link" を明記する
(#0066。`aug_params` が空なら "aug" キー自体を省略)。`last_g_swiid_val` は
直近に観測された g_swiid 値(まだ観測が無ければ `nothing`)。
"""
function g_diag_entry(t::Float64, phase::Symbol, update_type::Symbol,
                      E::AbstractMatrix{Float64},
                      last_g_swiid_val::Union{Nothing,Float64},
                      aug_params::Vector{AugmentedParam} = AugmentedParam[])
    gmean, gsd = _row_mean_sd(view(E, IX_G, :))
    entry = Dict{String,Any}(
        "t" => t,
        "phase" => String(phase),
        "update_type" => String(update_type),
        "g_mean" => gmean,
        "g_sd" => gsd,
        "last_g_swiid_obs" => last_g_swiid_val,
    )
    if !isempty(aug_params)
        augd = Dict{String,Any}()
        for (k, ap) in enumerate(aug_params)
            m, sd = _row_mean_sd(view(E, N_STATE + k, :))
            augd[String(ap.name)] = Dict{String,Any}(
                "link" => String(ap.link), "mean" => m, "sd" => sd)
        end
        entry["aug"] = augd
    end
    return entry
end

"""
    g_aug_series(ts, X, aug_params; stride = 2,
                 update_type = :forecast_free) -> Vector{Dict}

軌道配列 `X`(全行 × 時刻 × メンバー、`run_assimilation`/`forecast_ensemble`
の出力形式)から、`stride` グリッド刻み(既定 2 = dt 0.01 の週次窓 0.02 年)
で g・拡大パラメータのアンサンブル平均・sd の時系列を後処理で抽出する
(DECISIONS #0066 の予報窓計装)。完全に読み取り専用(乱数消費・副作用なし)
で、実行済み軌道の集計のみ — 判定数値には原理的に影響しない。レコード形式は
`g_diag_entry` と同一(phase = "series"、last_g_swiid_obs = nothing)。
末端時刻は stride に一致しなくても必ず含める。
"""
function g_aug_series(ts::AbstractVector{Float64}, X::AbstractArray{Float64,3},
                      aug_params::Vector{AugmentedParam}; stride::Int = 2,
                      update_type::Symbol = :forecast_free)
    stride >= 1 || throw(ArgumentError("stride must be >= 1"))
    out = Dict{String,Any}[]
    idxs = collect(1:stride:length(ts))
    isempty(idxs) && return out
    idxs[end] != length(ts) && push!(idxs, length(ts))
    for s in idxs
        push!(out, g_diag_entry(ts[s], :series, update_type,
                                view(X, :, s, :), nothing, aug_params))
    end
    return out
end

"""
    run_assimilation(params, E0, obs, event_times; cfg, seed,
                     augmented=false, augmented_params=AugmentedParam[],
                     obs_counts=nothing, count_scale=1.0)
        -> AssimResult

初期アンサンブル `E0`(n × N。n = N_STATE + 拡大パラメータ数)から §9 の
ハイブリッド同化(EnKF + ポアソン重み + イベント同期)を実行する。

駆動パラメータの L3 状態拡大(DECISIONS #0046): `augmented_params` に
`AugmentedParam` のリストを渡すと、その記述子順に状態行(N_STATE+1 以降)を
解釈し、各ステップでリンク座標のランダムウォークを加え、`ModelParameters`
への注入は名前で L2Params/L3Params/ConstantExogenous に振り分ける。
`augmented::Bool`(既定 false)は theta_sig 1個のみを状態拡大する従来経路
(#0010 系)で、内部的には `AugmentedParam(:theta_sig, :log, ..., rw_sd =
param_noise_sd)` 1個の記述子に変換して同じ機構で処理する(結果は従来と
同一。E1/既存テストの記録結果を保護)。`augmented` と `augmented_params` の
同時指定はエラー。

週次イベントカウントの扱い(DECISIONS #0031):
- 既定(`obs_counts = nothing`): E1 と同じく `event_times`(真値カタログ)を
  窓に集計して観測カウントとする(モデルジャンプ = 観測イベントが1対1)。
- 実データ(M8): `obs_counts` に窓別の観測カウント列(窓 k は区間
  [t0+(k−1)·event_window, t0+k·event_window))を渡す。**負値はデータなし**を
  意味し、その窓のポアソン重み更新をスキップする(#0031-3)。
  `count_scale` は報告率 ν(N_w 〜 Poisson(ν·Λ)、#0031-1)。`count_temper` は
  過分散カウントの尤度テンパリング係数(1/ν 推奨、#0033。既定 1 = 素のポアソン)。

カウント尤度モデルの切り替え(DECISIONS #0054):
- `count_model = :poisson`(既定): 従来どおり `count_temper` 付きポアソン重み。
- `count_model = :negbin`: NegBin(平均 ν·Λ_i、サイズ `count_dispersion`)の
  完全 log pmf を重みに使う(過分散を分布で表現するためテンパリングは
  適用しない — `count_temper` は無視される)。
診断: 各週次窓のカウントと1ステップ先予測 log スコア(重み付き混合予測
分布 log Σ_i w_i p(N_w | member i)。p は当該 `count_model` の素の pmf)を
`AssimResult.count_observed` / `count_logscore` に記録する。

g 計装診断(DECISIONS #0065/#0066): `instrument_g = true`(既定 false)に
すると、各週次カウント再重み付け(再抽選の前後、`update_type =
:count_weekly`)と各解析ステップ(`enks_analysis!` 前後、バッチに `:g_swiid`
観測を含めば `update_type = :g_swiid_annual`、それ以外は `:other_obs`)の
前後で、アンサンブルの `g`(`IX_G` 行)の平均・sd・直近の g_swiid 観測値と、
各拡大パラメータ行のメンバー平均・sd(#0066。**内部リンク座標**、各エントリに
"link" を明記)を `AssimResult.g_diag` に記録する。**読み取り専用**(乱数
消費・状態・重みへの副作用は一切ない)なので、既定 false のときは元コードと
完全に同一の計算(計装のオーバーヘッドすらない)、true でも同化の数値結果は
不変。
"""
function run_assimilation(params::ModelParameters, E0::Matrix{Float64},
                          obs::Vector{ObservationRecord},
                          event_times::Vector{Float64};
                          cfg::AssimConfig = AssimConfig(), seed::Integer,
                          augmented::Bool = false,
                          augmented_params::Vector{AugmentedParam} = AugmentedParam[],
                          obs_counts::Union{Nothing, Vector{Int}} = nothing,
                          count_scale::Float64 = 1.0,
                          count_temper::Float64 = 1.0,
                          count_model::Symbol = :poisson,
                          count_dispersion::Float64 = Inf,
                          instrument_g::Bool = false)
    count_model in (:poisson, :negbin) ||
        throw(ArgumentError("unknown count_model :$count_model (expected :poisson or :negbin)"))
    count_model === :negbin && !(isfinite(count_dispersion) && count_dispersion > 0) &&
        throw(ArgumentError("count_model = :negbin requires finite positive count_dispersion"))
    augmented && !isempty(augmented_params) &&
        throw(ArgumentError("augmented=true と augmented_params の同時指定はできません"))
    # 後方互換(#0010系): augmented=true は theta_sig 1個の記述子に変換する
    aug_params = augmented ?
        [AugmentedParam(name = :theta_sig, link = :log, rw_sd = cfg.param_noise_sd)] :
        augmented_params

    n, N = size(E0)
    n == N_STATE + length(aug_params) ||
        throw(DimensionMismatch("E0 has $n rows, expected $(N_STATE + length(aug_params)) " *
                                "(N_STATE=$N_STATE + $(length(aug_params)) augmented params)"))

    nsteps = round(Int, (cfg.t1 - cfg.t0) / cfg.dt)
    ts = collect(range(cfg.t0; step = cfg.dt, length = nsteps + 1))
    grid_index(t) = clamp(round(Int, (t - cfg.t0) / cfg.dt) + 1, 1, nsteps + 1)

    # 観測をグリッド点にグループ化(§9.2: 届いた時刻に届いた観測だけで解析)
    obs_at = Dict{Int,Vector{ObservationRecord}}()
    for o in obs
        push!(get!(obs_at, grid_index(o.t), ObservationRecord[]), o)
    end

    events = sort(event_times)
    next_ev = 1

    # 週次窓の境界グリッドと観測カウント
    wsteps = max(1, round(Int, cfg.event_window / cfg.dt))

    E = copy(E0)
    X = Array{Float64,3}(undef, n, nsteps + 1, N)
    X[:, 1, :] = E
    rngs = [Xoshiro(member_seed(seed, i)) for i in 1:N]
    # ランク計算専用の独立ストリーム(観測ノイズ抽選が力学の乱数列を乱さないように)
    rank_rng = Xoshiro(member_seed(seed, 10_000_019))
    f = Vector{Float64}(undef, N_STATE)
    sig = Vector{Float64}(undef, N_STATE)
    dW = Vector{Float64}(undef, N_STATE)
    Lambda = zeros(N)
    logw = zeros(N)          # 累積 log 重み(再抽選までウィンドウ間で持ち越す)
    window_count = 0
    ess_hist = Float64[]
    count_observed_hist = Int[]        # #0054 診断
    count_logscore_hist = Float64[]    # #0054 診断
    nresample = 0
    jump_since_analysis = false
    t_last_analysis = cfg.t0
    ranks = Dict{Symbol,Vector{Int}}()
    sqdt = sqrt(cfg.dt)
    g_diag = Dict{String,Any}[]              # #0065 診断(instrument_g = false なら空のまま)
    last_g_swiid = Ref{Union{Nothing,Float64}}(nothing)

    member_params(i) = isempty(aug_params) ? params :
        build_member_params(params, aug_params, E, N_STATE, i)

    # 固定ラグ EnKS(#0024): smoother_dt 刻みでスナップショットを保持し、
    # 解析のたびに現在時刻から smoother_lag 以内のものを同時更新する。
    smoothing = cfg.smoother_lag > 0
    snap_steps = max(1, round(Int, cfg.smoother_dt / cfg.dt))
    snap_ts = Float64[]
    snaps = Matrix{Float64}[]
    lag_start = 1                          # ラグ窓内の最初のスナップショット index
    if smoothing
        push!(snap_ts, ts[1])
        push!(snaps, copy(E))
    end

    for step in 1:nsteps
        t = ts[step]
        t_next = ts[step + 1]

        # (a) 強制ジャンプ(イベント同期、§9.3)
        while next_ev <= length(events) && events[next_ev] < t_next
            for i in 1:N
                xi = @view E[1:N_STATE, i]
                rho = draw_mark(rngs[i], params)
                apply_jump!(xi, rho, member_params(i))
            end
            window_count += 1
            jump_since_analysis = true
            next_ev += 1
        end

        # (b) Λ トラッキング(§9.3。直前状態で近似)と EM ステップ
        for i in 1:N
            p_i = member_params(i)
            xi = @view E[1:N_STATE, i]
            Lambda[i] += intensity(xi, p_i) * cfg.dt
            drift!(f, xi, p_i, t)
            guard_sigma_drift!(f)                    # σ_s ガード(#0032)
            diffusion!(sig, xi, p_i, t)
            randn!(rngs[i], dW)
            @. xi += cfg.dt * f + sqdt * sig * dW
            guard_sigma_state!(xi)
            for (k, ap) in enumerate(aug_params)   # d(param) = 0 + 微小ノイズ(§8.3/#0046)
                E[N_STATE + k, i] += ap.rw_sd * sqdt * randn(rngs[i])
            end
        end

        # (c) 週次窓の終端: ポアソン重みを累積し、ESS < N/2 で系統再抽選(§9.3)
        if step % wsteps == 0
            # 観測カウント: 既定はカタログ集計(E1)、実データでは窓別列(#0031)。
            # 負値 = データなし窓 → 重み更新スキップ(病的ガードは常時)。
            widx = step ÷ wsteps
            observed = obs_counts === nothing ? window_count :
                       (widx <= length(obs_counts) ? obs_counts[widx] : -1)
            push!(count_observed_hist, observed)
            if observed >= 0
                mus = count_scale .* Lambda
                # 1ステップ先予測 log スコア(#0054 診断): 更新前の重み付き
                # 混合予測分布 log Σ_i w_i p(N_w | member i)。p は当該
                # count_model の素の pmf(テンパリングは分布でないため不使用)
                wpre = normalize_weights(logw)
                lp = count_model === :negbin ?
                     negbin_logweights(observed, mus, count_dispersion) :
                     [poisson_logpmf(observed, mu) for mu in mus]
                lterms = log.(wpre) .+ lp
                mx = maximum(lterms)
                push!(count_logscore_hist,
                      isfinite(mx) ? mx + log(sum(exp, lterms .- mx)) : -Inf)
                # 重み更新: :negbin は完全 log pmf(#0054)、:poisson は従来の
                # count_temper 付き(#0033。既定1 = 素のポアソン)
                if count_model === :negbin
                    logw .+= lp
                else
                    logw .+= count_temper .* poisson_logweights(observed, mus)
                end
            else
                push!(count_logscore_hist, NaN)
            end
            # 病的メンバーは重みゼロ化して強制再抽選(#0011)。ESS は単一
            # 外れ値では下がらないため、暴走メンバーが強制ジャンプ
            # (m ∝ sigma_s^-)で数値爆発する前に淘汰する必要がある。
            npath = 0
            for i in 1:N
                if pathological(view(E, 1:N_STATE, i))
                    logw[i] = -Inf
                    npath += 1
                end
            end
            npath < N || error("filter diverged: all members pathological")
            w = normalize_weights(logw)
            essval = ess(w)
            push!(ess_hist, essval)
            instrument_g &&
                push!(g_diag, g_diag_entry(t_next, :pre, :count_weekly, E, last_g_swiid[], aug_params))
            if npath > 0 || essval < N * cfg.ess_ratio
                # 若返り統計(DECISIONS #0057)は再抽選**前**の重み付きアンサンブルで
                # 評価する。a == 1.0(既定・若返り無効)なら計算をスキップして
                # 従来動作と完全同一(ビット一致)にする。
                do_rejuvenate = cfg.rejuvenation_a < 1.0
                xbar, sigma = do_rejuvenate ? weighted_mean_std(E, w) :
                              (Float64[], Float64[])
                idx = systematic_resample(rngs[1], w)
                E .= E[:, idx]
                # メンバー対応を保つため、ラグ窓内のスナップショットも同じ
                # インデックスで再抽選する(EnKS、#0024)。ジッタは加えない
                # (過去断面への遡及ジッタはしない、#0057)。
                for s in lag_start:length(snaps)
                    snaps[s] .= snaps[s][:, idx]
                end
                if do_rejuvenate
                    liu_west_rejuvenate!(E, xbar, sigma, cfg.rejuvenation_a, rngs[1])
                end
                fill!(logw, 0.0)
                nresample += 1
            end
            instrument_g &&
                push!(g_diag, g_diag_entry(t_next, :post, :count_weekly, E, last_g_swiid[], aug_params))
            fill!(Lambda, 0.0)
            window_count = 0
        end

        # (d) 解析ステップ(この時刻に届いた観測のみ、§9.2)
        if haskey(obs_at, step + 1)
            batch = obs_at[step + 1]
            # tauA への緩い擬似観測(DECISIONS #0036、既定オフ)。batch を
            # コピーして追加するため obs_at 由来の元配列は変更しない。
            if cfg.tauA_pseudo_sd_mult > 0
                batch = augment_tauA_pseudo(batch, cfg.tauA_pseudo_sd_mult)
            end
            # g 計装診断(DECISIONS #0065): batch 確定後・E への一切の変更前に
            # 読み取り専用で記録する(instrument_g = false なら no-op)。
            if instrument_g
                has_g = any(o.spec.name === :g_swiid for o in batch)
                if has_g
                    last_g_swiid[] = first(o.value for o in batch if o.spec.name === :g_swiid)
                end
                utype = has_g ? :g_swiid_annual : :other_obs
                push!(g_diag, g_diag_entry(t_next, :pre, utype, E, last_g_swiid[], aug_params))
            end
            # 観測座標の事前スプレッド床(DECISIONS #0043、既定オフ)。batch
            # 確定後・yobs/R/hfun 組み立て前に E を直接摂動するため、以降の
            # ランク計算・spread_prior(RTPS)・enks_analysis! は全て床適用後の
            # 実効事前を見る(意図通り: #0043 は解析直前の実効事前の拡大)。
            if cfg.obs_spread_floor_frac > 0
                apply_obs_spread_floor!(E, batch, cfg.obs_spread_floor_frac, rngs[1])
            end
            # ランク(解析直前の事前アンサンブルに対する観測の順位)。
            # 観測 = 真値 + ノイズ のため、メンバー側にも観測ノイズ抽選を
            # 加えるのがランクヒストグラムの標準定義(Hamill 2001、#0017)。
            # 省くと高頻度観測変数で見かけの過小分散が生じる。
            for o in batch
                yj = [o.spec.h(view(E, 1:N_STATE, j)) + o.spec.sd * randn(rank_rng)
                      for j in 1:N]
                push!(get!(ranks, o.spec.name, Int[]),
                      count(<(o.value), yj) + 1)
            end
            yobs = [o.value for o in batch]
            R = Diagonal([o.spec.sd^2 for o in batch]) |> Matrix
            hfun = col -> [o.spec.h(view(col, 1:N_STATE)) for o in batch]

            # 現在時刻解析の変数局所化(#0040-(α))
            masked_rows = select_masked_rows(cfg, batch)

            # ラグ窓の前進(EnKS。窓外のスナップショットは確定)
            while smoothing && lag_start <= length(snap_ts) &&
                  snap_ts[lag_start] < t_next - cfg.smoother_lag
                lag_start += 1
            end
            # 平滑化は疎な観測を含む解析のみ(smoother_min_period)
            do_smooth = smoothing &&
                any(o.spec.period >= cfg.smoother_min_period for o in batch)
            window_snaps = do_smooth ? view(snaps, lag_start:length(snaps)) :
                           Matrix{Float64}[]

            # スプレッド注入(inflation_mode、DECISIONS #0013)。
            # 平滑化更新は現在状態と同一のイノベーションで行う(#0024)。
            if cfg.inflation_mode === :rtps
                spread_prior = ensemble_spread(E)
                enks_analysis!(E, window_snaps, yobs, hfun, R;
                               rng = rngs[1], rho_inf = 1.0,
                               smooth_rows = cfg.smoother_vars,
                               masked_rows)
                alpha = jump_since_analysis ? cfg.rtps_alpha_jump : cfg.rtps_alpha
                rtps!(E, spread_prior; alpha)
            else
                rho_base = jump_since_analysis ? cfg.rho_inf_jump : cfg.rho_inf
                rho = cfg.inflation_mode === :per_time ?
                    rho_base^((t_next - t_last_analysis) / cfg.tau_ref) : rho_base
                enks_analysis!(E, window_snaps, yobs, hfun, R;
                               rng = rngs[1], rho_inf = rho,
                               smooth_rows = cfg.smoother_vars,
                               masked_rows)
            end
            postprocess_analysis!(E)
            jump_since_analysis = false
            t_last_analysis = t_next
            instrument_g &&
                push!(g_diag, g_diag_entry(t_next, :post, utype, E, last_g_swiid[], aug_params))
        end

        # EnKS スナップショットの追加(解析後の状態、snap_dt 刻み)
        if smoothing && step % snap_steps == 0
            push!(snap_ts, t_next)
            push!(snaps, copy(E))
        end

        X[:, step + 1, :] = E
    end

    if smoothing
        Xs = Array{Float64,3}(undef, n, length(snaps), N)
        for (s, S) in enumerate(snaps)
            Xs[:, s, :] = S
        end
        return AssimResult(ts, X, ranks, ess_hist, nresample, snap_ts, Xs,
                           count_observed_hist, count_logscore_hist, g_diag)
    end
    return AssimResult(ts, X, ranks, ess_hist, nresample,
                       Float64[], Array{Float64,3}(undef, 0, 0, 0),
                       count_observed_hist, count_logscore_hist, g_diag)
end

"""
    free_ensemble(params, E0; cfg, seed) -> Array{Float64,3}

同化オフの自由ラン(同じ初期アンサンブル、内生 Hawkes、§13 手順5)。
戻り値は N_STATE × 時刻 × メンバー。
"""
function free_ensemble(params::ModelParameters, E0::Matrix{Float64};
                       cfg::AssimConfig = AssimConfig(), seed::Integer)
    N = size(E0, 2)
    nsteps = round(Int, (cfg.t1 - cfg.t0) / cfg.dt)
    X = Array{Float64,3}(undef, N_STATE, nsteps + 1, N)
    Threads.@threads for i in 1:N
        r = simulate_sde(params; seed = member_seed(seed, i),
                         t0 = cfg.t0, t1 = cfg.t1, dt = cfg.dt,
                         xi0 = E0[1:N_STATE, i])
        X[:, :, i] = r.traj.X
    end
    return X
end

"""
    simulate_sde_augmented(params, augmented_params; seed, t0, t1, dt=0.01, xi0)
        -> (; t, X, jumps)

拡大パラメータ込みの単一軌道前進(内生 Hawkes、DECISIONS #0053)。
`run_assimilation` の予測ステップと同一の生成則 — 各ステップで拡大行の
現在値から `build_member_params` でメンバーパラメータを構築して
ジャンプ・drift・diffusion を評価し、拡大行を `rw_sd·√dt` のランダム
ウォークで前進する — を、同化なしの純予報として実行する。M9 の1年先
予報アンサンブル(拡大事後値を初期値として持ち込み、予報区間中も RW を
継続する)向け。`augmented_params` が空なら `simulate_sde`(既定モード)と
同一の乱数消費順序になり、軌道はシード単位で一致する。

`xi0` は長さ `N_STATE + length(augmented_params)`(拡大行はリンク座標)。
戻り値 `X` も同じ行数(拡大行の軌道を含む)。
"""
function simulate_sde_augmented(params::ModelParameters,
                                augmented_params::Vector{AugmentedParam};
                                seed::Integer,
                                t0::Float64 = 0.0, t1::Float64 = 50.0,
                                dt::Float64 = 0.01,
                                xi0::AbstractVector{Float64})
    n = N_STATE + length(augmented_params)
    length(xi0) == n ||
        throw(DimensionMismatch("xi0 has $(length(xi0)) rows, expected $n " *
                                "(N_STATE=$N_STATE + $(length(augmented_params)) augmented params)"))
    rng = Xoshiro(seed)
    nsteps = round(Int, (t1 - t0) / dt)
    ts = collect(range(t0; step = dt, length = nsteps + 1))
    X = Matrix{Float64}(undef, n, nsteps + 1)
    jumps = JumpEvent[]
    x = collect(Float64, xi0)
    Ecol = reshape(x, :, 1)               # build_member_params 用(メモリ共有)
    xi = view(x, 1:N_STATE)
    f = Vector{Float64}(undef, N_STATE)
    sig = Vector{Float64}(undef, N_STATE)
    dW = Vector{Float64}(undef, N_STATE)
    X[:, 1] = x
    sqdt = sqrt(dt)

    for step in 1:nsteps
        t = ts[step]
        t_next = ts[step + 1]
        p = build_member_params(params, augmented_params, Ecol, N_STATE, 1)

        # (a) 内生ジャンプ(Ogata thinning。現在の拡大値のパラメータで評価)
        _jumps_in_interval!(xi, t, t_next, p, rng, jumps)

        # (b) EM ステップ(σ_s ガード #0032)+ 拡大行 RW(§8.3/#0046 と同一則)
        drift!(f, xi, p, t)
        guard_sigma_drift!(f)
        diffusion!(sig, xi, p, t)
        randn!(rng, dW)
        @. xi += dt * f + sqdt * sig * dW
        guard_sigma_state!(xi)
        for (k, ap) in enumerate(augmented_params)
            x[N_STATE + k] += ap.rw_sd * sqdt * randn(rng)
        end
        X[:, step + 1] = x
    end
    return (; t = ts, X, jumps)
end
