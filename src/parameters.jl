# パラメータ構造体・事前分布・双子実験セット(SPEC §8)
#
# L1: 文献値・外部推定で固定(同化しない)
# L2: 事前分布 + ヒンドキャストでオフライン較正(双子実験では仮値を真値として使用)
# L3: オンライン推定(状態拡大。最大3個まで)
#
# 分岐比 n = alpha0/beta < 1 の assert は parameters 側に置く(§0.5.1/§6.1)。
# De のレジーム assert(安定国 <1 / 変動国 >1)は §8.4 の指示による。

"L1: 文献値・外部推定で固定(SPEC §8.1)"
Base.@kwdef struct L1Params
    alpha::Float64 = 0.33     # 資本分配率
    delta::Float64 = 0.05     # 資本減耗率 /年
    s::Float64 = 0.28         # 貯蓄率
    beta::Float64 = 5.0       # Hawkes 減衰率 /年
    A0::Float64 = 1.0         # スケール定数(正規化)
    h0::Float64 = 1.0
    y_ref::Float64 = 1.0
    v0::Float64 = 1.0
    sigma_Y::Float64 = 1.0    # 降伏応力(定義により 1)
end

"L2: オフライン較正対象。既定値は安定国の仮値(SPEC §8.2)"
Base.@kwdef struct L2Params
    # 人口構造
    kappa_w::Float64 = 0.05
    sigma_w::Float64 = 0.01
    # 教育投資
    gamma_h::Float64 = 0.03
    chi::Float64 = 0.3
    sigma_h::Float64 = 0.01
    # 資本
    theta_tau::Float64 = 0.5
    sigma_k::Float64 = 0.02
    # 格差動学
    eta_gT::Float64 = 0.5
    kappa_g::Float64 = 0.1
    mu_gbar::Float64 = logit(0.33)
    sigma_g::Float64 = 0.02
    # 技術成長
    a_T::Float64 = 0.02
    c_T0::Float64 = 0.0
    eta_Ty::Float64 = 0.5
    # TFP 弾性
    theta_T::Float64 = 0.4
    theta_phi::Float64 = 0.2
    # 普及
    r_phi::Float64 = 0.15
    kappa_phiv::Float64 = 0.3
    kappa_phitau::Float64 = 0.3
    sigma_phi::Float64 = 0.05
    # 情報インフラ
    kappa_v::Float64 = 0.3
    c_v0::Float64 = 0.0
    sigma_v::Float64 = 0.05
    # 信頼動学
    kappa_tau::Float64 = 0.15
    eta_taup::Float64 = 0.1
    sigma_tau::Float64 = 0.02
    # アンカー時定数
    eps_A::Float64 = 0.02
    # 応力載荷
    eta_g::Float64 = 0.5
    g_c::Float64 = 0.35
    eta_p::Float64 = 0.3
    eta_y::Float64 = 2.0
    # 応力緩和・ゆらぎ
    delta_sig::Float64 = 0.10
    sigma_sig::Float64 = 0.05
    # 分極
    eta_pg::Float64 = 0.4
    kappa_p::Float64 = 0.3
    mu_p::Float64 = logit(0.3)
    sigma_p0::Float64 = 0.05
    kappa_pv::Float64 = 0.5
    # ジャンプ強度
    lam0::Float64 = 0.02
    alpha0::Float64 = 2.0
    kappa_alphav::Float64 = 0.5
    lam_bar::Float64 = 50.0
    theta_p::Float64 = 1.0
    theta_v::Float64 = 0.3
    # 残留率マーク
    a_rho::Float64 = 2.0
    b_rho::Float64 = 2.0
    # 衝撃係数
    c_tau::Float64 = 1.0
    c_star::Float64 = 0.3
    c_g::Float64 = 0.5
    c_k::Float64 = 0.1
end

"L3: オンライン推定対象の現在値(SPEC §8.3。初期値は事前分布の中央値)"
Base.@kwdef struct L3Params
    theta_sig::Float64 = 3.0  # 降伏鋭度(事前 LogNormal(log 3, 0.5))
end

"L2 事前分布の標準形 LogNormal(log(中央値), 0.5)(SPEC §8.2)"
prior_lognormal(median::Real) = LogNormal(log(median), 0.5)

"L3 の初期事前分布(SPEC §8.3)"
l3_priors() = (
    theta_sig = prior_lognormal(3.0),
    theta_p = prior_lognormal(1.0),   # E1b でのみ状態拡大に使う
)

"外生入力(SPEC §5/§8.4。双子実験では定数)"
Base.@kwdef struct ConstantExogenous
    netgrowth::Float64   # b(t) - d(t) + mig(t)
    wbar::Float64        # 人口構造ターゲット
end

"""
モデルパラメータ一式(双子実験の1カ国分)。

`x0_nat` は自然座標の初期条件(§8.4)、`x0` は保持座標(§3)。
"""
struct ModelParameters
    regime::Symbol             # :stable / :volatile
    l1::L1Params
    l2::L2Params
    l3::L3Params
    exo::ConstantExogenous
    x0_nat::Vector{Float64}
    x0::Vector{Float64}
end

# §8.4 双子実験の初期条件(自然座標。インデックスは §3 に対応)
function _initial_conditions(regime::Symbol)
    if regime === :stable
        # 安定国(日本風)
        return [1.0, 0.60, 1.0, 1.0, 0.33, 1.0, 0.5, 1.0,
                0.55, logit(0.55), 0.3, 0.30, 0.0]
    elseif regime === :volatile
        return [1.0, 0.62, 1.0, 1.0, 0.45, 1.0, 0.5, 1.0,
                0.35, logit(0.35), 0.8, 0.50, 0.0]
    end
    throw(ArgumentError("unknown regime $regime (:stable or :volatile)"))
end

# 変動国で L2 仮値が安定国と異なる項目(SPEC §8.2 右列)
const _VOLATILE_L2_OVERRIDES = (
    mu_gbar = logit(0.45),
    eta_g = 1.0,
    delta_sig = 0.05,
    mu_p = logit(0.5),
    lam0 = 0.10,
)

"""
    build_params(regime; l2_overrides...) -> ModelParameters

安定国(:stable)/変動国(:volatile)のパラメータセットを構築する。
構築時に分岐比 n = alpha0/beta < 1(§6.1)と、De のレジーム条件
(安定国 De < 1、変動国 De > 1、§8.4)を検証し、満たさなければ
ArgumentError を投げる。`l2_overrides` はテスト・感度分析用の上書き。
"""
function build_params(regime::Symbol; l2_overrides...)
    l1 = L1Params()
    base = regime === :volatile ? pairs(_VOLATILE_L2_OVERRIDES) : pairs((;))
    l2 = L2Params(; base..., l2_overrides...)
    l3 = L3Params()

    exo = regime === :stable ?
        ConstantExogenous(netgrowth = -0.003, wbar = 0.58) :
        ConstantExogenous(netgrowth = +0.015, wbar = 0.62)

    x0_nat = _initial_conditions(regime)
    x0 = to_state(x0_nat)

    # 安定条件: 分岐比 n < 1(§6.1。v=v0 基準の名目値、DECISIONS #0005)
    n = branching_ratio(l1, l2)
    n < 1 || throw(ArgumentError(
        "branching ratio n = alpha0/beta = $n must be < 1 (SPEC §6.1)"))

    # レジーム条件: De(初期状態で評価、§7/§8.4)
    De = deborah_number(l2, x0_nat[IX_G], x0_nat[IX_PP])
    if regime === :stable
        De < 1 || throw(ArgumentError(
            "stable regime requires De < 1, got De = $De (SPEC §8.4)"))
    else
        De > 1 || throw(ArgumentError(
            "volatile regime requires De > 1, got De = $De (SPEC §8.4)"))
    end

    return ModelParameters(regime, l1, l2, l3, exo, x0_nat, x0)
end
