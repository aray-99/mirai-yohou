# 診断量(SPEC §4・§7)
#
# 診断変数(A, y, mu_T 等)は状態には含めず、毎時刻代数的に計算する。
# yhat は「ノイズを除いたドリフト部分」として解析的に合成するため、
# ドリフト計算と一体で src/drift.jl 側で求める(§4)。

"全要素生産性 A(T, phi) = A0 * T^theta_T * exp(theta_phi * phi)(SPEC §4)"
tfp(l1::L1Params, l2::L2Params, T::Real, phi::Real) =
    l1.A0 * T^l2.theta_T * exp(l2.theta_phi * phi)

"一人当たり産出 y = A * k^alpha * (h * w)^(1 - alpha)(SPEC §4)"
output_y(l1::L1Params, A::Real, k::Real, h::Real, w::Real) =
    A * k^l1.alpha * (h * w)^(1 - l1.alpha)

"技術成長率 mu_T = a_T * softplus(c_T0 + eta_Ty * log(y / y_ref))(SPEC §5 式6)"
tech_growth(l1::L1Params, l2::L2Params, y::Real) =
    l2.a_T * softplus(l2.c_T0 + l2.eta_Ty * log(y / l1.y_ref))

"従属人口比率 dep_ratio = (1 - w) / w(SPEC §4)"
dep_ratio(w::Real) = (1 - w) / w

"""
    branching_ratio(l1, l2) -> Float64

分岐比 n = alpha0 / beta(SPEC §7)。イベント1件が誘発する子イベント数。
n ≥ 1 で連鎖爆発(禁止。v=v0 基準の名目値、DECISIONS #0005)。
"""
branching_ratio(l1::L1Params, l2::L2Params) = l2.alpha0 / l1.beta

"""
    deborah_number(l2, gbar, pbar) -> Float64

社会的デボラ数 De = ( eta_g*(gbar - g_c)_+ + eta_p*pbar ) / delta_sig(SPEC §7)。
gbar, pbar は評価点での代表値(初期状態で評価してよい)。
De < 1: 応力が忘却で緩和(安定国)。De > 1: 単調蓄積(変動国)。
"""
deborah_number(l2::L2Params, gbar::Real, pbar::Real) =
    (l2.eta_g * pluspart(gbar - l2.g_c) + l2.eta_p * pbar) / l2.delta_sig

"""
    hardening_ratio(l2) -> Float64

硬化比 H = c_star / c_tau(SPEC §7)。被害のうち恒久的な割合。
0 = 弾性、1 = 完全塑性。
"""
hardening_ratio(l2::L2Params) = l2.c_star / l2.c_tau

"モデルの「性格」を決める4つの無次元数(SPEC §7)"
struct DimensionlessNumbers
    De::Float64         # 社会的デボラ数
    n::Float64          # 分岐比
    H::Float64          # 硬化比
    theta_sig::Float64  # 降伏鋭度
end

"""
    dimensionless_numbers(params) -> DimensionlessNumbers

無次元数4量を計算する。De は初期状態(x0_nat)で評価する(SPEC §7)。
"""
function dimensionless_numbers(params::ModelParameters)
    return DimensionlessNumbers(
        deborah_number(params.l2, params.x0_nat[IX_G], params.x0_nat[IX_PP]),
        branching_ratio(params.l1, params.l2),
        hardening_ratio(params.l2),
        params.l3.theta_sig,
    )
end

function Base.show(io::IO, d::DimensionlessNumbers)
    print(io, "DimensionlessNumbers(De=", round(d.De, digits=4),
          ", n=", round(d.n, digits=4),
          ", H=", round(d.H, digits=4),
          ", theta_sig=", round(d.theta_sig, digits=4), ")")
end
