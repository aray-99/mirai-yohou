# ドリフト項(SPEC §5 式(1)〜(13))
#
# 力学は変換後の座標(§3)で直接定義されている。式(11)のみ自然座標の
# Maxwell 粘弾性を log 座標へ変換した最終形で、-0.5 sigma_sig^2 が伊藤補正項。
# yhat(産出成長率)はノイズを除いたドリフト部分として解析的に合成する(§4)。

"""
    drift_with_diagnostics!(f, xi, params, t) -> NamedTuple

状態 `xi`(保持座標)におけるドリフトを `f` に書き込み、診断変数
`(A, y, mu_T, yhat, load)` を返す(SPEC §4/§5)。
"""
function drift_with_diagnostics!(f::AbstractVector,
                                 xi::AbstractVector,
                                 params::ModelParameters, t::Real)
    l1, l2, exo = params.l1, params.l2, params.exo

    # 自然座標の値(診断・非線形項に使用)
    w = sigmoid(xi[IX_W])
    h = exp(xi[IX_H])
    k = exp(xi[IX_K])
    g = sigmoid(xi[IX_G])
    T = exp(xi[IX_T])
    phi = sigmoid(xi[IX_PHI])
    v = exp(xi[IX_V])
    tau = sigmoid(xi[IX_TAU])
    p_pol = sigmoid(xi[IX_PP])

    # 診断変数(§4)
    A = tfp(l1, l2, T, phi)
    y = output_y(l1, A, k, h, w)
    mu_T = tech_growth(l1, l2, y)

    ng = exo.netgrowth              # b(t) - d(t) + mig(t)(双子実験では定数、§8.4)

    # (1) 人口(拡散なし)
    f[IX_P] = ng
    # (2) 人口構造
    f[IX_W] = l2.kappa_w * (logit(exo.wbar) - xi[IX_W])
    # (3) 教育投資: hmax(y) = h0 * (y / y_ref)^chi
    hmax = l1.h0 * (y / l1.y_ref)^l2.chi
    f[IX_H] = l2.gamma_h * (1 - h / hmax)
    # (4) ソロー型+信頼摩擦
    f[IX_K] = l1.s * tau^l2.theta_tau * (y / k) - l1.delta - ng
    # (5) SBTC 載荷+平均回帰
    f[IX_G] = l2.eta_gT * mu_T - l2.kappa_g * (xi[IX_G] - l2.mu_gbar)
    # (6) 技術フロンティア(拡散ゼロ。mu_T > 0 で単調性を構成的に保証)
    f[IX_T] = mu_T
    # (7) 普及
    f[IX_PHI] = l2.r_phi * (v / l1.v0)^l2.kappa_phiv * tau^l2.kappa_phitau
    # (8) 情報インフラ
    f[IX_V] = l2.kappa_v * (xi[IX_T] + l2.c_v0 - xi[IX_V])
    # (9) 信頼: アンカー回帰 − 分極侵食
    f[IX_TAU] = l2.kappa_tau * (xi[IX_TAUA] - xi[IX_TAU]) - l2.eta_taup * p_pol
    # (10) アンカー(遅い移動平均。拡散なし)
    f[IX_TAUA] = l2.eps_A * (xi[IX_TAU] - xi[IX_TAUA])

    # yhat = d log y / dt のドリフト部分(§4。数値微分は使わない)
    yhat = l2.theta_T * mu_T +
           l2.theta_phi * phi * (1 - phi) * f[IX_PHI] +
           l1.alpha * f[IX_K] +
           (1 - l1.alpha) * (f[IX_H] + (1 - w) * f[IX_W])

    # (11) 社会的応力(log 座標。-0.5 sigma_sig^2 は伊藤補正項)
    load = l2.eta_g * pluspart(g - l2.g_c) + l2.eta_p * p_pol -
           l2.eta_y * pluspart(yhat)
    f[IX_SIG] = load * exp(-xi[IX_SIG]) - l2.delta_sig - 0.5 * l2.sigma_sig^2
    # (12) 分極(拡散側の連成 sigma_p0*(v/v0)^kappa_pv は diffusion.jl、M2)
    f[IX_PP] = l2.eta_pg * (g - l2.g_c) - l2.kappa_p * (xi[IX_PP] - l2.mu_p)
    # (13) Hawkes 興奮度(純減衰。増加はジャンプ時のみ、§6.3)
    f[IX_LAME] = -l1.beta * xi[IX_LAME]

    return (A = A, y = y, mu_T = mu_T, yhat = yhat, load = load)
end

"""
    drift!(f, xi, params, t)

ドリフトのみを `f` に書き込む(診断値は捨てる)。
"""
function drift!(f::AbstractVector, xi::AbstractVector,
                params::ModelParameters, t::Real)
    drift_with_diagnostics!(f, xi, params, t)
    return f
end
