# 拡散項(SPEC §5)
#
# dW は独立な標準ウィーナー過程(対角拡散)。式(1)(6)(10)(13)は拡散なし。
# 式(12)のみ状態依存: sigma_p0 * (v / v0)^kappa_pv(v が分極ノイズを増幅)。

"""
    diffusion!(sig, xi, params, t)

状態 `xi`(保持座標)における対角拡散係数を `sig` に書き込む。
"""
function diffusion!(sig::AbstractVector, xi::AbstractVector,
                    params::ModelParameters, t::Real)
    l1, l2 = params.l1, params.l2
    v = exp(xi[IX_V])

    sig[IX_P] = 0.0                     # (1) 拡散なし
    sig[IX_W] = l2.sigma_w
    sig[IX_H] = l2.sigma_h
    sig[IX_K] = l2.sigma_k
    sig[IX_G] = l2.sigma_g
    sig[IX_T] = 0.0                     # (6) 拡散ゼロ(単調性の構成的保証)
    sig[IX_PHI] = l2.sigma_phi
    sig[IX_V] = l2.sigma_v
    sig[IX_TAU] = l2.sigma_tau
    sig[IX_TAUA] = 0.0                  # (10) 拡散なし
    sig[IX_SIG] = l2.sigma_sig
    sig[IX_PP] = l2.sigma_p0 * (v / l1.v0)^l2.kappa_pv   # (12) ★拡散側の連成
    sig[IX_LAME] = 0.0                  # (13) 純減衰+ジャンプのみ
    return sig
end
