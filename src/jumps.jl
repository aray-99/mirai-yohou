# ジャンプ過程(SPEC §6): Hawkes 型・Markov 埋め込み
#
# 強度 lam = min(lam_bar, lam_b(X) + lam_e)。lam_e の力学は式(13)+
# ジャンプ時の加算で、指数カーネル Hawkes と等価(rate は現在状態のみで計算可)。
# ジャンプ処理は積分器から分離し、内生(thinning)/外生(強制発火)の
# 2モードを切り替えられるようにする(§6.4)。

"ジャンプ生成モード(SPEC §6.4)"
abstract type JumpMode end

"内生: thinning 法で発火(純予測用)"
struct EndogenousHawkes <: JumpMode end

"外生: 観測イベント時刻に強制発火(同化用)"
struct ExogenousEvents <: JumpMode
    times::Vector{Float64}
end

"""
    lam_b(xi, params) -> Float64

状態依存の基底強度(SPEC §6.1):
lam_b = lam0 * exp( theta_sig * (sigma_s - 1) + theta_p * p + theta_v * log(v / v0) )
"""
function lam_b(xi::AbstractVector, params::ModelParameters)
    l1, l2, l3 = params.l1, params.l2, params.l3
    sigma_s = exp(xi[IX_SIG])
    p = sigmoid(xi[IX_PP])
    return l2.lam0 * exp(l3.theta_sig * (sigma_s - 1) + l2.theta_p * p +
                         l2.theta_v * (xi[IX_V] - log(l1.v0)))
end

"""
    intensity(xi, params) -> Float64

全ジャンプ強度 lam = min(lam_bar, lam_b(X) + lam_e)(SPEC §6.1)。
"""
intensity(xi::AbstractVector, params::ModelParameters) =
    min(params.l2.lam_bar, lam_b(xi, params) + xi[IX_LAME])

"""
    draw_mark(rng, params) -> Float64

残留率 rho ~ Beta(a_rho, b_rho) を引く(SPEC §6.2)。
"""
draw_mark(rng::AbstractRNG, params::ModelParameters) =
    rand(rng, Beta(params.l2.a_rho, params.l2.b_rho))

"""
    apply_jump!(xi, rho, params) -> Float64

ジャンプ写像 Γ を発生直前の状態 `xi` に適用し、解放エネルギー
m = (1 - rho) * sigma_s^- を返す(SPEC §6.2/§6.3。全て加算)。
"""
function apply_jump!(xi::AbstractVector, rho::Real,
                     params::ModelParameters)
    l1, l2 = params.l1, params.l2
    sigma_s_minus = exp(xi[IX_SIG])
    m = (1 - rho) * sigma_s_minus
    v = exp(xi[IX_V])

    xi[IX_SIG] += log(rho)                              # 除荷
    xi[IX_TAU] -= l2.c_tau * m                          # 弾性衝撃
    xi[IX_TAUA] -= l2.c_star * m                        # 塑性変形
    xi[IX_G] -= l2.c_g * m                              # 格差圧縮
    xi[IX_K] -= l2.c_k * m                              # 資本破壊
    xi[IX_LAME] += l2.alpha0 * (v / l1.v0)^l2.kappa_alphav  # 自己励起
    return m
end

"1回のジャンプの記録(時刻・残留率・解放エネルギー)"
struct JumpEvent
    t::Float64
    rho::Float64
    m::Float64
end

"""
    simulate_hawkes(lam_b_const, params; t1, rng) -> Vector{Float64}

lam_b を定数に凍結した純 Hawkes 過程(状態力学なし)のイベント時刻列を
Ogata thinning(§10)で生成する。hawkes_stat テスト用
(平均イベント率の理論値は lam_b / (1 - n)、n = alpha0/beta)。
自己励起の増分は v = v0 のときの alpha0(§6.3)。
"""
function simulate_hawkes(lam_b_const::Real, params::ModelParameters;
                         t1::Real, rng::AbstractRNG)
    l1, l2 = params.l1, params.l2
    events = Float64[]
    t = 0.0
    lam_e = 0.0
    t_last = 0.0
    while true
        t += randexp(rng) / l2.lam_bar
        t >= t1 && break
        lam_e *= exp(-l1.beta * (t - t_last))   # 式(13)の厳密解で減衰
        t_last = t
        lam = min(l2.lam_bar, lam_b_const + lam_e)
        if rand(rng) < lam / l2.lam_bar
            push!(events, t)
            lam_e += l2.alpha0
        end
    end
    return events
end
