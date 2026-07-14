# Phase A 積分器(SPEC §10): 固定刻み Euler–Maruyama + Ogata thinning
#
# simulate_ode: ドリフトのみの決定論 ODE(M1)
# simulate_sde: 拡散 + ジャンプ(内生 thinning / 外生強制発火の2モード、M2)
# 格子点間の強度評価は直前状態で近似する(誤差 O(dt)、DECISIONS #0004/#0007)。

"シミュレーション結果(保持座標の軌道)"
struct Trajectory
    t::Vector{Float64}            # 時刻(年)
    X::Matrix{Float64}            # N_STATE × length(t)。保持座標
end

"""
    simulate_ode(params; t0=0.0, t1=50.0, dt=0.01, xi0=params.x0) -> Trajectory

ドリフトのみ(拡散・ジャンプなし)の決定論 ODE を前進 Euler
(dt = 0.01 年、§10 Phase A)で積分する。
"""
function simulate_ode(params::ModelParameters;
                      t0::Float64 = 0.0, t1::Float64 = 50.0, dt::Float64 = 0.01,
                      xi0::AbstractVector{Float64} = params.x0)
    nsteps = round(Int, (t1 - t0) / dt)
    ts = collect(range(t0; step = dt, length = nsteps + 1))
    X = Matrix{Float64}(undef, N_STATE, nsteps + 1)

    xi = copy(xi0)
    f = similar(xi)
    X[:, 1] = xi
    for step in 1:nsteps
        drift!(f, xi, params, ts[step])
        @. xi += dt * f
        X[:, step + 1] = xi
    end
    return Trajectory(ts, X)
end

"SDE + ジャンプのシミュレーション結果"
struct SDEResult
    traj::Trajectory
    jumps::Vector{JumpEvent}
end

# σ_s ステップガード(DECISIONS #0032): 式(11)の log 座標 drift は
# L·e^{-ξ_sig} 項をもち、L < 0 で σ_s → 0 に緩和する領域(実データの
# 安定国で恒常的に発生)では |drift| が指数発散して固定刻み EM が
# 破綻する。lam_bar・#0011 と同種の数値安全弁として、
#   (i) tame: |f_sig| ≤ SIG_DRIFT_MAX(双子実験の領域では |f_sig| ≲ 5
#       で不発火 — 力学・既存結果を変えない)
#   (ii) floor: ξ_sig ≥ XI_SIG_FLOOR(σ_s ≈ 6e-6。復帰時間を有界に保つ)
const SIG_DRIFT_MAX = 50.0    # /年(dt = 0.01 で最大 |Δξ| = 0.5)
const XI_SIG_FLOOR = -12.0

@inline guard_sigma_drift!(f::AbstractVector) =
    (f[IX_SIG] = clamp(f[IX_SIG], -SIG_DRIFT_MAX, SIG_DRIFT_MAX); f)
@inline guard_sigma_state!(xi::AbstractVector) =
    (xi[IX_SIG] = max(xi[IX_SIG], XI_SIG_FLOOR); xi)

# 格子区間 [t, t+dt) の内生ジャンプ(Ogata thinning、§10 擬似コード)。
# 受理のたびに Γ を適用し、以後の候補は更新後の状態で評価する。
#
# `gamma_thinning_p`(DECISIONS #0068、既定 1.0 = 現行動作): 予報窓の内生
# ジャンプに対する超過確率シンニング。候補時刻列の生成(lam_bar 上限レート)
# と候補判定用の乱数消費(rand(rng) の呼び出し回数・順序)は変えず、採択
# 判定 `u < intensity/lam_bar` の右辺に p を乗じることで採択確率のみを
# 下げる(u < p·intensity/lam_bar)。p=1 のとき既存条件と bitwise 同一。
function _jumps_in_interval!(xi::AbstractVector{Float64}, t::Float64, t_next::Float64,
                             params::ModelParameters, rng::AbstractRNG,
                             jumps::Vector{JumpEvent};
                             gamma_thinning_p::Float64 = 1.0)
    lam_bar = params.l2.lam_bar
    tj = t
    while true
        tj += randexp(rng) / lam_bar
        tj >= t_next && break
        if rand(rng) < gamma_thinning_p * intensity(xi, params) / lam_bar   # 直前状態で近似(O(dt))
            rho = draw_mark(rng, params)
            m = apply_jump!(xi, rho, params)
            push!(jumps, JumpEvent(tj, rho, m))
        end
    end
    return nothing
end

"""
    simulate_sde(params; seed, t0=0.0, t1=50.0, dt=0.01,
                 mode=EndogenousHawkes(), xi0=params.x0) -> SDEResult

拡散 + ジャンプ込みの Phase A 積分(固定刻み Euler–Maruyama + thinning、§10)。
`mode` により内生発火(EndogenousHawkes)と観測イベント時刻での強制発火
(ExogenousEvents、同化用。マーク rho は独立に引く)を切り替える(§6.4)。
乱数シードは明示的に受け取り、同一シードで完全再現される。
"""
function simulate_sde(params::ModelParameters;
                      seed::Integer,
                      t0::Float64 = 0.0, t1::Float64 = 50.0, dt::Float64 = 0.01,
                      mode::JumpMode = EndogenousHawkes(),
                      xi0::AbstractVector{Float64} = params.x0)
    rng = Xoshiro(seed)
    nsteps = round(Int, (t1 - t0) / dt)
    ts = collect(range(t0; step = dt, length = nsteps + 1))
    X = Matrix{Float64}(undef, N_STATE, nsteps + 1)
    jumps = JumpEvent[]

    exo_times = mode isa ExogenousEvents ? sort(mode.times) : Float64[]
    next_exo = 1

    xi = copy(xi0)
    f = similar(xi)
    sig = similar(xi)
    dW = similar(xi)
    X[:, 1] = xi
    sqdt = sqrt(dt)

    for step in 1:nsteps
        t = ts[step]
        t_next = ts[step + 1]

        # (a) 区間 [t, t+dt) のジャンプ処理(積分器から分離した関数、§6.4)
        if mode isa EndogenousHawkes
            _jumps_in_interval!(xi, t, t_next, params, rng, jumps)
        else
            while next_exo <= length(exo_times) && exo_times[next_exo] < t_next
                rho = draw_mark(rng, params)
                m = apply_jump!(xi, rho, params)
                push!(jumps, JumpEvent(exo_times[next_exo], rho, m))
                next_exo += 1
            end
        end

        # (b) Euler–Maruyama ステップ(σ_s ガード #0032)
        drift!(f, xi, params, t)
        guard_sigma_drift!(f)
        diffusion!(sig, xi, params, t)
        randn!(rng, dW)
        @. xi += dt * f + sqdt * sig * dW
        guard_sigma_state!(xi)
        X[:, step + 1] = xi
    end
    return SDEResult(Trajectory(ts, X), jumps)
end

"アンサンブル実行の結果(X は N_STATE × 時刻 × メンバー)"
struct EnsembleResult
    t::Vector{Float64}
    X::Array{Float64,3}
    jumps::Vector{Vector{JumpEvent}}
end

"メンバー別シードの導出(決定的。§10「メンバーごとに乱数シード固定」)"
member_seed(seed::Integer, i::Integer) = seed * 1_000_003 + i

"""
    simulate_ensemble(params; N=100, seed, kwargs...) -> EnsembleResult

N メンバーのアンサンブルを Threads 並列で実行する(メンバー間は独立、§10)。
各メンバーの乱数は member_seed(seed, i) で決定的に固定される。
`kwargs` は simulate_sde にそのまま渡す(t1, dt, mode, xi0 等)。
"""
function simulate_ensemble(params::ModelParameters;
                           N::Integer = 100, seed::Integer, kwargs...)
    results = Vector{SDEResult}(undef, N)
    Threads.@threads for i in 1:N
        results[i] = simulate_sde(params; seed = member_seed(seed, i), kwargs...)
    end
    nt = length(results[1].traj.t)
    X = Array{Float64,3}(undef, N_STATE, nt, N)
    for i in 1:N
        X[:, :, i] = results[i].traj.X
    end
    return EnsembleResult(results[1].traj.t, X, [r.jumps for r in results])
end
