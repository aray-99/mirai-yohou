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

# 格子区間 [t, t+dt) の内生ジャンプ(Ogata thinning、§10 擬似コード)。
# 受理のたびに Γ を適用し、以後の候補は更新後の状態で評価する。
function _jumps_in_interval!(xi::Vector{Float64}, t::Float64, t_next::Float64,
                             params::ModelParameters, rng::AbstractRNG,
                             jumps::Vector{JumpEvent})
    lam_bar = params.l2.lam_bar
    tj = t
    while true
        tj += randexp(rng) / lam_bar
        tj >= t_next && break
        if rand(rng) < intensity(xi, params) / lam_bar   # 直前状態で近似(O(dt))
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

        # (b) Euler–Maruyama ステップ
        drift!(f, xi, params, t)
        diffusion!(sig, xi, params, t)
        randn!(rng, dW)
        @. xi += dt * f + sqdt * sig * dW
        X[:, step + 1] = xi
    end
    return SDEResult(Trajectory(ts, X), jumps)
end
