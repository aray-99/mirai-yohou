# Phase A 積分器(SPEC §10): 固定刻み Euler–Maruyama + Ogata thinning
#
# M1 スコープ: ドリフトのみ(拡散・ジャンプなし)の決定論 ODE(前進 Euler)。
# 拡散(M2)・ジャンプ(M2)・アンサンブル(M3)は同じ骨格に追加する。

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
