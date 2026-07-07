# Phase B(SPEC §10 / DECISIONS #0020): SciML スタックによる純予測経路
#
# - simulate_sde_phaseb: JumpProcesses.jl の VariableRateJump + 適応 SDE
#   ソルバ(SOSRI)による内生 Hawkes ジャンプ拡散。Phase A(固定刻み EM +
#   thinning)との統計的一致は phaseb_agreement テストで検証する。
# - drift_jacobian_sparsity: Symbolics.jl によるドリフトヤコビアンの
#   疎性パターン自動検出(§5.1 の見積り ≈28% と照合。手書き禁止)。
#
# 同化(ExogenousEvents + EnKF)は Phase A の固定刻み経路を維持する(#0020)。

"""
    simulate_sde_phaseb(params; seed, t0=0.0, t1=50.0, saveat=0.1) -> SDEResult

Phase B 積分: 適応 SDE ソルバ(SOSRI)+ VariableRateJump による
内生 Hawkes ジャンプ(§10 Phase B)。保存点は `saveat` 刻み。
マーク・ジャンプ写像は Phase A と同一の draw_mark / apply_jump! を用いる。

!!! warning
    Threads.@threads によるアンサンブル並列化ではジャンプが発火しない
    (JumpProcesses の VariableRateJump がスレッド安全でない、#0021)。
    アンサンブルは逐次実行すること(1ラン ≈ 13ms で十分速い)。
"""
function simulate_sde_phaseb(params::ModelParameters;
                             seed::Integer,
                             t0::Float64 = 0.0, t1::Float64 = 50.0,
                             saveat::Float64 = 0.1)
    rng = Xoshiro(seed)
    fdrift = (du, u, p, t) -> (drift!(du, u, params, t); nothing)
    gdiff = (du, u, p, t) -> (diffusion!(du, u, params, t); nothing)
    prob = SDEProblem(fdrift, gdiff, copy(params.x0), (t0, t1))

    jumps = JumpEvent[]
    rate(u, p, t) = intensity(u, params)
    function affect!(integrator)
        rho = draw_mark(rng, params)
        m = apply_jump!(integrator.u, rho, params)
        push!(jumps, JumpEvent(integrator.t, rho, m))
        return nothing
    end
    jump = JumpProcesses.VariableRateJump(rate, affect!)
    jprob = JumpProcesses.JumpProblem(prob, jump; rng)

    sol = SciMLBase.solve(jprob, StochasticDiffEq.SOSRI(); saveat)
    ts = collect(sol.t)
    X = Matrix{Float64}(undef, N_STATE, length(ts))
    for (k, u) in enumerate(sol.u)
        for i in 1:N_STATE
            X[i, k] = u[i]
        end
    end
    return SDEResult(Trajectory(ts, X), jumps)
end

"""
    drift_jacobian_sparsity(params) -> SparseMatrixCSC{Bool}

ドリフト場のヤコビアン疎性パターンを Symbolics.jacobian_sparsity で
自動検出する(§5.1/§10。手書きしない)。
"""
function drift_jacobian_sparsity(params::ModelParameters)
    f! = (du, u) -> drift!(du, u, params, 0.0)
    return Symbolics.jacobian_sparsity(f!, zeros(N_STATE), zeros(N_STATE))
end
