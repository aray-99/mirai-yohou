# Phase B テスト(SPEC §10 / DECISIONS #0020。トレランスは #0020 で事前凍結)
#
# phaseb_agreement: Phase A(固定刻み EM + thinning)と Phase B(SOSRI +
# VariableRateJump)のアンサンブル統計の一致。
# sparsity: Symbolics によるドリフトヤコビアン疎性が §5.1 の見積りと整合。

@testset "phaseb_agreement" begin
    # N=400: 稀ジャンプレジームでは分散がジャンプ有無の重い裾に支配され、
    # N=100 では分散比の標本雑音が大きい(実測 3.5 → N=400 で 1.7 に収束)
    N = 400
    t1 = 30.0
    check_vars = [(:xi_g, IX_G), (:xi_tau, IX_TAU), (:xi_k, IX_K), (:xi_sig, IX_SIG)]

    for regime in (:stable, :volatile)
        params = build_params(regime)

        ensA = simulate_ensemble(params; N, seed = 33, t1)
        finalA = ensA.X[:, end, :]
        jumpsA = [length(j) for j in ensA.jumps]

        # 逐次実行(VariableRateJump はスレッド並列で発火しない、#0021)
        rsB = [simulate_sde_phaseb(params; seed = member_seed(77, i), t1)
               for i in 1:N]
        finalB = reduce(hcat, [r.traj.X[:, end] for r in rsB])
        jumpsB = [length(r.jumps) for r in rsB]

        @testset "$regime: final-time moments agree" begin
            for (name, ix) in check_vars
                a = @view finalA[ix, :]
                b = @view finalB[ix, :]
                ma, mb = sum(a) / N, sum(b) / N
                va = sum(abs2, a .- ma) / (N - 1)
                vb = sum(abs2, b .- mb) / (N - 1)
                se = sqrt(va / N + vb / N)
                @test abs(ma - mb) <= 4 * se        # 平均差 ≤ 4×結合SE(#0020)
                @test 0.5 <= va / vb <= 2.0          # 分散比 ∈ [0.5, 2](#0020)
            end
        end

        @testset "$regime: event counts agree" begin
            ma, mb = sum(jumpsA) / N, sum(jumpsB) / N
            pooled = (ma + mb) / 2
            @test abs(ma - mb) <= 4 * sqrt(max(pooled, 0.1) / N)  # #0020
        end

        @testset "$regime: phase B invariants" begin
            @test all(all(isfinite, r.traj.X) for r in rsB)
            @test all(all(diff(r.traj.X[IX_T, :]) .>= -1e-12) for r in rsB)
        end
    end
end

@testset "sparsity" begin
    params = build_params(:volatile)
    S = drift_jacobian_sparsity(params)
    @test size(S) == (N_STATE, N_STATE)
    density = count(!iszero, S) / (N_STATE * N_STATE)
    # §5.1 の実効密度見積り ≈ 28%(48/169)。検出は 0.20〜0.35 の帯で整合(#0020)
    @test 0.20 <= density <= 0.35
    # 式(11)の載荷則(yhat 経由でほぼ全域に依存)が最密の行
    row_nnz = [count(!iszero, S[i, :]) for i in 1:N_STATE]
    @test row_nnz[IX_SIG] == maximum(row_nnz)
    @test row_nnz[IX_SIG] >= 6
    # 外生のみの式(1)は自分自身にも依存しない
    @test row_nnz[IX_P] == 0
    # 式(13)は純減衰(lam_e 自身のみ)
    @test row_nnz[IX_LAME] == 1 && S[IX_LAME, IX_LAME]
end
