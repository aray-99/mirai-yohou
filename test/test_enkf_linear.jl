# enkf_linear テスト(SPEC §11):
# 2次元線形ガウス玩具問題で EnKF(N=10^4)≈ 厳密カルマンフィルタ
# (平均・分散の一致)

@testset "enkf_linear" begin
    rng = Xoshiro(123)

    # 玩具モデル: x_{k+1} = M x_k + q,  q ~ N(0, Q);  y = H x + r,  r ~ N(0, R)
    M = [0.9 0.1; -0.05 0.8]
    Q = [0.01 0.0; 0.0 0.01]
    H = [1.0 0.0]
    R = fill(0.04, 1, 1)
    m0 = [1.0, -0.5]
    P0 = [0.25 0.0; 0.0 0.25]

    nsteps = 20
    Nens = 10_000

    # 真値と観測
    x_true = copy(m0)
    ys = Vector{Float64}[]
    truth_rng = Xoshiro(7)
    for _ in 1:nsteps
        x_true = M * x_true + cholesky(Q).L * randn(truth_rng, 2)
        push!(ys, H * x_true + cholesky(R).L * randn(truth_rng, 1))
    end

    # 厳密カルマンフィルタ
    mkf, Pkf = copy(m0), copy(P0)
    # EnKF アンサンブル(インフレーションなしで厳密解と比較)
    E = m0 .+ cholesky(P0).L * randn(rng, 2, Nens)

    Qchol = cholesky(Q).L
    for k in 1:nsteps
        # 予報
        mkf = M * mkf
        Pkf = M * Pkf * M' + Q
        E .= M * E + Qchol * randn(rng, 2, Nens)
        # 解析(厳密 KF)
        S = H * Pkf * H' + R
        Kk = (Pkf * H') / S
        mkf = mkf + vec(Kk * (ys[k] - H * mkf))
        Pkf = (I - Kk * H) * Pkf
        # 解析(EnKF、rho_inf = 1 で純粋比較)
        enkf_analysis!(E, ys[k], xi -> H * xi, R; rng, rho_inf = 1.0)
    end

    m_ens = vec(sum(E; dims = 2)) ./ Nens
    Xp = E .- m_ens
    P_ens = (Xp * Xp') ./ (Nens - 1)

    # 平均: 事後標準偏差に対して十分近い(N=10^4 のサンプリング誤差 ~1/100)
    for i in 1:2
        @test abs(m_ens[i] - mkf[i]) < 5 * sqrt(Pkf[i, i]) / sqrt(Nens) * 10
    end
    # 分散・共分散: 相対誤差 10% 以内
    for i in 1:2, j in 1:2
        @test abs(P_ens[i, j] - Pkf[i, j]) < 0.1 * sqrt(Pkf[i, i] * Pkf[j, j])
    end

    @testset "analysis reduces variance toward observation" begin
        rng2 = Xoshiro(9)
        E2 = randn(rng2, 2, 2000) .* 2.0
        var_before = sum(abs2, E2[1, :] .- sum(E2[1, :]) / 2000) / 1999
        enkf_analysis!(E2, [0.0], xi -> [xi[1]], fill(0.01, 1, 1); rng = rng2,
                       rho_inf = 1.0)
        var_after = sum(abs2, E2[1, :] .- sum(E2[1, :]) / 2000) / 1999
        @test var_after < var_before
    end

    @testset "postprocess clamps lam_e" begin
        Emod = repeat(build_params(:stable).x0, 1, 4)
        Emod[IX_LAME, :] .= [-0.5, -1e-9, 0.0, 0.3]
        postprocess_analysis!(Emod)
        @test all(Emod[IX_LAME, :] .>= 0)
        @test Emod[IX_LAME, 4] ≈ 0.3
    end
end
