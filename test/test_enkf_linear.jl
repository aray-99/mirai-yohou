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

    @testset "rtps" begin
        rng3 = Xoshiro(31)
        Eprior = randn(rng3, 3, 5000) .* [2.0, 1.0, 0.5]
        spread_prior = ensemble_spread(Eprior)
        Epost = copy(Eprior)
        # 成分1のみ観測で解析(成分3は無相関 → スプレッド不変のはず)
        enkf_analysis!(Epost, [0.0], xi -> [xi[1]], fill(0.01, 1, 1);
                       rng = rng3, rho_inf = 1.0)
        m_before = vec(sum(Epost; dims = 2)) ./ 5000
        Ertps = copy(Epost)
        rtps!(Ertps, spread_prior; alpha = 0.5)

        # (i) 平均を変えない
        @test vec(sum(Ertps; dims = 2)) ./ 5000 ≈ m_before atol = 1e-12
        # (ii) 緩和後スプレッド = σ_post + α(σ_prior − σ_post)
        s_post = ensemble_spread(Epost)
        s_rtps = ensemble_spread(Ertps)
        @test s_rtps ≈ s_post .+ 0.5 .* (spread_prior .- s_post) rtol = 1e-10
        # 観測成分は事前と事後の間、非観測成分はほぼ不変
        @test s_post[1] < s_rtps[1] < spread_prior[1]
        @test isapprox(s_rtps[3], s_post[3]; rtol = 0.05)
        # (iii) α = 0 は無効化
        E0rt = copy(Epost)
        rtps!(E0rt, spread_prior; alpha = 0.0)
        @test E0rt == Epost
    end

    @testset "enks fixed-lag smoother" begin
        # 2次元 AR(1): x_{k+1} = 0.9 x_k + q。時刻1のスナップショットを保持し、
        # 時刻2の観測で平滑化(#0024)
        rng4 = Xoshiro(41)
        Nens = 5000
        S1 = randn(rng4, 2, Nens)                       # 時刻1の状態
        E2 = 0.9 .* S1 .+ 0.3 .* randn(rng4, 2, Nens)   # 時刻2へ伝播
        S1b = copy(S1)
        yobs = [0.5]
        R1 = fill(0.04, 1, 1)

        # (i) スナップショットなしの enks == enkf(同一 rng シード)
        Ea = copy(E2); Eb = copy(E2)
        enkf_analysis!(Ea, yobs, xi -> [xi[1]], R1; rng = Xoshiro(5), rho_inf = 1.0)
        enks_analysis!(Eb, Matrix{Float64}[], yobs, xi -> [xi[1]], R1;
                       rng = Xoshiro(5), rho_inf = 1.0)
        @test Ea == Eb

        # (ii) 平滑化は過去状態の相関成分の分散を減らし、平均を観測方向へ更新
        var1_before = ensemble_spread(S1b) .^ 2
        enks_analysis!(E2, [S1b], yobs, xi -> [xi[1]], R1;
                       rng = Xoshiro(5), rho_inf = 1.0)
        var1_after = ensemble_spread(S1b) .^ 2
        @test var1_after[1] < var1_before[1]            # 観測と相関する成分
        @test isapprox(var1_after[2], var1_before[2]; rtol = 0.05)  # 無相関成分は不変
        m1 = vec(sum(S1b; dims = 2)) ./ Nens
        @test 0 < m1[1] < 0.5                            # 事前平均0から観測0.5方向へ
    end

    @testset "enks analysis masked_rows (#0040-(α))" begin
        # 2次元 AR(1) 玩具問題(強く相関: row2 は row1 の観測に強く反応する
        # ように相関を作る)。masked_rows で行2の現在時刻更新を無効化できるか。
        rng5 = Xoshiro(51)
        Nens = 5000
        S0 = randn(rng5, 2, Nens)
        E = copy(S0)
        E[2, :] .= 0.9 .* E[1, :] .+ 0.1 .* randn(rng5, Nens)   # row2 は row1 と強相関
        yobs = [0.5]
        R1 = fill(0.04, 1, 1)

        # (i) masked_rows 既定(空)= 従来動作: 行2も更新される
        Edefault = copy(E)
        enks_analysis!(Edefault, Matrix{Float64}[], yobs, xi -> [xi[1]], R1;
                       rng = Xoshiro(9), rho_inf = 1.0)
        rowmean(A, i) = sum(@view A[i, :]) / size(A, 2)
        mbar2_default = rowmean(Edefault, 2)
        @test mbar2_default > rowmean(E, 2) + 0.05   # 観測方向へ有意に動く

        # (ii) masked_rows = [2]: 行2の現在時刻更新はゼロ化され平均不変
        Emasked = copy(E)
        enks_analysis!(Emasked, Matrix{Float64}[], yobs, xi -> [xi[1]], R1;
                       rng = Xoshiro(9), rho_inf = 1.0, masked_rows = [2])
        @test rowmean(Emasked, 2) ≈ rowmean(E, 2) atol = 1e-9
        # 行1(マスクなし)は従来どおり更新される
        @test rowmean(Emasked, 1) ≈ rowmean(Edefault, 1) atol = 1e-9

        # (iii) 過去平滑化(snapshots)側は masked_rows の影響を受けない
        Ssnap_a = copy(S0); Ssnap_b = copy(S0)
        Ea = copy(E); Eb = copy(E)
        enks_analysis!(Ea, [Ssnap_a], yobs, xi -> [xi[1]], R1;
                       rng = Xoshiro(9), rho_inf = 1.0)
        enks_analysis!(Eb, [Ssnap_b], yobs, xi -> [xi[1]], R1;
                       rng = Xoshiro(9), rho_inf = 1.0, masked_rows = [2])
        @test Ssnap_a ≈ Ssnap_b atol = 1e-9
    end

    @testset "postprocess clamps lam_e" begin
        Emod = repeat(build_params(:stable).x0, 1, 4)
        Emod[IX_LAME, :] .= [-0.5, -1e-9, 0.0, 0.3]
        postprocess_analysis!(Emod)
        @test all(Emod[IX_LAME, :] .>= 0)
        @test Emod[IX_LAME, 4] ≈ 0.3
    end
end
