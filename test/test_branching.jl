# branching テスト(SPEC §11):
# パラメータ構築時 n = alpha0/beta < 1 を assert(n≥1 でエラーになることも試験)

@testset "branching" begin
    @testset "valid parameter sets" begin
        for regime in (:stable, :volatile)
            params = build_params(regime)
            d = dimensionless_numbers(params)
            @test d.n ≈ 0.4          # alpha0/beta = 2.0/5.0(§8.2)
            @test d.n < 1
            @test d.H ≈ 0.3          # c_star/c_tau = 0.3/1.0
            @test d.theta_sig ≈ 3.0  # L3 事前中央値(§8.3)
        end
        # De のレジーム条件(§8.4、DECISIONS #0006)
        @test dimensionless_numbers(build_params(:stable)).De ≈ 0.9 atol = 1e-12
        @test dimensionless_numbers(build_params(:stable)).De < 1
        @test dimensionless_numbers(build_params(:volatile)).De ≈ 5.0 atol = 1e-12
        @test dimensionless_numbers(build_params(:volatile)).De > 1
    end

    @testset "n >= 1 is rejected" begin
        # alpha0 = beta = 5 → n = 1(ちょうど1も禁止)
        @test_throws ArgumentError build_params(:stable; alpha0 = 5.0)
        # n = 1.2
        @test_throws ArgumentError build_params(:volatile; alpha0 = 6.0)
    end

    @testset "De regime assertion" begin
        # 安定国で応力載荷を強めると De ≥ 1 になり拒否される
        @test_throws ArgumentError build_params(:stable; eta_p = 0.5)   # De = 1.5
        # 変動国で緩和を強めると De ≤ 1 になり拒否される
        @test_throws ArgumentError build_params(:volatile; delta_sig = 0.5)  # De = 0.5
    end

    @testset "misc construction" begin
        @test_throws ArgumentError build_params(:unknown)
        p = build_params(:stable)
        @test p.exo.netgrowth ≈ -0.003
        @test p.exo.wbar ≈ 0.58
        v = build_params(:volatile)
        @test v.exo.netgrowth ≈ +0.015
        @test v.exo.wbar ≈ 0.62
        @test v.l2.lam0 ≈ 0.10
        @test v.l2.eta_g ≈ 1.0
        @test sigmoid(v.l2.mu_gbar) ≈ 0.45
    end
end
