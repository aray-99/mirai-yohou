# coordinates テスト(SPEC §11): to_state ∘ from_state = 恒等(round-trip、境界近傍含む)

@testset "coordinates" begin
    @testset "basic transforms" begin
        @test logit(0.5) ≈ 0.0
        @test sigmoid(0.0) ≈ 0.5
        @test sigmoid(logit(0.3)) ≈ 0.3
        @test softplus(0.0) ≈ log(2)
        @test softplus(-100.0) > 0          # 常に正
        @test softplus(100.0) ≈ 100.0       # 大きな z で z に漸近(オーバーフローなし)
        @test softplus(-750.0) ≥ 0.0        # 極端な負値でも有限・非負
        @test isfinite(softplus(750.0))
        @test pluspart(-1.5) == 0.0
        @test pluspart(2.5) == 2.5
    end

    @testset "per-variable round-trip" begin
        # log 座標変数(正値)
        for i in (IX_P, IX_H, IX_K, IX_T, IX_V, IX_SIG), x in (1e-8, 0.3, 1.0, 1e6)
            @test from_state_var(i, to_state_var(i, x)) ≈ x rtol = 1e-12
        end
        # logit 座標変数(0〜1)、境界近傍を含む
        for i in (IX_W, IX_G, IX_PHI, IX_TAU, IX_PP),
            x in (1e-12, 1e-6, 0.01, 0.5, 0.99, 1 - 1e-6, 1 - 1e-12)
            @test from_state_var(i, to_state_var(i, x)) ≈ x rtol = 1e-6
        end
        # 恒等変数: xi_tauA(非有界)と lam_e(≥0)
        for x in (-5.0, 0.0, 3.7)
            @test to_state_var(IX_TAUA, x) == x
            @test from_state_var(IX_TAUA, x) == x
        end
        for x in (0.0, 0.1, 42.0)
            @test to_state_var(IX_LAME, x) == x
            @test from_state_var(IX_LAME, x) == x
        end
        @test_throws ArgumentError to_state_var(14, 1.0)
        @test_throws ArgumentError from_state_var(0, 1.0)
    end

    @testset "13-variable bulk round-trip" begin
        # 双子実験の初期条件(§8.4)
        for regime in (:stable, :volatile)
            params = build_params(regime)
            @test length(params.x0) == N_STATE
            @test from_state(params.x0) ≈ params.x0_nat rtol = 1e-10
            @test to_state(from_state(params.x0)) ≈ params.x0 atol = 1e-10
        end
        # 境界近傍を含む合成ベクトル
        x_nat = [0.5, 1e-9, 2.0, 0.7, 1 - 1e-9, 3.0, 0.5, 1.2,
                 1e-6, -2.5, 0.01, 1 - 1e-6, 0.3]
        @test from_state(to_state(x_nat)) ≈ x_nat rtol = 1e-6
        # ξ 側の往復は logit の飽和を避けた中庸値で検証
        xi = to_state([0.5, 0.6, 2.0, 0.7, 0.33, 3.0, 0.5, 1.2,
                       0.55, -2.5, 0.3, 0.3, 0.3])
        @test to_state(from_state(xi)) ≈ xi atol = 1e-10
        @test_throws DimensionMismatch to_state(ones(12))
        @test_throws DimensionMismatch from_state(ones(14))
    end
end
