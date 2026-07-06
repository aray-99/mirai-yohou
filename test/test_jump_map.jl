# jump_map テスト(SPEC §11): Γ 適用後に
# sigma_s が rho 倍 / tau, tauA, g, k 減少 / lam_e 増加 / 他7変数不変

@testset "jump_map" begin
    params = build_params(:volatile)
    xi0 = copy(params.x0)

    rho = 0.6
    xi = copy(xi0)
    m = apply_jump!(xi, rho, params)

    sigma_s_before = exp(xi0[IX_SIG])

    @testset "mark energy" begin
        @test m ≈ (1 - rho) * sigma_s_before   # m = (1-rho) * sigma_s^-(§6.2)
    end

    @testset "unloading: sigma_s -> rho * sigma_s" begin
        @test exp(xi[IX_SIG]) ≈ rho * sigma_s_before
        @test xi[IX_SIG] ≈ xi0[IX_SIG] + log(rho)
    end

    @testset "impact: tau, tauA, g, k decrease" begin
        @test xi[IX_TAU] ≈ xi0[IX_TAU] - params.l2.c_tau * m
        @test xi[IX_TAUA] ≈ xi0[IX_TAUA] - params.l2.c_star * m
        @test xi[IX_G] ≈ xi0[IX_G] - params.l2.c_g * m
        @test xi[IX_K] ≈ xi0[IX_K] - params.l2.c_k * m
        @test xi[IX_TAU] < xi0[IX_TAU]
        @test xi[IX_TAUA] < xi0[IX_TAUA]
        @test xi[IX_G] < xi0[IX_G]
        @test xi[IX_K] < xi0[IX_K]
    end

    @testset "self-excitation: lam_e increases" begin
        v = exp(xi0[IX_V])
        expected = params.l2.alpha0 * (v / params.l1.v0)^params.l2.kappa_alphav
        @test xi[IX_LAME] ≈ xi0[IX_LAME] + expected
        @test xi[IX_LAME] > xi0[IX_LAME]
    end

    @testset "other 7 variables unchanged" begin
        for i in (IX_P, IX_W, IX_H, IX_T, IX_PHI, IX_V, IX_PP)
            @test xi[i] == xi0[i]
        end
    end

    @testset "rho = 1 releases nothing" begin
        xi1 = copy(xi0)
        m1 = apply_jump!(xi1, 1.0, params)
        @test m1 ≈ 0.0 atol = 1e-15
        @test xi1[IX_SIG] ≈ xi0[IX_SIG]        # log(1) = 0
        @test xi1[IX_TAU] ≈ xi0[IX_TAU]
        @test xi1[IX_LAME] > xi0[IX_LAME]      # 自己励起はマークによらず入る
    end
end
