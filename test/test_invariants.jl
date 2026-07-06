# invariants テスト(SPEC §11):
# 50年シミュレーションで T 単調非減少 / sigma_s 有限・正 / 全 logit 座標 |ξ|<10
#
# M1 時点はドリフトのみの決定論 ODE で検証する。M2 で拡散+ジャンプに拡張する。

@testset "invariants" begin
    for regime in (:stable, :volatile)
        params = build_params(regime)
        traj = simulate_ode(params; t1 = 50.0)

        @testset "$regime: no divergence" begin
            @test all(isfinite, traj.X)
        end

        @testset "$regime: T monotone nondecreasing" begin
            xiT = @view traj.X[IX_T, :]
            @test all(diff(xiT) .>= 0)
        end

        @testset "$regime: sigma_s finite and positive" begin
            sig = exp.(@view traj.X[IX_SIG, :])
            @test all(isfinite, sig)
            @test all(sig .> 0)
        end

        @testset "$regime: logit coordinates bounded" begin
            for i in (IX_W, IX_G, IX_PHI, IX_TAU, IX_PP)
                @test maximum(abs, @view traj.X[i, :]) < 10
            end
        end

        @testset "$regime: lam_e stays zero without jumps" begin
            @test all(traj.X[IX_LAME, :] .== 0)
        end
    end
end
