# invariants テスト(SPEC §11):
# 50年シミュレーションで T 単調非減少 / sigma_s 有限・正 / 全 logit 座標 |ξ|<10
#
# 決定論 ODE(M1)と拡散+ジャンプ込みの SDE(M2)の両方で検証する。

function check_invariants(traj)
    @testset "no divergence" begin
        @test all(isfinite, traj.X)
    end
    @testset "T monotone nondecreasing" begin
        @test all(diff(@view traj.X[IX_T, :]) .>= 0)
    end
    @testset "sigma_s finite and positive" begin
        sig = exp.(@view traj.X[IX_SIG, :])
        @test all(isfinite, sig)
        @test all(sig .> 0)
    end
    @testset "logit coordinates bounded" begin
        for i in (IX_W, IX_G, IX_PHI, IX_TAU, IX_PP)
            @test maximum(abs, @view traj.X[i, :]) < 10
        end
    end
end

@testset "invariants" begin
    @testset "SDE with jumps: $regime (seed $seed)" for regime in (:stable, :volatile),
                                                        seed in (1, 2)
        params = build_params(regime)
        res = simulate_sde(params; seed, t1 = 50.0)
        check_invariants(res.traj)
        @testset "lam_e nonnegative" begin
            @test all(res.traj.X[IX_LAME, :] .>= 0)
        end
    end

    @testset "drift-only ODE: $regime" for regime in (:stable, :volatile)
        params = build_params(regime)
        traj = simulate_ode(params; t1 = 50.0)
        check_invariants(traj)
        @testset "lam_e stays zero without jumps" begin
            @test all(traj.X[IX_LAME, :] .== 0)
        end
    end
end
