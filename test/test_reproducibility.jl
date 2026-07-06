# reproducibility テスト(SPEC §11): 同一シードで完全再現

@testset "reproducibility" begin
    for regime in (:stable, :volatile)
        params = build_params(regime)

        r1 = simulate_sde(params; seed = 42, t1 = 20.0)
        r2 = simulate_sde(params; seed = 42, t1 = 20.0)
        r3 = simulate_sde(params; seed = 43, t1 = 20.0)

        @testset "$regime: identical with same seed" begin
            @test r1.traj.X == r2.traj.X
            @test length(r1.jumps) == length(r2.jumps)
            @test all(e1.t == e2.t && e1.rho == e2.rho && e1.m == e2.m
                      for (e1, e2) in zip(r1.jumps, r2.jumps))
        end

        @testset "$regime: different seed diverges" begin
            @test r1.traj.X != r3.traj.X
        end
    end

    @testset "exogenous mode is reproducible" begin
        params = build_params(:volatile)
        mode = ExogenousEvents([3.0, 7.5, 12.25])
        r1 = simulate_sde(params; seed = 1, t1 = 20.0, mode)
        r2 = simulate_sde(params; seed = 1, t1 = 20.0, mode)
        @test r1.traj.X == r2.traj.X
        @test [e.t for e in r1.jumps] == [3.0, 7.5, 12.25]   # 強制発火時刻に一致
    end
end
