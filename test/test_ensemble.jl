# アンサンブル実行のテスト(§10: メンバー独立・シード固定・Threads 並列)

@testset "ensemble" begin
    params = build_params(:volatile)
    ens = simulate_ensemble(params; N = 8, seed = 5, t1 = 10.0)

    @test size(ens.X) == (N_STATE, 1001, 8)
    @test length(ens.jumps) == 8
    @test all(isfinite, ens.X)

    @testset "members differ" begin
        @test ens.X[:, :, 1] != ens.X[:, :, 2]
    end

    @testset "reproducible and thread-count independent" begin
        ens2 = simulate_ensemble(params; N = 8, seed = 5, t1 = 10.0)
        @test ens.X == ens2.X
    end

    @testset "member i equals standalone run with member_seed" begin
        r3 = simulate_sde(params; seed = member_seed(5, 3), t1 = 10.0)
        @test ens.X[:, :, 3] == r3.traj.X
    end
end
