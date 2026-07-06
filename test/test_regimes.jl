# regimes テスト(SPEC §11): 各20回反復で判定
# 安定国セット: 50年でジャンプ ≈ 0(高確率)。変動国セット: 少なくとも1回発火(高確率)。
# 「高確率」の定量化は DECISIONS #0008(固定シードで決定的に実行)。

@testset "regimes" begin
    nreps = 20
    counts = Dict{Symbol,Vector{Int}}()
    for regime in (:stable, :volatile)
        params = build_params(regime)
        counts[regime] = [length(simulate_sde(params; seed = 1000 + r,
                                              t1 = 50.0).jumps)
                          for r in 1:nreps]
    end

    @testset "stable: jumps ≈ 0 with high probability" begin
        zero_runs = count(==(0), counts[:stable])
        @test zero_runs >= 14                              # 実測 16/20(#0008)
        @test sum(counts[:stable]) / nreps < 1.0           # 平均1回未満
    end

    @testset "volatile: at least one jump with high probability" begin
        fired = count(>=(1), counts[:volatile])
        @test fired >= 19                                  # 実測 20/20(#0008)
    end

    @testset "regime separation" begin
        @test sum(counts[:volatile]) > 10 * max(sum(counts[:stable]), 1)
    end
end
