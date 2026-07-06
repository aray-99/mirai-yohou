# hawkes_stat テスト(SPEC §11):
# lam_b を定数に凍結して長時間シミュレーション。平均イベント率が
# 理論値 lam_b / (1 - n) の信頼区間内(指数カーネル Hawkes の既知公式)。

@testset "hawkes_stat" begin
    params = build_params(:volatile)
    n = branching_ratio(params.l1, params.l2)   # 0.4
    lam_b_const = 0.5
    T = 3000.0

    rate_theory = lam_b_const / (1 - n)
    # 定常 Hawkes のカウント分散率 ≈ lam_bar_rate / (1-n)^2
    sd_count = sqrt(T * rate_theory / (1 - n)^2)

    events = simulate_hawkes(lam_b_const, params; t1 = T,
                             rng = Xoshiro(20260707))
    N = length(events)

    @test all(diff(events) .> 0)                       # 時刻は狭義単調増加
    @test abs(N - T * rate_theory) < 4 * sd_count      # 理論平均の 4σ 帯

    @testset "no self-excitation reduces to Poisson" begin
        p0 = build_params(:volatile; alpha0 = 0.0)     # n = 0
        ev0 = simulate_hawkes(lam_b_const, p0; t1 = T, rng = Xoshiro(7))
        @test abs(length(ev0) - T * lam_b_const) < 4 * sqrt(T * lam_b_const)
        # 自己励起ありは平均レートが厳密に増える(同一 T・十分長い系列)
        @test N > length(ev0)
    end
end
