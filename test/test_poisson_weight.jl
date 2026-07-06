# poisson_weight テスト(SPEC §11): 既知強度で重み・ESS・リサンプリングの正しさ

@testset "poisson_weight" begin
    @testset "weights match analytic Poisson likelihood ratio" begin
        Lambdas = [1.0, 2.0, 4.0]
        N_obs = 3
        w = normalize_weights(poisson_logweights(N_obs, Lambdas))
        # w_i ∝ Λ_i^N exp(-Λ_i)
        raw = [L^N_obs * exp(-L) for L in Lambdas]
        @test w ≈ raw ./ sum(raw) rtol = 1e-12
        @test sum(w) ≈ 1.0
    end

    @testset "ESS" begin
        @test ess(fill(0.25, 4)) ≈ 4.0          # 均等重み → N
        @test ess([1.0, 0.0, 0.0, 0.0]) ≈ 1.0   # 退化 → 1
        w = normalize_weights(poisson_logweights(5, [1.0, 5.0]))
        @test 1.0 <= ess(w) <= 2.0
    end

    @testset "systematic resampling" begin
        rng = Xoshiro(11)
        # 均等重みでは各メンバーがちょうど1回ずつ残る
        @test sort(systematic_resample(rng, fill(0.25, 4))) == [1, 2, 3, 4]
        @test sort(systematic_resample(rng, fill(1e-3, 1000))) == collect(1:1000)
        # 退化重みでは全複製
        @test all(systematic_resample(rng, [0.0, 1.0, 0.0]) .== 2)
        # 複製数は floor(N*w_i) と ceil(N*w_i) の間(系統リサンプリングの性質)
        w = [0.5, 0.3, 0.2]
        idx = systematic_resample(rng, w)
        @test length(idx) == 3
        for (i, wi) in enumerate(w)
            @test floor(3wi) <= count(==(i), idx) <= ceil(3wi)
        end
    end

    @testset "resample_if_needed!" begin
        rng = Xoshiro(21)
        # 均等な Λ → ESS = N → リサンプリングなし
        E = randn(rng, 2, 8)
        E0 = copy(E)
        w, essval, resampled = resample_if_needed!(E, fill(2.0, 8), 2; rng)
        @test !resampled
        @test E == E0
        @test essval ≈ 8.0
        # 極端な Λ 差 → ESS < N/2 → リサンプリングされ、支配メンバーが複製される
        E2 = Float64.(reshape(1:16, 2, 8))
        Lambdas = [1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 1e-6, 5.0]
        w2, ess2, res2 = resample_if_needed!(E2, Lambdas, 5; rng)
        @test res2
        @test ess2 < 4.0
        @test all(E2[:, j] == [15.0, 16.0] for j in 1:8)   # メンバー8が全複製
        @test argmax(w2) == 8
    end
end
