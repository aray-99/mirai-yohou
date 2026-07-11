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

    @testset "NegBin 尤度オプション (#0054)" begin
        # 既知パラメータでの log pmf(Γ 関数を使わない初等的同一性で検証):
        # P(0|mu,r) = (r/(r+mu))^r、P(N+1)/P(N) = (N+r)/(N+1) · mu/(r+mu)
        for (mu, r) in ((2.0, 5.0), (1.5, 0.7), (8.0, 2.0))
            @test negbin_logpmf(0, mu, r) ≈ r * log(r / (r + mu)) rtol = 1e-10
            for N_obs in 0:10
                ratio = log((N_obs + r) / (N_obs + 1) * mu / (r + mu))
                @test negbin_logpmf(N_obs + 1, mu, r) - negbin_logpmf(N_obs, mu, r) ≈
                      ratio atol = 1e-9
            end
            # 正規化: Σ_N P(N) = 1(裾は幾何減衰なので有限和で十分)
            @test sum(exp(negbin_logpmf(N_obs, mu, r)) for N_obs in 0:2000) ≈ 1.0 rtol = 1e-6
        end

        # r → ∞ でポアソンに収束
        @test negbin_logpmf(4, 3.0, 1e8) ≈ poisson_logpmf(4, 3.0) atol = 1e-5

        # poisson_logpmf は poisson_logweights + log N! 定数
        @test poisson_logpmf(4, 3.0) ≈
              poisson_logweights(4, [3.0])[1] - log(factorial(4)) rtol = 1e-12

        # 重みの正規化: 定数項はキャンセルし、メンバー間比だけが残る
        mus = [1.0, 2.0, 4.0]
        w = normalize_weights(negbin_logweights(3, mus, 5.0))
        raw = [exp(negbin_logpmf(3, mu, 5.0)) for mu in mus]
        @test w ≈ raw ./ sum(raw) rtol = 1e-10
        @test sum(w) ≈ 1.0

        # プロファイル MLE の回復: 既知 r から生成したカウントで r̂ ≈ r
        rng = Xoshiro(3054)
        r_true = 3.0
        mus_sim = 2.0 .+ 8.0 .* rand(rng, 400)
        counts = [rand(rng, MiraiYohou.NegativeBinomial(r_true, r_true / (r_true + mu)))
                  for mu in mus_sim]
        r_hat = negbin_profile_r(counts, mus_sim)
        @test 0.5 * r_true < r_hat < 2.0 * r_true

        # ポアソン生成データでは r̂ が大きく(過分散なし)、その NegBin は
        # 実質ポアソン(有限標本の揺らぎで r̂ は有限になりうるため、pmf の
        # 近さで判定する)
        counts_p = [rand(rng, MiraiYohou.Poisson(mu)) for mu in mus_sim]
        r_hat_p = negbin_profile_r(counts_p, mus_sim)
        @test r_hat_p > 10 * r_true
        @test abs(negbin_logpmf(5, 5.0, r_hat_p) - poisson_logpmf(5, 5.0)) < 0.05

        # 入力検証
        @test_throws DimensionMismatch negbin_profile_r([1, 2], [1.0])
        @test_throws ArgumentError negbin_profile_r(Int[], Float64[])
    end

    @testset "run_assimilation の count_model 切り替え (#0054)" begin
        params = build_params(:volatile)
        t1 = 2.0
        truth = simulate_sde(params; seed = 5401, t1 = t1)
        event_times = [e.t for e in truth.jumps]
        obs = synthesize_observations(truth.traj, standard_observations(params);
                                      rng = Xoshiro(5402))
        N = 30
        E0 = params.x0 .+ 0.3 .* randn(Xoshiro(5403), MiraiYohou.N_STATE, N)
        postprocess_analysis!(E0)
        cfg = AssimConfig(t1 = t1)
        nwin = round(Int, t1 / cfg.event_window)

        # 既定(:poisson)は従来動作 + 診断フィールドが窓数分埋まる
        r1 = run_assimilation(params, copy(E0), obs, event_times;
                              cfg, seed = 5404)
        @test length(r1.count_observed) == nwin
        @test length(r1.count_logscore) == nwin
        @test all(r1.count_observed .>= 0)          # カタログ集計は常にデータあり
        @test all(isfinite, r1.count_logscore)
        @test all(r1.count_logscore .<= 0)          # log pmf ≤ 0

        # :negbin(大 r)は :poisson と実質同一の重み → 軌道もほぼ同一
        # (完全一致は乱数消費順序が同じことに依存。ここでは同一シードで
        #  X が一致することを確認する — 重み計算は乱数を消費しない)
        r2 = run_assimilation(params, copy(E0), obs, event_times;
                              cfg, seed = 5404,
                              count_model = :negbin, count_dispersion = 1e10)
        @test r2.X ≈ r1.X rtol = 1e-6

        # 小さい r(強い過分散許容)では重みが均され、ESS が下がらない方向
        r3 = run_assimilation(params, copy(E0), obs, event_times;
                              cfg, seed = 5404,
                              count_model = :negbin, count_dispersion = 0.5)
        @test length(r3.ess) == nwin
        @test all(isfinite, r3.count_logscore)

        # 不正な指定はエラー
        @test_throws ArgumentError run_assimilation(params, copy(E0), obs,
                                                    event_times; cfg, seed = 5404,
                                                    count_model = :negbin)
        @test_throws ArgumentError run_assimilation(params, copy(E0), obs,
                                                    event_times; cfg, seed = 5404,
                                                    count_model = :bogus)
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
