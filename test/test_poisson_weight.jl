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

    @testset "Liu-West rejuvenation (#0057)" begin
        @testset "weighted_mean_std" begin
            E = [1.0 2.0 3.0; 10.0 20.0 30.0]
            w = [0.5, 0.25, 0.25]
            xbar, sigma = MiraiYohou.weighted_mean_std(E, w)
            @test xbar ≈ [0.5*1 + 0.25*2 + 0.25*3, 0.5*10 + 0.25*20 + 0.25*30]
            expected_var1 = 0.5*(1-xbar[1])^2 + 0.25*(2-xbar[1])^2 + 0.25*(3-xbar[1])^2
            @test sigma[1] ≈ sqrt(expected_var1)
            # 定数行(分散ゼロ)は sigma = 0
            Econst = [5.0 5.0 5.0]
            _, sigma_c = MiraiYohou.weighted_mean_std(Econst, [1/3, 1/3, 1/3])
            @test sigma_c[1] ≈ 0.0 atol = 1e-12
        end

        @testset "liu_west_rejuvenate!: a=1 は不変" begin
            rng = Xoshiro(41)
            E = randn(rng, 3, 20)
            E0 = copy(E)
            xbar, sigma = MiraiYohou.weighted_mean_std(E, fill(1/20, 20))
            MiraiYohou.liu_west_rejuvenate!(E, xbar, sigma, 1.0, rng)
            @test E ≈ E0   # b=√(1-1²)=0 → ジッタなし、xbar項も a=1 で無効
        end

        @testset "liu_west_rejuvenate!: モーメント近似保存(a<1)" begin
            rng = Xoshiro(42)
            n, N = 4, 20000
            E = randn(rng, n, N)   # 平均0・分散1のi.i.d.列(等重み想定)
            w = fill(1.0 / N, N)
            xbar, sigma = MiraiYohou.weighted_mean_std(E, w)
            a = 0.9
            MiraiYohou.liu_west_rejuvenate!(E, xbar, sigma, a, rng)
            _mean(v) = sum(v) / length(v)
            _std(v) = (m = _mean(v); sqrt(sum(abs2, v .- m) / (length(v) - 1)))
            # 1次モーメント: 平均は保存されるはず(理論上 xbar 自身に収束)
            for k in 1:n
                @test isapprox(_mean(E[k, :]), xbar[k]; atol = 0.05)
                # 2次モーメント: a²σ² + (1-a²)σ² = σ²(縮小核は分散を保存)
                @test isapprox(_std(E[k, :]), sigma[k]; rtol = 0.05)
            end
        end

        @testset "sigma=0 の座標はジッタしない" begin
            rng = Xoshiro(43)
            n, N = 2, 10
            E = zeros(n, N)
            E[1, :] .= 3.0                 # 定数行(sigma=0)
            E[2, :] .= randn(rng, N)
            xbar = [3.0, 0.0]
            sigma = [0.0, 1.0]
            a = 0.9
            E_before_row1 = copy(E[1, :])
            MiraiYohou.liu_west_rejuvenate!(E, xbar, sigma, a, rng)
            # row 1: a*3 + (1-a)*3 + 0 = 3(不変)
            @test E[1, :] ≈ E_before_row1
        end

        @testset "run_assimilation: rejuvenation_a=1.0 は従来動作とビット一致" begin
            params = build_params(:volatile)
            t1 = 2.0
            truth = simulate_sde(params; seed = 5501, t1 = t1)
            event_times = [e.t for e in truth.jumps]
            obs = synthesize_observations(truth.traj, standard_observations(params);
                                          rng = Xoshiro(5502))
            N = 30
            E0 = params.x0 .+ 0.3 .* randn(Xoshiro(5503), MiraiYohou.N_STATE, N)
            postprocess_analysis!(E0)
            cfg_legacy = AssimConfig(t1 = t1)
            cfg_explicit = AssimConfig(t1 = t1, rejuvenation_a = 1.0)
            r_legacy = run_assimilation(params, copy(E0), obs, event_times;
                                       cfg = cfg_legacy, seed = 5504)
            r_explicit = run_assimilation(params, copy(E0), obs, event_times;
                                         cfg = cfg_explicit, seed = 5504)
            @test r_legacy.X == r_explicit.X
            @test r_legacy.nresample == r_explicit.nresample
            @test r_legacy.ess == r_explicit.ess
        end

        @testset "run_assimilation: rejuvenation_a<1 は再抽選後にスプレッドを回復する" begin
            params = build_params(:volatile)
            t1 = 2.0
            truth = simulate_sde(params; seed = 5601, t1 = t1)
            event_times = [e.t for e in truth.jumps]
            obs = synthesize_observations(truth.traj, standard_observations(params);
                                          rng = Xoshiro(5602))
            N = 30
            E0 = params.x0 .+ 0.3 .* randn(Xoshiro(5603), MiraiYohou.N_STATE, N)
            postprocess_analysis!(E0)
            cfg_none = AssimConfig(t1 = t1)
            cfg_rej = AssimConfig(t1 = t1, rejuvenation_a = 0.9)
            r_none = run_assimilation(params, copy(E0), obs, event_times;
                                     cfg = cfg_none, seed = 5604)
            r_rej = run_assimilation(params, copy(E0), obs, event_times;
                                    cfg = cfg_rej, seed = 5604)
            @test r_none.nresample == r_rej.nresample   # 再抽選判定自体は不変
            if r_none.nresample > 0
                @test r_rej.X != r_none.X   # 若返りが軌道に反映されている
            end
        end
    end

    @testset "g 計装診断 instrument_g (#0065)" begin
        # 標準観測(twin)に :g_swiid 年次観測を1本だけ追加した簡易シナリオ。
        # 実データ層(build_observations.jl)を経由せず、ObservationSpec/
        # ObservationRecord を直接組み立てる(パッケージテストは experiments
        # 側のデータファイルに依存させない、test_augmentation.jl と同じ流儀)。
        params = build_params(:volatile)
        t1 = 3.0
        truth = simulate_sde(params; seed = 7101, t1 = t1)
        event_times = [e.t for e in truth.jumps]
        obs = synthesize_observations(truth.traj, standard_observations(params);
                                      rng = Xoshiro(7102))
        g_swiid_spec = ObservationSpec(:g_swiid, 1.0, 0.05, xi -> xi[IX_G], IX_G)
        g_swiid_obs = ObservationRecord(2.0, g_swiid_spec, 0.3)
        obs_all = sort(vcat(obs, [g_swiid_obs]); by = o -> o.t)
        N = 30
        E0 = params.x0 .+ 0.3 .* randn(Xoshiro(7103), N_STATE, N)
        postprocess_analysis!(E0)
        cfg = AssimConfig(t1 = t1)
        seed = 7104

        r_off = run_assimilation(params, copy(E0), obs_all, event_times; cfg, seed)
        r_on = run_assimilation(params, copy(E0), obs_all, event_times; cfg, seed,
                                instrument_g = true)

        # 判定数値へ影響ゼロ(絶対条件、#0065): 既定無効時と有効時で
        # g_diag 以外の全フィールドがビット一致する。
        @test r_off.X == r_on.X
        @test r_off.ranks == r_on.ranks
        @test r_off.ess == r_on.ess
        @test r_off.nresample == r_on.nresample
        @test r_off.ts_snap == r_on.ts_snap
        @test r_off.Xs == r_on.Xs
        @test r_off.count_observed == r_on.count_observed
        @test isequal(r_off.count_logscore, r_on.count_logscore)   # NaN 対応

        # 既定無効時は計装ログが空(計算自体を行わない)
        @test isempty(r_off.g_diag)

        # 有効時はログが記録され、週次カウント更新と観測解析の両方を含む
        @test !isempty(r_on.g_diag)
        types = Set(d["update_type"] for d in r_on.g_diag)
        @test "count_weekly" in types
        @test "g_swiid_annual" in types
        # pre/post が対で記録される
        @test count(d -> d["phase"] == "pre", r_on.g_diag) ==
              count(d -> d["phase"] == "post", r_on.g_diag)
        # g_swiid_annual の解析後、以降のログの last_g_swiid_obs が観測値に更新される
        idx_g = findfirst(d -> d["update_type"] == "g_swiid_annual" && d["phase"] == "post",
                          r_on.g_diag)
        @test idx_g !== nothing
        @test r_on.g_diag[idx_g]["last_g_swiid_obs"] ≈ 0.3
        for d in r_on.g_diag[(idx_g + 1):end]
            @test d["last_g_swiid_obs"] ≈ 0.3
        end
        # 解析前は g_swiid 未観測(nothing)
        @test r_on.g_diag[1]["last_g_swiid_obs"] === nothing
        # g_mean/g_sd は有限な実数(病的値でない)
        @test all(isfinite(d["g_mean"]) && isfinite(d["g_sd"]) && d["g_sd"] >= 0
                  for d in r_on.g_diag)
        # t は単調非減少(ログ順序が時系列どおり)
        @test issorted(d["t"] for d in r_on.g_diag)
        # 拡大なしのランでは "aug" キーは省略される(#0066)
        @test all(!haskey(d, "aug") for d in r_on.g_diag)
    end

    @testset "g 計装診断: 拡大パラメータ統計と予報窓系列 (#0066)" begin
        # #0065 のシナリオに L3 拡大(mu_gbar identity + theta_sig log)を加え、
        # (i) 拡大ラン有効/無効の bitwise 同一性、(ii) 各レコードの "aug"
        # 統計がリンク座標のアンサンブル実測と一致、(iii) g_aug_series の
        # 後処理抽出が X の該当断面と厳密一致、を検証する。
        params = build_params(:volatile)
        t1 = 3.0
        truth = simulate_sde(params; seed = 7201, t1 = t1)
        event_times = [e.t for e in truth.jumps]
        obs = synthesize_observations(truth.traj, standard_observations(params);
                                      rng = Xoshiro(7202))
        g_swiid_spec = ObservationSpec(:g_swiid, 1.0, 0.05, xi -> xi[IX_G], IX_G)
        obs_all = sort(vcat(obs, [ObservationRecord(2.0, g_swiid_spec, 0.3)]);
                       by = o -> o.t)
        N = 30
        aug = [AugmentedParam(name = :mu_gbar, link = :identity,
                              init = params.l2.mu_gbar, init_sd = 0.1, rw_sd = 0.02),
               AugmentedParam(name = :theta_sig, link = :log,
                              init = params.l3.theta_sig, init_sd = 0.1, rw_sd = 0.01)]
        E0s = params.x0 .+ 0.3 .* randn(Xoshiro(7203), N_STATE, N)
        postprocess_analysis!(E0s)
        E0 = MiraiYohou.augment_ensemble(E0s, aug; rng = Xoshiro(7204))
        cfg = AssimConfig(t1 = t1)
        seed = 7205

        r_off = run_assimilation(params, copy(E0), obs_all, event_times;
                                 cfg, seed, augmented_params = aug)
        r_on = run_assimilation(params, copy(E0), obs_all, event_times;
                                cfg, seed, augmented_params = aug,
                                instrument_g = true)

        # (i) 判定数値へ影響ゼロ(拡大ランでも bitwise 同一)
        @test r_off.X == r_on.X
        @test r_off.ranks == r_on.ranks
        @test r_off.ess == r_on.ess
        @test r_off.nresample == r_on.nresample
        @test isequal(r_off.count_logscore, r_on.count_logscore)
        @test isempty(r_off.g_diag)

        # (ii) 全レコードに拡大パラメータ統計("aug")が付く
        @test !isempty(r_on.g_diag)
        for d in r_on.g_diag
            @test haskey(d, "aug")
            @test Set(keys(d["aug"])) == Set(["mu_gbar", "theta_sig"])
            @test d["aug"]["mu_gbar"]["link"] == "identity"
            @test d["aug"]["theta_sig"]["link"] == "log"
            @test all(isfinite(d["aug"][k]["mean"]) && d["aug"][k]["sd"] >= 0
                      for k in ("mu_gbar", "theta_sig"))
        end

        # (iii) g_aug_series: 実行済み軌道の後処理集計が X の断面と厳密一致
        ser = MiraiYohou.g_aug_series(r_on.t, r_on.X, aug; stride = 2)
        @test !isempty(ser)
        @test ser[1]["t"] == r_on.t[1]
        @test ser[end]["t"] == r_on.t[end]   # 末端は stride 不一致でも含む
        @test all(d["phase"] == "series" && d["update_type"] == "forecast_free"
                  for d in ser)
        for (d, s) in ((ser[1], 1), (ser[end], length(r_on.t)))
            gv = vec(r_on.X[IX_G, s, :])
            @test d["g_mean"] ≈ sum(gv) / N atol = 1e-12
            @test d["g_sd"] ≈ sqrt(sum(abs2, gv .- sum(gv) / N) / (N - 1)) atol = 1e-12
            mv = vec(r_on.X[N_STATE + 1, s, :])   # mu_gbar 行(リンク座標)
            @test d["aug"]["mu_gbar"]["mean"] ≈ sum(mv) / N atol = 1e-12
        end
        # 元配列は変更されない(読み取り専用)
        Xc = copy(r_on.X)
        MiraiYohou.g_aug_series(r_on.t, r_on.X, aug)
        @test r_on.X == Xc
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
