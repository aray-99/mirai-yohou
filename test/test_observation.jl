# tauA 擬似観測テスト(DECISIONS #0036)
#
# augment_tauA_pseudo は既定オフ(mult=0)で従来動作を保ち、mult>0 では
# tau 観測に同期した tauA 擬似観測を batch に追加する(元配列は不変)。
# 追加した擬似観測を実際に enkf_analysis! に通し、tauA アンサンブル平均が
# tau 観測値方向へ引かれることも確認する(方向性の検証で十分)。

@testset "tauA pseudo observation (#0036)" begin
    spec_tau = ObservationSpec(:tau, 3.0, 0.15, xi -> xi[IX_TAU])
    spec_other = ObservationSpec(:p, 2.0, 0.10, xi -> xi[IX_PP])
    batch = [ObservationRecord(1.0, spec_tau, 0.2),
             ObservationRecord(1.0, spec_other, -0.3)]

    @testset "mult = 0 は従来動作(no-op)" begin
        out = MiraiYohou.augment_tauA_pseudo(batch, 0.0)
        @test out === batch
        @test length(out) == 2
    end

    @testset "mult > 0 は tauA 擬似観測を追加し元配列を破壊しない" begin
        out = MiraiYohou.augment_tauA_pseudo(batch, 3.0)
        @test length(batch) == 2                 # 元配列は不変
        @test length(out) == 3
        pseudo = out[end]
        @test pseudo.spec.name === :tauA_pseudo
        @test pseudo.spec.sd ≈ 3.0 * spec_tau.sd
        @test pseudo.value == 0.2
        @test pseudo.t == 1.0
        xi = zeros(N_STATE)
        xi[IX_TAUA] = 1.23
        @test pseudo.spec.h(xi) == 1.23

        # tau 観測を含まない batch には擬似観測を追加しない
        out2 = MiraiYohou.augment_tauA_pseudo([ObservationRecord(1.0, spec_other, -0.3)], 3.0)
        @test length(out2) == 1
    end

    @testset "解析後の tauA アンサンブル平均が tau 観測値方向へ引かれる" begin
        rng = Xoshiro(2026)
        N = 4000
        E = zeros(N_STATE, N)
        E[IX_TAU, :] .= 0.1 .* randn(rng, N)          # tau ≈ 0(観測値と整合)
        E[IX_TAUA, :] .= 5.0 .+ 0.1 .* randn(rng, N)  # tauA は遠方に漂流(#0035想定)

        batch2 = [ObservationRecord(1.0, spec_tau, 0.0)]
        aug = MiraiYohou.augment_tauA_pseudo(batch2, 3.0)
        @test length(aug) == 2

        yobs = [o.value for o in aug]
        R = Diagonal([o.spec.sd^2 for o in aug]) |> Matrix
        hfun = col -> [o.spec.h(col) for o in aug]

        mean_before = sum(E[IX_TAUA, :]) / N
        enkf_analysis!(E, yobs, hfun, R; rng, rho_inf = 1.0)
        mean_after = sum(E[IX_TAUA, :]) / N

        @test mean_after < mean_before    # 5.0 から観測値 0.0 方向へ
        # カルマンゲイン K = Pxy/(Pyy+R) ≈ 0.01/(0.01+0.45^2) ≈ 0.047 なので
        # 単回解析の引き幅は ΔK ≈ 0.047×5 ≈ 0.23。有意に引かれていることのみ確認する。
        @test mean_after < 5.0 - 0.1      # 有意に引かれていること
    end
end

@testset "select_masked_rows (#0040-(α))" begin
    spec_tau = ObservationSpec(:tau, 3.0, 0.15, xi -> xi[IX_TAU])
    spec_pseudo = ObservationSpec(:tauA_pseudo, 3.0, 0.45, xi -> xi[IX_TAUA])
    spec_other = ObservationSpec(:p, 2.0, 0.10, xi -> xi[IX_PP])

    @testset "既定(analysis_masked_vars 空)は常に無マスク" begin
        cfg = AssimConfig()
        @test isempty(cfg.analysis_masked_vars)
        batch = [ObservationRecord(1.0, spec_other, -0.3)]
        @test MiraiYohou.select_masked_rows(cfg, batch) == Int[]
        batch_tau = [ObservationRecord(1.0, spec_tau, 0.2)]
        @test MiraiYohou.select_masked_rows(cfg, batch_tau) == Int[]
    end

    @testset "analysis_masked_vars 指定 + unmask 観測なし → マスク適用" begin
        cfg = AssimConfig(analysis_masked_vars = [IX_TAUA],
                          analysis_unmask_names = [:tau, :tauA_pseudo])
        batch = [ObservationRecord(1.0, spec_other, -0.3)]
        @test MiraiYohou.select_masked_rows(cfg, batch) == [IX_TAUA]
    end

    @testset "batch に unmask 観測名が含まれれば解除" begin
        cfg = AssimConfig(analysis_masked_vars = [IX_TAUA],
                          analysis_unmask_names = [:tau, :tauA_pseudo])
        batch_tau = [ObservationRecord(1.0, spec_tau, 0.2),
                     ObservationRecord(1.0, spec_other, -0.3)]
        @test MiraiYohou.select_masked_rows(cfg, batch_tau) == Int[]
        batch_pseudo = [ObservationRecord(1.0, spec_pseudo, 0.2)]
        @test MiraiYohou.select_masked_rows(cfg, batch_pseudo) == Int[]
        # unmask 対象外の観測のみでは解除されない
        batch_other = [ObservationRecord(1.0, spec_other, -0.3)]
        @test MiraiYohou.select_masked_rows(cfg, batch_other) == [IX_TAUA]
    end
end

@testset "obs spread floor (#0043)" begin
    _sd(v) = sqrt(sum(abs2, v .- (sum(v) / length(v))) / (length(v) - 1))
    N = 8000
    spec_tau = ObservationSpec(:tau, 3.0, 0.15, xi -> xi[IX_TAU], IX_TAU)   # sd=0.15
    spec_pp  = ObservationSpec(:p, 2.0, 0.10, xi -> xi[IX_PP], IX_PP)       # sd=0.10
    spec_synth = ObservationSpec(:log_y, 0.25, 0.01, xi -> xi[IX_T])       # target_ix既定=0

    make_E() = begin
        rng = Xoshiro(1)
        E = zeros(N_STATE, N)
        E[IX_TAU, :] .= 0.01 .* randn(rng, N)   # sd ≈ 0.01 ≪ floor=0.075
        E[IX_PP, :]  .= 0.5 .+ 0.5 .* randn(rng, N)  # sd ≈ 0.5 ≫ floor=0.05
        E
    end

    @testset "既定(0.0)は従来同一(no-op)" begin
        rng = Xoshiro(2)
        E = make_E()
        E0 = copy(E)
        batch = [ObservationRecord(1.0, spec_tau, 0.0)]
        MiraiYohou.apply_obs_spread_floor!(E, batch, 0.0, rng)
        @test E == E0
    end

    @testset "床未満の行は sd ≈ floor まで拡大される" begin
        rng = Xoshiro(3)
        E = make_E()
        batch = [ObservationRecord(1.0, spec_tau, 0.0)]
        MiraiYohou.apply_obs_spread_floor!(E, batch, 0.5, rng)
        floor = 0.5 * spec_tau.sd
        sd_after = _sd(E[IX_TAU, :])
        @test isapprox(sd_after, floor; rtol = 0.08)
    end

    @testset "床以上の行は不変" begin
        rng = Xoshiro(4)
        E = make_E()
        row_before = copy(E[IX_PP, :])
        batch = [ObservationRecord(1.0, spec_pp, 0.0)]
        MiraiYohou.apply_obs_spread_floor!(E, batch, 0.5, rng)
        @test E[IX_PP, :] == row_before
    end

    @testset "target_ix = 0(合成観測)は対象外" begin
        rng = Xoshiro(5)
        E = zeros(N_STATE, N)
        E[IX_T, :] .= 0.001 .* randn(rng, N)   # 極小スプレッド
        row_before = copy(E[IX_T, :])
        batch = [ObservationRecord(1.0, spec_synth, 0.0)]
        MiraiYohou.apply_obs_spread_floor!(E, batch, 0.5, rng)
        @test E[IX_T, :] == row_before
    end

    @testset "同一行を対象にする複数観測は床の大きい方1回のみ" begin
        rng = Xoshiro(6)
        E = zeros(N_STATE, N)
        E[IX_TAU, :] .= 0.01 .* randn(rng, N)
        spec_small = ObservationSpec(:tau_small, 3.0, 0.05, xi -> xi[IX_TAU], IX_TAU)
        spec_big   = ObservationSpec(:tau_big, 3.0, 0.15, xi -> xi[IX_TAU], IX_TAU)
        batch = [ObservationRecord(1.0, spec_small, 0.0),
                 ObservationRecord(1.0, spec_big, 0.0)]
        MiraiYohou.apply_obs_spread_floor!(E, batch, 0.5, rng)
        sd_after = _sd(E[IX_TAU, :])
        @test isapprox(sd_after, 0.5 * spec_big.sd; rtol = 0.08)
    end

    @testset "augmented 構成(14行)でも theta_sig 行に触れない" begin
        rng = Xoshiro(7)
        E = zeros(N_STATE + 1, N)
        E[IX_TAU, :] .= 0.01 .* randn(rng, N)
        E[end, :] .= 1.23 .+ 0.001 .* randn(rng, N)
        row_before = copy(E[end, :])
        batch = [ObservationRecord(1.0, spec_tau, 0.0)]
        MiraiYohou.apply_obs_spread_floor!(E, batch, 0.5, rng)
        @test E[end, :] == row_before
        @test size(E, 1) == N_STATE + 1
    end
end
