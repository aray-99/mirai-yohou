# 駆動パラメータの L3 状態拡大の汎用化テスト(DECISIONS #0046)
#
# augmented_params(AugmentedParam のリスト)による拡大が、(i) 空リストで
# 従来動作と同一、(ii) augmented=true の後方互換経路が新機構経由でも
# 従来と一致、(iii) 複数記述子で行数・注入・RW 拡散が正しい、(iv) log
# リンクの正値性、(v) identity リンクの負値許容、(vi) スプレッド床/
# マスキングとの共存、を検証する。

_mean(v) = sum(v) / length(v)
_std(v) = sqrt(sum(abs2, v .- _mean(v)) / (length(v) - 1))

@testset "AugmentedParam L3 状態拡大の汎用化 (#0046)" begin

    @testset "リンク座標の相互変換" begin
        @test MiraiYohou._link_to(:log, 3.0) ≈ log(3.0)
        @test MiraiYohou._link_from(:log, log(3.0)) ≈ 3.0
        @test MiraiYohou._link_to(:identity, -0.01) == -0.01
        @test MiraiYohou._link_from(:identity, -0.01) == -0.01
        @test_throws ArgumentError MiraiYohou._link_to(:bogus, 1.0)
        @test_throws ArgumentError MiraiYohou._link_from(:bogus, 1.0)
    end

    @testset "_aug_location: 名前だけで所属構造体を判定" begin
        @test MiraiYohou._aug_location(:theta_sig) === :l3
        @test MiraiYohou._aug_location(:netgrowth) === :exo
        @test MiraiYohou._aug_location(:wbar) === :exo
        @test MiraiYohou._aug_location(:a_T) === :l2
        @test MiraiYohou._aug_location(:r_phi) === :l2
        @test MiraiYohou._aug_location(:mu_gbar) === :l2
        @test_throws ArgumentError MiraiYohou._aug_location(:not_a_field)
    end

    @testset "_with_field / _inject_param: 他フィールドを変えず1個だけ差し替える" begin
        params = build_params(:stable)
        p2 = MiraiYohou._inject_param(params, :a_T, 0.099)
        @test p2.l2.a_T == 0.099
        @test p2.l2.r_phi == params.l2.r_phi           # 他 L2 フィールドは不変
        @test p2.l3 === params.l3
        @test p2.exo === params.exo

        p3 = MiraiYohou._inject_param(params, :theta_sig, 4.5)
        @test p3.l3.theta_sig == 4.5
        @test p3.l2 === params.l2

        p4 = MiraiYohou._inject_param(params, :netgrowth, -0.02)
        @test p4.exo.netgrowth == -0.02
        @test p4.exo.wbar == params.exo.wbar           # ConstantExogenous の他フィールドは不変
        @test p4.l2 === params.l2
    end

    @testset "build_member_params: 複数記述子を順に注入する" begin
        params = build_params(:stable)
        N = 3
        E = zeros(N_STATE + 2, N)
        aug = [AugmentedParam(name = :netgrowth, link = :identity),
               AugmentedParam(name = :a_T, link = :log)]
        E[N_STATE + 1, :] .= [-0.01, 0.0, 0.02]        # netgrowth(自然単位そのまま)
        E[N_STATE + 2, :] .= log.([0.02, 0.05, 0.1])   # log a_T

        for i in 1:N
            p = MiraiYohou.build_member_params(params, aug, E, N_STATE, i)
            @test p.exo.netgrowth ≈ E[N_STATE + 1, i]
            @test p.l2.a_T ≈ exp(E[N_STATE + 2, i])
            @test p.exo.wbar == params.exo.wbar
            @test p.l2.r_phi == params.l2.r_phi
        end

        # 空リストは params をそのまま返す(後方互換)
        @test MiraiYohou.build_member_params(params, AugmentedParam[], E, N_STATE, 1) === params
    end

    @testset "augment_ensemble: 初期アンサンブル行の追加" begin
        N = 4000
        E0_state = zeros(N_STATE, N)
        aug = [AugmentedParam(name = :netgrowth, link = :identity,
                              init = -0.01, init_sd = 0.05),
               AugmentedParam(name = :a_T, link = :log, init = 0.02, init_sd = 0.1)]
        rng = Xoshiro(42)
        E0 = MiraiYohou.augment_ensemble(E0_state, aug; rng)
        @test size(E0) == (N_STATE + 2, N)
        @test E0[1:N_STATE, :] == E0_state
        @test isapprox(_mean(E0[N_STATE + 1, :]), -0.01; atol = 0.01)
        @test isapprox(_mean(E0[N_STATE + 2, :]), log(0.02); atol = 0.03)
        @test _std(E0[N_STATE + 1, :]) > 0                 # 実際に摂動されている

        # 元の状態行列は変更されない
        E0_state_copy = copy(E0_state)
        MiraiYohou.augment_ensemble(E0_state, aug; rng = Xoshiro(1))
        @test E0_state == E0_state_copy

        # 空リストは E0_state のコピーを返す(後方互換)
        E0_empty = MiraiYohou.augment_ensemble(E0_state, AugmentedParam[]; rng)
        @test E0_empty == E0_state
        @test size(E0_empty) == (N_STATE, N)
    end

    # ---- run_assimilation 統合テスト(小規模・短時間で高速化) ----

    function _twin_scenario(; seed = 7001, N = 40, t1 = 3.0)
        params = build_params(:volatile)
        truth = simulate_sde(params; seed = seed, t1 = t1)
        event_times = [e.t for e in truth.jumps]
        obs = synthesize_observations(truth.traj, standard_observations(params);
                                      rng = Xoshiro(seed + 1))
        E0 = params.x0 .+ 0.3 .* randn(Xoshiro(seed + 2), N_STATE, N)
        postprocess_analysis!(E0)
        cfg = AssimConfig(t1 = t1)
        return (; params, obs, event_times, E0, cfg, seed, N)
    end

    @testset "空の augmented_params は従来動作と同一" begin
        s = _twin_scenario()
        r1 = run_assimilation(s.params, copy(s.E0), s.obs, s.event_times;
                              cfg = s.cfg, seed = s.seed + 2)
        r2 = run_assimilation(s.params, copy(s.E0), s.obs, s.event_times;
                              cfg = s.cfg, seed = s.seed + 2,
                              augmented_params = AugmentedParam[])
        @test r1.X == r2.X
    end

    @testset "augmented=true の後方互換経路は新機構経由でも従来結果と一致" begin
        s = _twin_scenario()
        prior = l3_priors().theta_sig
        E0b = vcat(copy(s.E0),
                   reshape(log.(rand(Xoshiro(s.seed + 4), prior, s.N)), 1, :))

        r_legacy = run_assimilation(s.params, copy(E0b), s.obs, s.event_times;
                                    cfg = s.cfg, seed = s.seed + 4, augmented = true)
        aug = [AugmentedParam(name = :theta_sig, link = :log,
                              rw_sd = s.cfg.param_noise_sd)]
        r_new = run_assimilation(s.params, copy(E0b), s.obs, s.event_times;
                                 cfg = s.cfg, seed = s.seed + 4,
                                 augmented_params = aug)
        @test r_legacy.X == r_new.X
        @test r_legacy.ess == r_new.ess
        @test r_legacy.nresample == r_new.nresample
    end

    @testset "augmented=true と augmented_params の同時指定はエラー" begin
        s = _twin_scenario()
        E0b = vcat(copy(s.E0), zeros(1, s.N))
        aug = [AugmentedParam(name = :theta_sig, link = :log)]
        @test_throws ArgumentError run_assimilation(s.params, E0b, s.obs, s.event_times;
                                                    cfg = s.cfg, seed = s.seed,
                                                    augmented = true,
                                                    augmented_params = aug)
    end

    @testset "行数不一致は DimensionMismatch" begin
        s = _twin_scenario()
        aug = [AugmentedParam(name = :netgrowth), AugmentedParam(name = :a_T, link = :log)]
        E0_wrong = copy(s.E0)   # N_STATE 行のまま(2個の記述子には足りない)
        @test_throws DimensionMismatch run_assimilation(s.params, E0_wrong, s.obs,
                                                        s.event_times;
                                                        cfg = s.cfg, seed = s.seed,
                                                        augmented_params = aug)
    end

    @testset "複数記述子: 行数・注入・RW 拡散・リンクの正値/負値許容" begin
        s = _twin_scenario(; N = 60, t1 = 4.0)
        aug = [AugmentedParam(name = :netgrowth, link = :identity,
                              init = s.params.exo.netgrowth, init_sd = 0.02,
                              rw_sd = 0.02),
               AugmentedParam(name = :a_T, link = :log,
                              init = s.params.l2.a_T, init_sd = 0.1, rw_sd = 0.1)]
        rng0 = Xoshiro(s.seed + 2)
        E0 = MiraiYohou.augment_ensemble(s.E0, aug; rng = rng0)
        @test size(E0, 1) == N_STATE + 2

        res = run_assimilation(s.params, E0, s.obs, s.event_times;
                               cfg = s.cfg, seed = s.seed + 2, augmented_params = aug)
        @test size(res.X, 1) == N_STATE + 2

        # log リンク(a_T)行: 自然単位への逆変換は常に正
        @test all(>(0), exp.(res.X[N_STATE + 2, :, :]))

        # identity リンク(netgrowth)行: 状態がそのまま自然単位(負値も表現可能)
        # であることを、負の値を直接注入した member_params で確認する。
        final_netgrowth = res.X[N_STATE + 1, end, :]
        @test _std(final_netgrowth) > 0
        E_neg = copy(res.X[:, end, :])
        E_neg[N_STATE + 1, :] .= -0.05          # 明示的に負値を注入
        for i in 1:5
            p = MiraiYohou.build_member_params(s.params, aug, E_neg, N_STATE, i)
            @test p.exo.netgrowth == -0.05      # identity リンクは値をそのまま通す(負値OK)
        end
    end

    @testset "simulate_sde_augmented: 拡大込み前進(#0053)" begin
        params = build_params(:volatile)
        t1 = 2.0

        # (i) 空リストは simulate_sde(既定 EndogenousHawkes)と同一軌道
        ref = simulate_sde(params; seed = 501, t1 = t1)
        r0 = simulate_sde_augmented(params, AugmentedParam[];
                                    seed = 501, t1 = t1, xi0 = params.x0)
        @test r0.X == ref.traj.X
        @test r0.t == ref.traj.t
        @test length(r0.jumps) == length(ref.jumps)

        # (ii) 行数チェック: xi0 の長さ不一致は DimensionMismatch
        aug = [AugmentedParam(name = :netgrowth, link = :identity,
                              init = params.exo.netgrowth, rw_sd = 0.02),
               AugmentedParam(name = :a_T, link = :log,
                              init = params.l2.a_T, rw_sd = 0.1)]
        @test_throws DimensionMismatch simulate_sde_augmented(params, aug;
                                                              seed = 502, t1 = t1,
                                                              xi0 = params.x0)

        # (iii) rw_sd = 0 なら拡大行は全期間一定(RW ノイズなし)
        aug0 = [AugmentedParam(name = :netgrowth, link = :identity, rw_sd = 0.0),
                AugmentedParam(name = :a_T, link = :log, rw_sd = 0.0)]
        xi0 = vcat(params.x0, [-0.015, log(0.05)])
        rc = simulate_sde_augmented(params, aug0; seed = 503, t1 = t1, xi0)
        @test size(rc.X, 1) == N_STATE + 2
        @test all(rc.X[N_STATE + 1, :] .== -0.015)
        @test all(rc.X[N_STATE + 2, :] .== log(0.05))

        # (iv) rw_sd > 0 なら拡大行が実際に RW 拡散する
        rw = simulate_sde_augmented(params, aug; seed = 504, t1 = t1,
                                    xi0 = vcat(params.x0,
                                               [params.exo.netgrowth,
                                                log(params.l2.a_T)]))
        @test _std(rw.X[N_STATE + 1, :]) > 0
        @test _std(rw.X[N_STATE + 2, :]) > 0

        # (v) シード再現性
        rw2 = simulate_sde_augmented(params, aug; seed = 504, t1 = t1,
                                     xi0 = vcat(params.x0,
                                                [params.exo.netgrowth,
                                                 log(params.l2.a_T)]))
        @test rw.X == rw2.X

        # (vi) 拡大値が力学に効く: netgrowth(identity)の初期値を大きく変える
        # と P(IX_P)の軌道が変わる(rw_sd=0 で決定的に注入・同一シード)
        xa = vcat(params.x0, [0.05, log(params.l2.a_T)])
        xb = vcat(params.x0, [-0.05, log(params.l2.a_T)])
        aug_fix = [AugmentedParam(name = :netgrowth, link = :identity, rw_sd = 0.0),
                   AugmentedParam(name = :a_T, link = :log, rw_sd = 0.0)]
        ra = simulate_sde_augmented(params, aug_fix; seed = 505, t1 = t1, xi0 = xa)
        rb = simulate_sde_augmented(params, aug_fix; seed = 505, t1 = t1, xi0 = xb)
        @test ra.X[IX_P, end] != rb.X[IX_P, end]
        @test ra.X[IX_P, end] > rb.X[IX_P, end]   # 高い純増加率 → 大きい log P
    end

    @testset "スプレッド床・マスキングとの共存(拡大行は非対象のまま)" begin
        N = 2000
        rng = Xoshiro(11)
        E = zeros(N_STATE + 2, N)
        E[IX_TAU, :] .= 0.01 .* randn(rng, N)
        E[N_STATE + 1, :] .= -0.01 .+ 0.02 .* randn(rng, N)
        E[N_STATE + 2, :] .= log(0.02) .+ 0.05 .* randn(rng, N)
        row1_before = copy(E[N_STATE + 1, :])
        row2_before = copy(E[N_STATE + 2, :])
        spec_tau = ObservationSpec(:tau, 3.0, 0.15, xi -> xi[IX_TAU], IX_TAU)
        batch = [ObservationRecord(1.0, spec_tau, 0.0)]
        MiraiYohou.apply_obs_spread_floor!(E, batch, 0.5, rng)
        @test E[N_STATE + 1, :] == row1_before
        @test E[N_STATE + 2, :] == row2_before
        @test size(E, 1) == N_STATE + 2

        # 現在時刻マスキング(select_masked_rows)は状態行のみを対象にする既存
        # 実装であり、拡大行の行番号(N_STATE+1 以降)を含めても不変(構成の
        # 独立性を確認)。
        cfg = AssimConfig(analysis_masked_vars = [IX_TAUA])
        masked = MiraiYohou.select_masked_rows(cfg, batch)
        @test masked == [IX_TAUA]
        @test !(N_STATE + 1 in masked) && !(N_STATE + 2 in masked)
    end
end

@testset "simulate_sde_augmented: gamma_thinning_p 超過確率シンニング (#0068)" begin
    params = build_params(:volatile)
    t1 = 2.0
    aug = AugmentedParam[]

    @testset "p=1.0 は既定と bitwise 同一(乱数消費不変)" begin
        ref = simulate_sde_augmented(params, aug; seed = 601, t1 = t1, xi0 = params.x0)
        r1 = simulate_sde_augmented(params, aug; seed = 601, t1 = t1, xi0 = params.x0,
                                    gamma_thinning_p = 1.0)
        @test r1.X == ref.X
        @test r1.t == ref.t
        @test length(r1.jumps) == length(ref.jumps)
        @test all(a.t == b.t && a.rho == b.rho && a.m == b.m
                  for (a, b) in zip(r1.jumps, ref.jumps))
    end

    @testset "p=0.0 は内生ジャンプが一切発生しない" begin
        r0 = simulate_sde_augmented(params, aug; seed = 602, t1 = t1, xi0 = params.x0,
                                    gamma_thinning_p = 0.0)
        @test isempty(r0.jumps)
        # λ_e(IX_LAME)加算もなし: apply_jump! を一切呼ばないので、λ_e は
        # ジャンプ以外の要因(drift/diffusion)でのみ変化する。既定シナリオ
        # では IX_LAME の drift/diffusion 寄与は無いため、初期値のまま。
        @test all(r0.X[IX_LAME, :] .== params.x0[IX_LAME])
    end

    @testset "0<p<1 でジャンプ数が単調に減る(固定シード)" begin
        seed = 603
        n_jumps(p) = length(simulate_sde_augmented(params, aug; seed, t1 = t1,
                                                    xi0 = params.x0,
                                                    gamma_thinning_p = p).jumps)
        counts = [n_jumps(p) for p in (1.0, 0.75, 0.5, 0.25, 0.0)]
        @test issorted(counts; rev = true)
        @test counts[end] == 0
        @test counts[1] > 0
    end
end
