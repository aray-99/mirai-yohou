# 観測モデル(SPEC §9.1): H・R・観測スケジュール
#
# 全て変換座標で定義する。H は「状態→観測空間の関数」として書けば十分
# (接線行列は不要。EnKF はアンサンブル共分散で代用する、§9.1)。
# イベント数 N(週次バッチ)は EnKF ではなくポアソン重み(weights.jl)で扱う。

"1種類の観測の定義(スカラー観測)"
struct ObservationSpec
    name::Symbol
    period::Float64                 # 観測間隔(年)
    sd::Float64                     # 観測ノイズ標準偏差
    h::Function                     # xi(保持座標)→ 観測値
end

"""
    standard_observations(params) -> Vector{ObservationSpec}

§9.1 の観測演算子一覧(双子実験用の合成観測)。
log y のみ弱い非線形(logit 経由の w)を含む。
"""
function standard_observations(params::ModelParameters)
    l1, l2 = params.l1, params.l2
    h_logy = xi -> begin
        # log y = log A + alpha*xi_k + (1-alpha)*(xi_h + log w)
        phi = sigmoid(xi[IX_PHI])
        logA = log(l1.A0) + l2.theta_T * xi[IX_T] + l2.theta_phi * phi
        logA + l1.alpha * xi[IX_K] +
            (1 - l1.alpha) * (xi[IX_H] + log(sigmoid(xi[IX_W])))
    end
    return [
        ObservationSpec(:log_P, 1.0, 0.002, xi -> xi[IX_P]),
        ObservationSpec(:logit_w, 1.0, 0.01, xi -> xi[IX_W]),
        ObservationSpec(:log_y, 0.25, 0.01, h_logy),
        ObservationSpec(:logit_g, 5.0, 0.02, xi -> xi[IX_G]),
        ObservationSpec(:log_T, 1.0, 0.05, xi -> xi[IX_T]),
        ObservationSpec(:logit_phi, 1.0, 0.02, xi -> xi[IX_PHI]),
        ObservationSpec(:log_v, 1 / 12, 0.05, xi -> xi[IX_V]),
        ObservationSpec(:logit_tau, 3.0, 0.15, xi -> xi[IX_TAU]),
        ObservationSpec(:logit_p, 2.0, 0.10, xi -> xi[IX_PP]),
    ]
end

"1個の観測レコード(時刻・種類・値)"
struct ObservationRecord
    t::Float64
    spec::ObservationSpec
    value::Float64
end

"""
    observation_times(spec, t0, t1) -> Vector{Float64}

観測時刻の列(period 刻み。t0 は含まない: 初期時刻の観測はなし)。
"""
observation_times(spec::ObservationSpec, t0::Real, t1::Real) =
    collect((t0 + spec.period):spec.period:t1)

"""
    synthesize_observations(truth, specs; rng) -> Vector{ObservationRecord}

真値軌道 `truth::Trajectory` から §9.1 の頻度・ノイズで合成観測列を生成し、
時刻順に返す(双子実験 §13 手順2)。観測時刻は真値格子に丸める。
"""
function synthesize_observations(truth::Trajectory,
                                 specs::Vector{ObservationSpec};
                                 rng::AbstractRNG)
    t0, t1 = truth.t[1], truth.t[end]
    dt = truth.t[2] - truth.t[1]
    obs = ObservationRecord[]
    for spec in specs
        for t in observation_times(spec, t0, t1)
            k = clamp(round(Int, (t - t0) / dt) + 1, 1, length(truth.t))
            xi = @view truth.X[:, k]
            push!(obs, ObservationRecord(truth.t[k], spec,
                                         spec.h(xi) + spec.sd * randn(rng)))
        end
    end
    sort!(obs; by = o -> o.t)
    return obs
end
