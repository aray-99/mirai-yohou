using MiraiYohou
using Random
using Test
using TOML

# 合格しきい値(§13/§0.5.8)。M0 では読み込み可能性のみ検証し、E1(M5)で使用する。
const THRESHOLDS = TOML.parsefile(joinpath(@__DIR__, "acceptance_thresholds.toml"))

@testset "MiraiYohou" begin
    @testset "acceptance thresholds file" begin
        @test haskey(THRESHOLDS, "E1")
        @test haskey(THRESHOLDS, "E1b")
        @test THRESHOLDS["E1"]["rmse_ratio_max"] == 0.5
        @test THRESHOLDS["E1"]["sigma_s_time_correlation_min"] == 0.6
        @test THRESHOLDS["E1"]["coverage_min"] == 0.80
        @test THRESHOLDS["E1"]["coverage_max"] == 0.98
        @test THRESHOLDS["E1b"]["theta_sig_relative_error_max"] == 0.5
    end

    include("test_coordinates.jl")
    include("test_branching.jl")
    include("test_invariants.jl")
    include("test_jump_map.jl")
    include("test_hawkes_stat.jl")
    include("test_reproducibility.jl")
    include("test_regimes.jl")
    include("test_ensemble.jl")
end
