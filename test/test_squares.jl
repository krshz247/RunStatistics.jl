# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

using RunStatistics
using Test

@testset "squares" begin
    T_obs = 3.4
    N = 30

    @test squares_pvalue(T_obs, N) ≈ 0.818208032939994
    @test squares_cdf(T_obs, N) ≈ 0.18179196706000597
end

@testset "squares_paper" begin
    alpha = 0.001

    @test abs(squares_pvalue(15.5, 5) - alpha) <= 0.00002
    @test abs(squares_pvalue(23.8, 50) - alpha) <= 0.00002
    @test abs(squares_pvalue(25.6, 100) - alpha) <= 0.00008
end

@testset "squares_mathematica" begin
    n = 20
    T = [2, 5, 10, 20, 50]
    P = [0.8936721808595665, 0.42457437866154357, 0.06934906413527009, 0.0014159488909252227, 9.575188641974819e-9]
    eps = [1e-15, 1e-15, 1e-15, 1e-15, 1e-15]

    for i in 1:5
        @test abs(squares_pvalue(T[i], n) - P[i]) <= eps[i]
    end
end

@testset "squares_critical" begin
    N = 50 

    @test abs(squares_pvalue(19.645, 2 * N) - 0.01) <= 3e-5
    @test abs(squares_pvalue(15.34, 2 * N) - 0.05) <= 1e-4
end