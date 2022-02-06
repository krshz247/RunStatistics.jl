# This file is a part of RunStatistics.jl, licensed under the MIT License (MIT).

using RunStatistics
using Test
using Distributions

# I measured the timing for N=96 with the exact and approximate solution. See if julia is really faster. It seems some experiment with threading was done but without threads, I'd be surprised if it's faster than the C++ code built in release mode

@testset "squares_approx" begin

    T_obs = 20
    L = 800
    N = 80
    n = 10
    Ns = [N, n]
    epsp_1 = 10^(-2)
    epsp_2 = 0
        

    @test RunStatistics.squares_pvalue_approx(T_obs, Ns, epsp_1) ≈ 0.07082230261169509
    @test RunStatistics.squares_cdf_approx(T_obs, Ns, epsp_1) ≈ 0.9291776973883049

    @test RunStatistics.squares_cdf_approx(T_obs, L, epsp_1) ≈ 0.9291776973883049
    @test RunStatistics.squares_pvalue_approx(T_obs, L, epsp_1) ≈ 0.07082230261169509


    @test RunStatistics.squares_pvalue_approx(T_obs, Ns, epsp_2) == 0.0708223271587789
    @test RunStatistics.squares_cdf_approx(T_obs, Ns, epsp_2) == 0.9291776728412211

    @test RunStatistics.squares_cdf_approx(T_obs, L, epsp_2) == 0.9291776728412211
    @test RunStatistics.squares_pvalue_approx(T_obs, L, epsp_2) == 0.0708223271587789
end

@testset "squares_approx_split" begin

    Tobs = 15.5
    N = 12
    n = 2

    F = squares_cdf(Tobs, N)
    approx = F*F / (1.0 + RunStatistics.Delta(Tobs, N, N))

    @test abs(approx - squares_cdf_approx(Tobs, [N, n])) <= 1e-15
end 

function test_on_grid(K::Integer, N::Integer, n::Real, eps::Real)

    for i in 1:K
        Tobs = 20 + 2 * i
        @test abs(squares_cdf_approx(Tobs, [N, n]) - squares_cdf(Tobs, n * N)) <= eps
    end 
end

@testset "squares_approx_grid" begin

    test_on_grid(15, 40, 2, 2e-7)
    # test_on_grid(15, 40, 3, 4e-7) expensive
end 

@testset "squares_approx_hH" begin
    
    Tobs = 15.5 
    x = 3.3
    N = 12 

    @test abs(RunStatistics.h(Tobs, N) - 0.000373964) <= 1e-8
    @test abs(RunStatistics.H(Tobs-x, Tobs, N) -  0.00245352) <= 1e-8
    @test abs(RunStatistics.Delta(Tobs, N, N) - 0.00175994) <= 1e-8
end

@testset "squares_approx_cdf" begin
    
    @test abs(cdf(Chisq(12), 15.5) - 0.784775) <= 1e-6
end

@testset "squares_approx_interpolate" begin
    Tobs = 32

    F71 = squares_cdf_approx(Tobs, [71, 5])
    F88 = squares_cdf_approx(Tobs, [88, 4])
    F89 = squares_cdf_approx(Tobs, [89, 4])
    F100 = squares_cdf_approx(Tobs, [Int64(100), 3.55])

    @test abs(F71 - F100) <= 3e-9
    @test abs(F71 - (F88 + 3*F89)/4) <= 2e-9
end

@testset "squares_approx_bound_error" begin

    Tobs = 8
    @test abs(squares_cdf(Tobs, 100) - (squares_cdf(Tobs, 50) ^ 2 - squares_cdf(Tobs, 100) * RunStatistics.Delta(Tobs, 100, 100))) <= 1e-3
end