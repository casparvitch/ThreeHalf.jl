using FromFile: @from
@from "../src/ThreeHalf.jl" import ThreeHalf
using Test

@testset "spinop" begin
    @test ThreeHalf.spinop(1 // 2)[1] ≈ (1 / 2) * [0 1; 1 0] atol = 1e-5
    @test ThreeHalf.spinop(1 // 2)[2] ≈ (1 / 2) * [0 -1im; 1im 0] atol = 1e-5
    @test ThreeHalf.spinop(1 // 2)[3] ≈ (1 / 2) * [1 0; 0 -1] atol = 1e-5
    @test ThreeHalf.spinop(1 // 1)[1] ≈ (1 / √(2)) * [0 1 0; 1 0 1; 0 1 0] atol =
        1e-5
    @test ThreeHalf.spinop(1 // 1)[2] ≈
          (1 / √(2)) * [0 -1im 0; 1im 0 -1im; 0 1im 0] atol = 1e-5
    @test ThreeHalf.spinop(1 // 1)[3] ≈ [1 0 0; 0 0 0; 0 0 -1] atol = 1e-5
end

@testset "calc_esr_freqs" begin
    esr_freqs = ThreeHalf.calc_esr_freqs(-2e-3, 2.87e9, 0, 1 // 1)[1]
    sesr = sort(esr_freqs, rev = true)
    @test (sesr[1] - sesr[2]) / (2 * 28e9) ≈ 2e-3 atol = 1e-2

    esr_freqs = ThreeHalf.calc_esr_freqs(-10e-3, 2e7, 0, 1 // 1)[1]
    sesr = sort(esr_freqs, rev = true)
    @test (sesr[1] - sesr[2]) / (2 * 28e9) ≈ 10e-3 atol = 1e-2
end

@testset "NV single PL" begin
    𝔹 = (0, 0, -1e-3) # T 
    Ds = (2.87e9, 1.42e9)
    Es = (0, 0)
    S = 1 // 1
    𝕜₀ = Dict(
        # pump & radiatiative spin-conserving rates
        ((:G, 1 // 1), (:E, 1 // 1)) => 1e9,
        ((:E, 1 // 1), (:G, 1 // 1)) => 65e6, # 1/(11e-6),
        ((:G, 0 // 1), (:E, 0 // 1)) => 1e9,
        ((:E, 0 // 1), (:G, 0 // 1)) => 65e6, # 1/(11e-6),
        ((:G, -1 // 1), (:E, -1 // 1)) => 1e9,
        ((:E, -1 // 1), (:G, -1 // 1)) => 65e6, # 1/(11e-6),

        # ISC transitions
        ((:E, 1 // 1), (:M, 1 // 0)) => 80e6, # 2.03e9,
        ((:E, 0 // 1), (:M, 1 // 0)) => 11e6, # 1.01e9,
        ((:E, -1 // 1), (:M, 1 // 0)) => 80e6, # 2.03e9,
        ((:M, 1 // 0), (:G, 1 // 1)) => 2.6e6, # 20e6 * 0.34/2,
        ((:M, 1 // 0), (:G, 0 // 1)) => 3.0e6, # 20e6,
        ((:M, 1 // 0), (:G, -1 // 1)) => 2.6e6, # 20e6 * 0.34/2, 
    )
    𝕜ᵢ = Dict(
        Set(((:G, 0 // 1), (:G, 1 // 1))) => 1.0e10,
        Set(((:G, 0 // 1), (:G, -1 // 1))) => 0.0e9,
        Set(((:G, -1 // 1), (:G, 1 // 1))) => 0.0e9,
    )
    PL = ThreeHalf.pl(𝔹, Ds, Es, S, 𝕜₀, 𝕜ᵢ)
    @test PL ≈ 1.044e7 rtol=1e-3
end

@testset "3//2 contrast" begin
    𝔹 = (0, 0, -50e-3) # T 
    Ds = (10e6, 5e6) 
    Es = (50e6, 65e6)
    S = 3 // 2
    𝕜₀ = Dict(
        # pump & radiatiative spin-conserving rates
        # copy vb
        ((:G, 3 // 2), (:E, 3 // 2)) => 1e9,
        ((:E, 3 // 2), (:G, 3 // 2)) => 0.7e9,
        ((:G, 1 // 2), (:E, 1 // 2)) => 1e9,
        ((:E, 1 // 2), (:G, 1 // 2)) => 0.7e9,
        ((:G, -1 // 2), (:E, -1 // 2)) => 1e9,
        ((:E, -1 // 2), (:G, -1 // 2)) => 0.7e9,
        ((:G, -3 // 2), (:E, -3 // 2)) => 1e9,
        ((:E, -3 // 2), (:G, -3 // 2)) => 0.7e9,
        
        # ISC transitions CHOOSE 3 // 2 AS FASTER (init in 3 // 2)
        ((:E, 3 // 2), (:M, 1 // 0)) => 1.3e9,
        ((:E, 1 // 2), (:M, 1 // 0)) => 0.3e9,
        ((:E, -1 // 2), (:M, 1 // 0)) => 0.3e9,
        ((:E, -3 // 2), (:M, 1 // 0)) => 1.3e9,
        
        # once again, ≈ vb
        ((:M, 1 // 0), (:G, 3 // 2)) => 20e6,
        ((:M, 1 // 0), (:G, 1 // 2)) => 20e6 * 0.34/2,
        ((:M, 1 // 0), (:G, -1 // 2)) => 20e6 * 0.34/2,
        ((:M, 1 // 0), (:G, -3 // 2)) => 20e6, 
    )
    𝕜on = Dict(
        Set(((:G, 3 // 2), (:G, 1 // 2))) => 1.0e8,
    )
    𝕜off = Dict()
    PLon = ThreeHalf.pl(𝔹, Ds, Es, S, 𝕜₀, 𝕜on)
    PLoff = ThreeHalf.pl(𝔹, Ds, Es, S, 𝕜₀, 𝕜off)
    contrast = 100*(PLon - PLoff) / PLoff
    @test contrast ≈ 5.52 rtol=1e-2
end

@testset "NV GS only" begin
    𝔹 = (0, 0, -1e-3) # T 
    D = 2.87e9
    E = 0
    S = 1 // 1
    𝕜₀ = Dict(
        # pump
        (1 // 1, 0 // 1) => 1e9,
        (-1 // 1, 0 // 1) => 1e9,
    )
    ℙ = Dict(
        (0 // 1) => 1,
        (1 // 1) => 0,
        (-1 // 1) => 0
    )
    𝕜₁ = Dict(Set((0 // 1, 1 // 1)) => 1e20)

    PL1 = ThreeHalf.pl_simple(𝔹, D, E, S, 𝕜₀, 𝕜₁, ℙ)
    @test PL1 ≈ 0.5 atol=1e-3
    𝕜₂ = Dict()
    PL2 = ThreeHalf.pl_simple(𝔹, D, E, S, 𝕜₀, 𝕜₂, ℙ)
    @test PL2 ≈ 1.0 atol=1e-3
end

@testset "3 // 2 GS only" begin
    𝔹 = (0, 0, -50e-3) # T 
    D = 10e6
    E = 0
    S = 3 // 2
    𝕜₀ = Dict(
        # pump
        (1 // 2, 3 // 2) => 1e9,
        (1 // 2, -3 // 2) => 1e9,
        (-1 // 2, 3 // 2) => 1e9,
        (-1 // 2, -3 // 2) => 1e9,
    )
    ℙ = Dict(
        (1 // 2) => 0,
        (-1 // 2) => 0,
        (3 // 2) => 1,
        (-3 // 2) => 1
    )
    𝕜₁ = Dict()
    PL1 = ThreeHalf.pl_simple(𝔹, D, E, S, 𝕜₀, 𝕜₁, ℙ)
    @test PL1 ≈ 1.0 atol=1e-6
    
    𝕜₂ = Dict(
        Set((-1 // 2, -3 // 2)) => 1e9,
        Set((1 // 2, -3 // 2)) => 1e9,
        Set((1 // 2, 3 // 2)) => 1e9,
        Set((-1 // 2, 3 // 2)) => 1e9,
    )
    PL2 = ThreeHalf.pl_simple(𝔹, D, E, S, 𝕜₀, 𝕜₂, ℙ)
    @test PL2 ≈ 2/3 atol=1e-6 # not entirely sure this is exactly what we expect
end
    
