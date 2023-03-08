using FromFile: @from
@from "../src/ThreeHalf.jl" import ThreeHalf
using Test
import LinearAlgebra: eigen

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

    esr_freqs2 = ThreeHalf.calc_esr_freqs(2e-3, 2.87e9, 0, 1 // 1)[1]
    sesr2 = sort(esr_freqs2, rev = true)
    @test (sesr2[1] - sesr2[2]) / (2 * 28e9) ≈ 2e-3 atol = 1e-2

    esr_freqs3 = ThreeHalf.calc_esr_freqs(-10e-3, 2.87e9, 0, 1 // 1)[1]
    sesr3 = sort(esr_freqs3, rev = true)
    @test (sesr3[1] - sesr3[2]) / (2 * 28e9) ≈ 10e-3 atol = 1e-2

    esr_freqs4 = ThreeHalf.calc_esr_freqs(0.2e-3, 10e6, 50e6, 3 // 2)[1]
    sesr4 = sort(esr_freqs4, rev = true)
    @test (sesr4[1] - sesr4[2]) / (2 * 28e9) ≈ 0.2e-3 atol = 1e-2

    esr_freqs5 = ThreeHalf.calc_esr_freqs(-0.2e-3, 10e6, 50e6, 3 // 2)[1]
    sesr5 = sort(esr_freqs5, rev = true)
    @test (sesr5[1] - sesr5[2]) / (2 * 28e9) ≈ 0.2e-3 atol = 1e-2
end

@testset "NV single PL" begin
    # NOTE this guy is outdated (pl fn is outdated)
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
    @test PL ≈ 1.044e7 rtol = 1e-3
end

@testset "zf. op. trans" begin
    kvb = Dict(
        # pump
        (1 // 2, 3 // 2) => 1e9,
        (1 // 2, -3 // 2) => 1e9,
        (-1 // 2, 3 // 2) => 1e9,
        (-1 // 2, -3 // 2) => 1e9,
        # need some decay to cycle...
        (3 // 2, 1 // 2) => 1e7,
        (-3 // 2, 1 // 2) => 1e7,
        (3 // 2, -1 // 2) => 1e7,
        (-3 // 2, -1 // 2) => 1e7,
    )
    resvb = ThreeHalf.zero_field_optical_transitions(kvb, 3 // 2)
    @test resvb == [0 1e7 1e7 0; 1e9 0 0 1e9; 1e9 0 0 1e9; 0 1e7 1e7 0]

    knv = Dict(
        # pump
        (1 // 1, 0 // 1) => 1e9,
        (-1 // 1, 0 // 1) => 1e9,
        # small decay
        (0 // 1, 1 // 1) => 5e8,
        (0 // 1, -1 // 1) => 5e8,
    )
    resnv = ThreeHalf.zero_field_optical_transitions(knv, 1 // 1)
    @test resnv == [0 1e9 0; 5e8 0 5e8; 0 1e9 0]
end

@testset "mix. opt. trans !" begin
    knv = Dict(
        # pump
        (1 // 1, 0 // 1) => 1e9,
        (-1 // 1, 0 // 1) => 1e9,
        # small decay
        (0 // 1, 1 // 1) => 5e8,
        (0 // 1, -1 // 1) => 5e8,
    )
    nstates = Int(2 * 1 // 1 + 1)
    k = zeros((nstates, nstates))

    _, evecs = eigen(ThreeHalf._ham_sd(1e-3, 2870e6, 0, 1 // 1))
    ThreeHalf.add_mixed_optical_transitions!(
        k,
        ThreeHalf.zero_field_optical_transitions(knv, 1 // 1),
        evecs,
    )
    @test k ≈ [280e3 5e8 5e8; 1e9 140e3 140e3; 1e9 140e3 130e3] rtol = 1e-3

    kvb = Dict(
        # pump
        (1 // 2, 3 // 2) => 1e9,
        (1 // 2, -3 // 2) => 1e9,
        (-1 // 2, 3 // 2) => 1e9,
        (-1 // 2, -3 // 2) => 1e9,
        # need some decay to cycle...
        (3 // 2, 1 // 2) => 1e7,
        (-3 // 2, 1 // 2) => 1e7,
        (3 // 2, -1 // 2) => 1e7,
        (-3 // 2, -1 // 2) => 1e7,
    )
    nstates = Int(2 * 3 // 2 + 1)
    k = zeros((nstates, nstates))

    _, evecs = eigen(ThreeHalf._ham_sd(1e-3, 5e6, 0, 3 // 2))
    ThreeHalf.add_mixed_optical_transitions!(
        k,
        ThreeHalf.zero_field_optical_transitions(kvb, 3 // 2),
        evecs,
    )
    @test k ≈
          [
        2.45 3.23 2.66 3.43
        1.90 2.50 2.05 2.66
        2.31 3.04 2.50 3.23
        1.76 2.31 1.90 2.45
    ] .* 1e8 rtol = 1e-2
end

@testset "unmix. mw trans !" begin
    k = zeros((Int(2 * 1 // 1 + 1), Int(2 * 1 // 1 + 1)))

    _, evecs = eigen(ThreeHalf._ham_sd(1e-3, 2870e6, 0, 1 // 1))

    kdrive = Dict(Set((0 // 1, 1 // 1)) => 1e10)
    ThreeHalf.add_unmixed_mw_transtitions!(k, evecs, kdrive, 1 // 1)

    @test k == [0 1e10 0; 1e10 0 0; 0 0 0]

    k = zeros((Int(2 * 3 // 2 + 1), Int(2 * 3 // 2 + 1)))

    _, evecs = eigen(ThreeHalf._ham_sd(1e-3, 5e6, 0, 3 // 2))

    kdrive = Dict(Set((1 // 2, 3 // 2)) => 1e10)
    ThreeHalf.add_unmixed_mw_transtitions!(k, evecs, kdrive, 3 // 2)

    @test k == [0 1e10 0 0; 1e10 0 0 0; 0 0 0 0; 0 0 0 0]
end

@testset "mix. mw trans !" begin
    k = zeros((Int(2 * 1 // 1 + 1), Int(2 * 1 // 1 + 1)))
    evals, evecs = eigen(ThreeHalf._ham_sd(1e-3, 2870e6, 0, 1 // 1))
    ω = 2870e6
    ls = (Ω, Δ, Γ) -> Ω^2 / (Δ^2 + Γ^2)
    ls_params = [1e10, 2e6]

    ThreeHalf.add_mixed_mw_transitions!(
        k,
        evecs,
        ω,
        1 // 1,
        evals,
        evecs,
        ls,
        ls_params,
    )
    @test isapprox(
        k,
        [0 7.4e-6 7.4e-6; 7.4e-6 0 1.5e-5; 7.4e-6 1.5e-5 0],
        rtol = 1e-1,
    )

    k = zeros((Int(2 * 3 // 2 + 1), Int(2 * 3 // 2 + 1)))
    evals, evecs = eigen(ThreeHalf._ham_sd(1e-3, 5e6, 0, 3 // 2))
    ω = 1e9
    ls = (Ω, Δ, Γ) -> Ω^2 / (Δ^2 + Γ^2)
    ls_params = [1e10, 2e6]

    ThreeHalf.add_mixed_mw_transitions!(
        k,
        evecs,
        ω,
        3 // 2,
        evals,
        evecs,
        ls,
        ls_params,
    )

    @test isapprox(
        k,
        [0 3.9 3.1 2.0; 3.9 0 4.5 2.7; 3.1 4.4 0 3.1; 2.0 2.7 3.1 0] .* 1e-3,
        rtol = 1e-1,
    )
end

@testset "solve_steady_state" begin
    knv = Dict(
        # pump
        (1 // 1, 0 // 1) => 1e9,
        (-1 // 1, 0 // 1) => 1e9,
        # small decay
        (0 // 1, 1 // 1) => 5e8,
        (0 // 1, -1 // 1) => 5e8,
    )
    k = zeros(Int(2 * 1 // 1 + 1), Int(2 * 1 // 1 + 1))
    _, evecs = eigen(ThreeHalf._ham_sd(1e-3, 2870e6, 0, 1 // 1))
    k0 = ThreeHalf.zero_field_optical_transitions(knv, 1 // 1)
    ThreeHalf.add_mixed_optical_transitions!(k, k0, evecs)

    res = ThreeHalf.solve_steady_state(k)
    @test isapprox(res, [0.50, 0.25, 0.25], rtol = 1e-3)
end

@testset "pl unmixed" begin
    knv = Dict(
        # pump
        (1 // 1, 0 // 1) => 1e9,
        (-1 // 1, 0 // 1) => 1e9,
        # small decay
        (0 // 1, 1 // 1) => 5e8,
        (0 // 1, -1 // 1) => 5e8,
    )
    k = zeros(Int(2 * 1 // 1 + 1), Int(2 * 1 // 1 + 1))
    _, evecs = eigen(ThreeHalf._ham_sd(1e-3, 2870e6, 0, 1 // 1))
    k0 = ThreeHalf.zero_field_optical_transitions(knv, 1 // 1)
    ThreeHalf.add_mixed_optical_transitions!(k, k0, evecs)

    pops = ThreeHalf.solve_steady_state(k)
    ℙ = Dict((0 // 1) => 1, (1 // 1) => 0, (-1 // 1) => 0)
    res = ThreeHalf.pl_from_unmixed_pops(pops, ℙ, 1 // 1, evecs)
    @test isapprox(res, 0.5, atol = 1e-3)
end

@testset "pl mixed" begin
    knv = Dict(
        # pump
        (1 // 2, 3 // 2) => 1e9,
        (1 // 2, -3 // 2) => 1e9,
        (-1 // 2, 3 // 2) => 1e9,
        (-1 // 2, -3 // 2) => 1e9,
        # need some decay to cycle...
        (3 // 2, 1 // 2) => 1e7,
        (-3 // 2, 1 // 2) => 1e7,
        (3 // 2, -1 // 2) => 1e7,
        (-3 // 2, -1 // 2) => 1e7,
    )
    k = zeros(Int(2 * 3 // 2 + 1), Int(2 * 3 // 2 + 1))
    _, evecs = eigen(ThreeHalf._ham_sd(1e-3, 5e6, 0, 3 // 2))
    k0 = ThreeHalf.zero_field_optical_transitions(knv, 3 // 2)
    ThreeHalf.add_mixed_optical_transitions!(k, k0, evecs)

    pops = ThreeHalf.solve_steady_state(k)
    ℙ = Dict((1 // 2) => 0, (-1 // 2) => 0, (3 // 2) => 1, (-3 // 2) => 1)
    res = ThreeHalf.pl_from_mixed_pops(pops, ℙ, 3 // 2, evecs)
    @test isapprox(res, 0.518, atol = 1e-3)
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
        ((:M, 1 // 0), (:G, 1 // 2)) => 20e6 * 0.34 / 2,
        ((:M, 1 // 0), (:G, -1 // 2)) => 20e6 * 0.34 / 2,
        ((:M, 1 // 0), (:G, -3 // 2)) => 20e6,
    )
    𝕜on = Dict(Set(((:G, 3 // 2), (:G, 1 // 2))) => 1.0e8)
    𝕜off = Dict()
    PLon = ThreeHalf.pl(𝔹, Ds, Es, S, 𝕜₀, 𝕜on)
    PLoff = ThreeHalf.pl(𝔹, Ds, Es, S, 𝕜₀, 𝕜off)
    contrast = 100 * (PLon - PLoff) / PLoff
    @test contrast ≈ 5.52 rtol = 1e-2
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
    ℙ = Dict((0 // 1) => 1, (1 // 1) => 0, (-1 // 1) => 0)
    𝕜₁ = Dict(Set((0 // 1, 1 // 1)) => 1e20)

    PL1 = ThreeHalf.pl_simple(𝔹, D, E, S, 𝕜₀, 𝕜₁, ℙ)
    @test PL1 ≈ 0.5 atol = 1e-3
    𝕜₂ = Dict()
    PL2 = ThreeHalf.pl_simple(𝔹, D, E, S, 𝕜₀, 𝕜₂, ℙ)
    @test PL2 ≈ 1.0 atol = 1e-3
end

@testset "NV GS only +𝔹" begin
    𝔹 = (0, 0, 1e-3) # T 
    D = 2.87e9
    E = 0
    S = 1 // 1
    𝕜₀ = Dict(
        # pump
        (1 // 1, 0 // 1) => 1e9,
        (-1 // 1, 0 // 1) => 1e9,
    )
    ℙ = Dict((0 // 1) => 1, (1 // 1) => 0, (-1 // 1) => 0)
    𝕜₁ = Dict(Set((0 // 1, 1 // 1)) => 1e20)

    PL1 = ThreeHalf.pl_simple(𝔹, D, E, S, 𝕜₀, 𝕜₁, ℙ)
    @test PL1 ≈ 0.5 atol = 1e-3
    𝕜₂ = Dict()
    PL2 = ThreeHalf.pl_simple(𝔹, D, E, S, 𝕜₀, 𝕜₂, ℙ)
    @test PL2 ≈ 1.0 atol = 1e-3
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

        # decay...
        (3 // 2, 1 // 2) => 1e7,
        (-3 // 2, 1 // 2) => 1e7,
        (3 // 2, -1 // 2) => 1e7,
        (-3 // 2, -1 // 2) => 1e7,
    )
    ℙ = Dict((1 // 2) => 0, (-1 // 2) => 0, (3 // 2) => 1, (-3 // 2) => 1)
    𝕜₁ = Dict()
    PL1 = ThreeHalf.pl_simple(𝔹, D, E, S, 𝕜₀, 𝕜₁, ℙ)
    @test PL1 ≈ 0.990 atol = 1e-3

    𝕜₂ = Dict(Set((1 // 2, 3 // 2)) => 1e9)
    PL2 = ThreeHalf.pl_simple(𝔹, D, E, S, 𝕜₀, 𝕜₂, ℙ)
    @test PL2 ≈ 0.981 atol = 1e-3 # not entirely sure this is exactly what we expect
end

@testset "NV GS mixed" begin
    D = 2.87e9
    E = 0
    S = 1 // 1
    𝕜₀ = Dict(
        # pump
        (1 // 1, 0 // 1) => 1e9,
        (-1 // 1, 0 // 1) => 1e9,
        # small decay
        (0 // 1, 1 // 1) => 1e3,
        (0 // 1, -1 // 1) => 1e3,
    )
    ℙ = Dict((0 // 1) => 1, (1 // 1) => 0, (-1 // 1) => 0)
    ls = (Δ, Ω, Γ) -> Ω^2 / (Δ^2 + Γ^2)
    Ω, Γ = 1e10, 1e6
    𝔹norm = 20e-4
    ω = 2.926e9
    𝔹 = ThreeHalf.rad2cart(𝔹norm, 0, 0)
    esr_freqs, esr_old_ids, _ = ThreeHalf.calc_esr_freqs(𝔹, D, E, S)
    PL_on = ThreeHalf.pl_mw(𝔹, D, E, S, 𝕜₀, ℙ, ω, ls, (Ω, Γ))
    PL_off = ThreeHalf.pl_mw(𝔹, D, E, S, 𝕜₀, ℙ, ω, ls, (0, Γ))
    contrast = (PL_on - PL_off) / PL_off
    @test isapprox(contrast, -0.083, rtol = 1e-2)
end

@testset "3 // 2 GS mixed" begin
    D = 5e6
    E = 0
    S = 3 // 2
    𝕜₀ = Dict(
        # pump
        (1 // 2, 3 // 2) => 1e9,
        (1 // 2, -3 // 2) => 1e9,
        (-1 // 2, 3 // 2) => 1e9,
        (-1 // 2, -3 // 2) => 1e9,
        # need some decay to cycle...
        (3 // 2, 1 // 2) => 1e7,
        (-3 // 2, 1 // 2) => 1e7,
        (3 // 2, -1 // 2) => 1e7,
        (-3 // 2, -1 // 2) => 1e7,
    )
    ℙ = Dict((1 // 2) => 0, (-1 // 2) => 0, (3 // 2) => 1, (-3 // 2) => 1)
    ls = (Δ, Ω, Γ) -> Ω^2 / (Δ^2 + Γ^2)
    Ω, Γ = 1e10, 1e6
    𝔹norm = 35e-3
    𝔹 = ThreeHalf.rad2cart(𝔹norm, π / 4, 0)
    ω = 9.8e8
    esr_freqs, esr_old_ids, _ = ThreeHalf.calc_esr_freqs(𝔹, D, E, S)
    PL_on = ThreeHalf.pl_mw(𝔹, D, E, S, 𝕜₀, ℙ, ω, ls, (Ω, Γ))
    PL_off = ThreeHalf.pl_mw(𝔹, D, E, S, 𝕜₀, ℙ, ω, ls, (0, Γ))
    contrast = (PL_on - PL_off) / PL_off
    @test isapprox(contrast, -0.00106, rtol = 1e-3)

    ω = 2 * 9.8e8 # 2γB resonance
    esr_freqs, esr_old_ids, _ = ThreeHalf.calc_esr_freqs(𝔹, D, E, S)
    PL_on = ThreeHalf.pl_mw(𝔹, D, E, S, 𝕜₀, ℙ, ω, ls, (Ω, Γ))
    PL_off = ThreeHalf.pl_mw(𝔹, D, E, S, 𝕜₀, ℙ, ω, ls, (0, Γ))
    contrast = (PL_on - PL_off) / PL_off
    @test isapprox(contrast, -0.000957, rtol = 1e-3)
end
