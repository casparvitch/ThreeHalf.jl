module Output

using FromFile: @from
@from "./ThreeHalf.jl" import ThreeHalf
using Plots
using PrettyPrint
using Formatting: printfmtln
using Random: Random
using Distributions: Distributions
import LinearAlgebra: eigen

# just for testing output of some functions, to keep run.jl somewhat clean

# Rabi plotting
if false
    taus = LinRange(0, 100e-9, 200)
    rabi_ω = π / 50e-9 # t -> rad/s
    # rabi_ω = π / 231e-9 # t -> rad/s
    detuning_mhz = [-12.5e6, 0, 12.5e6]
    # detuning_mhz = [50e6 - 12.5e6, 50e6, 50e6 + 12.5e6]
    # detuning_mhz = [50e6 - 30e6, 50e6, 50e6 + 30e6]
    # something to do with average derivative/rate? Δ c.f. Ω (natural...) 
    labels = ["L:  -1/2  ↔ -3/2", "M: -1/2  ↔  1/2", "R:   1/2  ↔  3/2"]

    detuning_ω = 2 .* π .* detuning_mhz # MHz -> ang freq

    # populations = [ThreeHalf.rabi.(taus, rabi_ω, Δ) for Δ in detuning_ω]
    populations = [ThreeHalf.dt_rabi.(taus, rabi_ω, Δ) for Δ in detuning_ω]
    plt = ThreeHalf.series_plotter(taus, populations, labels)
    title!(plt, "Driving on LHS")
    vline!(plt, [50e-9], ls = :dot, label = "π", color = :black)
    display(plt)
end

# NV ESR test
if false
    𝔹 = [0, 0, 50e-3] # T
    D = 10e6
    E = 50e6
    S = 3 // 2
    esr_freqs, esr_old_ids, esr_new_ids = ThreeHalf.calc_esr_freqs(𝔹, D, E, S)
    print("esr_freqs: ")
    pprintln(esr_freqs)
    print("esr_old_ids: ")
    pprintln(esr_old_ids)
    print("esr_new_ids: ")
    pprintln(esr_new_ids)
    sesr = sort(esr_freqs, rev = true)
    (sesr[1] - sesr[2]) |> display
    (sesr[1] - sesr[2]) / (2 * 28e9) |> display
    (sesr[1] - sesr[2]) / (2 * 28e9) - abs(𝔹[3]) |> display
end

# NV +𝔹 ESR test
if false
    𝔹 = [0, 0, 0.2e-3] # T
    D = 2870e6 # Hz
    E = 100e3
    S = 1 // 1
    esr_freqs, esr_old_ids, esr_new_ids = ThreeHalf.calc_esr_freqs(𝔹, D, E, S)
    print("esr_freqs: ")
    pprintln(esr_freqs)
    print("esr_old_ids: ")
    pprintln(esr_old_ids)
    print("esr_new_ids: ")
    pprintln(esr_new_ids)
    sesr = sort(esr_freqs, rev = true)
    (sesr[1] - sesr[2]) |> display
    (sesr[1] - sesr[2]) / (2 * 28e9) |> display
    (sesr[1] - sesr[2]) / (2 * 28e9) - abs(𝔹[3]) |> display
end

# 3 // 2 ESR test
if false
    𝔹 = [0, 0, 0.2e-3] # T
    D = 10e6 # Hz
    E = 50e6
    S = 3 // 2
    esr_freqs = ThreeHalf.calc_esr_freqs(𝔹, D, E, S)[1]
    sesr = sort(esr_freqs, rev = true)
    (sesr[1] - sesr[2]) / (2 * 28e9) |> display
end

# NV :G, :E, :M
if false
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
    display(PL)
end

# 3//2 :G, :E, :M
if false
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
    printfmtln("𝐶: {:.2f}% ", contrast)
end

end # module

nothing
