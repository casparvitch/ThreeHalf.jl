module Run

using FromFile: @from
@from "./ThreeHalf.jl" import ThreeHalf
using Plots
using PrettyPrint
using Formatting: printfmtln
using Random: Random
using Distributions: Distributions
import LinearAlgebra: eigen

gr()
# unicodeplots()

# NV :G only
if false
    𝔹 = (0, 0, 1e-3) # T 
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
    # 1e20 and 1e9 drive the same transition -> should be in same element?
    𝕜ᵢ = Dict(Set((0 // 1, 1 // 1)) => 1e10)
    ℙ = Dict((0 // 1) => 1, (1 // 1) => 0, (-1 // 1) => 0)
    PL = ThreeHalf.pl_simple(𝔹, D, E, S, 𝕜₀, 𝕜ᵢ, ℙ)
    pprintln(PL)
    @assert ≈(PL, 1.0)
end

# 3 // 2 :G only
if false
    𝔹 = (0, 0, 50e-3) # T 
    D = 10e6
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

    𝕜ᵢ = Dict(
        #         Set((-1 // 2, -3 // 2)) => 1e9,
        Set((1 // 2, 3 // 2)) => 1e9,
        # Set((1 // 2, -3 // 2)) => 1e9,
        #         Set((-1 // 2, 3 // 2)) => 1e9,
    )
    ℙ = Dict((1 // 2) => 0, (-1 // 2) => 0, (3 // 2) => 1, (-3 // 2) => 1)
    PLon = ThreeHalf.pl_simple(𝔹, D, E, S, 𝕜₀, 𝕜ᵢ, ℙ)
    PLoff = ThreeHalf.pl_simple(𝔹, D, E, S, 𝕜₀, Dict(), ℙ)
    contrast = 100 * (PLon - PLoff) / PLoff
    printfmtln("𝐶: {:.2f}% ", contrast)
end

# NV :G Δ sweep
if false
    𝔹 = (0, 1e-3, 1e-3) # T 
    D = 2.87e9
    E = 0
    S = 1 // 1
    𝕜₀ = Dict(
        # pump
        (1 // 1, 0 // 1) => 1e9,
        (-1 // 1, 0 // 1) => 1e9,
    )
    ℙ = Dict((0 // 1) => 1, (1 // 1) => 0, (-1 // 1) => 0)

    # simples for now
    function lineshape(Δ, Ω, Γ)
        return Ω^2 / (Δ^2 + Γ^2)
    end

    NVALS = 250
    𝐶s = Vector{Float64}(undef, NVALS)
    ωs = collect(LinRange(2.87e9 - 0.1e9, 2.87e9 + 0.1e9, NVALS))

    Γ = 1.5e6
    Ω = 1e9

    esr_freqs, esr_old_ids, _ = ThreeHalf.calc_esr_freqs(𝔹, D, E, S)

    for (i, ω) in pairs(IndexLinear(), ωs)
        𝕜ᵢ = Dict()
        for (res_id, ω₀) in zip(esr_old_ids, esr_freqs)
            𝕜ᵢ[res_id] = lineshape(ω - ω₀, Ω, Γ)
        end

        PL_on = ThreeHalf.pl_simple(𝔹, D, E, S, 𝕜₀, 𝕜ᵢ, ℙ)
        PL_off = ThreeHalf.pl_simple(𝔹, D, E, S, 𝕜₀, Dict(), ℙ)
        𝐶s[i] = (PL_on - PL_off) / (PL_off == 0 ? 1 : PL_off)
    end
    plt = plot(ωs, 𝐶s, seriestype = :scatter)
    vline!(plt, [2.87e9], linestyle = :dash)
    display(plt)
end

# 3 // 2 :G Δ sweep
if false
    # 𝔹 = (0, 0, -1e-3) # T 
    𝔹 = (0, 0, -10e-3) # T 
    #     defined as half, to get correct splitting
    #     D = 10e6
    D = 10e6 / 2
    #     D = 1e9 / 2
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
    # ℙ = Dict((1 // 2) => 1, (-1 // 2) => 1, (3 // 2) => 0, (-3 // 2) => 0)

    function lineshape(Δ, Ω, Γ)
        #         return Ω * exp(-1 * Δ^2 / (2 * Γ^2))
        return Ω^2 / (Δ^2 + Γ^2)
        #         return  Ω * exp(-abs(Δ) / Ω)
    end

    esr_freqs, esr_old_ids, _ = ThreeHalf.calc_esr_freqs(𝔹, D, E, S)
    print("esr_freqs: ")
    println(esr_freqs)
    print("esr_old_ids: ")
    println(esr_old_ids)
    print("<esr_freqs>: ")
    println(sum(esr_freqs) / length(esr_freqs))

    NVALS = 1000
    #     WIDTH = 2e9
    #         WIDTH = 1e9
    # WIDTH = 0.05e9
    WIDTH = 0.5e9
    #     WIDTH = 0.0005e9
    ωs = collect(
        LinRange(
            sum(esr_freqs) / length(esr_freqs) - WIDTH,
            sum(esr_freqs) / length(esr_freqs) + WIDTH,
            NVALS,
        ),
    )
    𝐶s = Vector{Float64}(undef, NVALS)

    #         Γ = 0
    #     Γ = 1e3
    #  Γ = 5e6
    Γ = 1e7
    #     Γ = 1.5e7
    Ω = 1e9
    #     Ω = 1e7
    #     Ω = 1e11
    for (i, ω) in pairs(IndexLinear(), ωs)
        𝕜ᵢ = Dict()

        for (res_id, ω₀) in zip(esr_old_ids, esr_freqs)
            𝕜ᵢ[res_id] = lineshape(ω - ω₀, Ω, Γ)
        end

        PL_on = ThreeHalf.pl_simple(𝔹, D, E, S, 𝕜₀, 𝕜ᵢ, ℙ)
        PL_off = ThreeHalf.pl_simple(𝔹, D, E, S, 𝕜₀, Dict(), ℙ)
        𝐶s[i] = (PL_on - PL_off) / (PL_off == 0 ? 1 : PL_off)
    end
    plt = plot(ωs, 𝐶s, seriestype = :scatter)
    _, idx_min = findmin(𝐶s)
    display(plt)
end

# NV :G Δ sweep, <θ>
if false
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

    function lineshape(Δ, Ω, Γ)
        #         return Ω * exp(-1 * Δ^2 / (2 * Γ^2))
        #         return  Ω * exp(-abs(Δ) / Ω)
        return Ω^2 / (Δ^2 + Γ^2)
    end

    𝔹norm = 20e-4

    Nθs = 50
    Nωs = 250
    # WIDTH = 5e9
    # WIDTH = 1e9
    WIDTH = 0.2e9
    ωs = collect(LinRange(D - WIDTH / 2, D + WIDTH / 2, Nωs))

    cs = zeros(Nωs)
    sig = zeros(Nωs)
    ref = zeros(Nωs)
    # ISO 80000-2:2019 convention, θ polar, ϕ azimuthal
    θs = LinRange(0, π, Nθs)
    # θs = LinRange(-π/4, π/4, Nθs)
    # θs = [0, pi/2]

    Γ = 1e6
    Ω = 1e10
    # Ω = 5e11

    for θ in θs
        # u = ThreeHalf.rad2cart(1, θ, 0)
        # 𝔹 = sum([0, 0, 𝔹norm] .* u) .* u
        # 𝔹 = 𝔹norm .* u
        # 𝔹 = sum([0, 0, 𝔹norm] .* u) .* [0, 0, 1]
        𝔹 = ThreeHalf.rad2cart(𝔹norm, θ, 0)
        # FIXME cleanup?

        esr_freqs, esr_old_ids, _ = ThreeHalf.calc_esr_freqs(𝔹, D, E, S)

        for (i, ω) in pairs(IndexLinear(), ωs)
            𝕜ᵢ = Dict()
            for (res_id, ω₀) in zip(esr_old_ids, esr_freqs)
                𝕜ᵢ[res_id] = lineshape(ω - ω₀, Ω, Γ)
            end
            PL_on = ThreeHalf.pl_simple(𝔹, D, E, S, 𝕜₀, 𝕜ᵢ, ℙ)
            PL_off = ThreeHalf.pl_simple(𝔹, D, E, S, 𝕜₀, Dict(), ℙ)
            sig[i] += PL_on
            ref[i] += PL_off
        end
    end
    cs .+= (sig .- ref) ./ ref
    plt = plot(ωs, cs, seriestype = :scatter)
    vline!(plt, [2.87e9], linestyle = :dash)
    display(plt)
end

# NV :G Δ sweep, <𝔹>
if false
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

    function lineshape(Δ, Ω, Γ)
        #         return Ω * exp(-1 * Δ^2 / (2 * Γ^2))
        #         return  Ω * exp(-abs(Δ) / Ω)
        return Ω^2 / (Δ^2 + Γ^2)
    end

    𝔹norm = 20e-4

    NBs = 250
    Nωs = 250
    # WIDTH = 5e9
    # WIDTH = 1e9
    WIDTH = 0.2e9
    ωs = collect(LinRange(D - WIDTH / 2, D + WIDTH / 2, Nωs))

    cs = zeros(Nωs)
    sig = zeros(Nωs)
    ref = zeros(Nωs)
    # ISO 80000-2:2019 convention, θ polar, ϕ azimuthal
    Random.seed!(1234)
    d = Distributions.Normal(0, 20e-4 / 3)
    # Bnorms = LinRange(0, 20e-4, NBs)
    Bnorms = Random.rand(d, NBs)
    θ = 0

    Γ = 1e6
    Ω = 1e10
    # Ω = 5e11

    for Bnorm in Bnorms
        𝔹 = ThreeHalf.rad2cart(Bnorm, θ, 0)

        esr_freqs, esr_old_ids, _ = ThreeHalf.calc_esr_freqs(𝔹, D, E, S)

        for (i, ω) in pairs(IndexLinear(), ωs)
            𝕜ᵢ = Dict()
            for (res_id, ω₀) in zip(esr_old_ids, esr_freqs)
                𝕜ᵢ[res_id] = lineshape(ω - ω₀, Ω, Γ)
            end
            PL_on = ThreeHalf.pl_simple(𝔹, D, E, S, 𝕜₀, 𝕜ᵢ, ℙ)
            PL_off = ThreeHalf.pl_simple(𝔹, D, E, S, 𝕜₀, Dict(), ℙ)
            sig[i] += PL_on
            ref[i] += PL_off
        end
    end
    cs .+= (sig .- ref) ./ ref
    plt = plot(ωs, cs, seriestype = :scatter)
    vline!(plt, [2.87e9], linestyle = :dash)
    display(plt)
end

# 3 // 2 :G Δ sweep, <θ>
if false
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
    #     ℙ = Dict((1 // 2) => 1, (-1 // 2) => 1, (3 // 2) => 0, (-3 // 2) => 0)

    function lineshape(Δ, Ω, Γ)
        #         return Ω * exp(-1 * Δ^2 / (2 * Γ^2))
        #         return  Ω * exp(-abs(Δ) / Ω)
        return Ω^2 / (Δ^2 + Γ^2)
    end

    D = 10e6 / 2
    # D = 25e6 / 2
    𝔹norm = 35e-3
    mean_freqs, _, _ = ThreeHalf.calc_esr_freqs([0, 0, 𝔹norm], D, E, S)

    Nθs = 100
    Nωs = 500
    # WIDTH = 5e9
    WIDTH = 2.5e9
    # WIDTH = 1e9
    # WIDTH = 0.5e9
    # WIDTH = 0.25e9
    # ωs = collect(LinRange(mean_freqs[2] - WIDTH/2, mean_freqs[2]+ WIDTH /2, Nωs))
    # ωs = collect(LinRange(1e9 - WIDTH / 2, 1e9 + WIDTH / 2, Nωs))
    ωs = collect(LinRange(0.5e9, 2.5e9, Nωs))
    # ωs = [2e9]

    cs = zeros(Nωs)
    sig = zeros(Nωs)
    ref = zeros(Nωs)
    # ISO 80000-2:2019 convention, θ polar, ϕ azimuthal
    θs = LinRange(0, π, Nθs)
    # θs = LinRange(-π/4, π/4, Nθs)
    # θs = LinRange(-π/5, π/5, Nθs)
    # θs = [0]
    # θs = [pi / 4]
    # θs = [pi/2]
    # θs = [3*pi/5]

    # Γ = 1e7
    Γ = 5e6
    Ω = 1e10
    # Ω = 5e11

    for θ in θs
        # u = ThreeHalf.rad2cart(1, θ, 0)
        # 𝔹 = sum([0, 0, 𝔹norm] .* u) .* [0, 0, 1]
        # 𝔹 = [0, 0, 𝔹norm * u[3]]
        𝔹 = ThreeHalf.rad2cart(𝔹norm, θ, 0)

        esr_freqs, esr_old_ids, _ = ThreeHalf.calc_esr_freqs(𝔹, D, E, S)
        # println(esr_freqs)
        # println(esr_old_ids)

        for (i, ω) in pairs(IndexLinear(), ωs)
            𝕜ᵢ = Dict()
            for (res_id, ω₀) in zip(esr_old_ids, esr_freqs)
                𝕜ᵢ[res_id] = lineshape(ω - ω₀, Ω, Γ)
            end
            PL_on = ThreeHalf.pl_simple(𝔹, D, E, S, 𝕜₀, 𝕜ᵢ, ℙ)
            PL_off = ThreeHalf.pl_simple(𝔹, D, E, S, 𝕜₀, Dict(), ℙ)
            sig[i] += PL_on
            ref[i] += PL_off
        end
    end
    cs .+= (sig .- ref) ./ ref
    plt = plot(ωs, cs, seriestype = :scatter)
    display(plt)
end

# NV :G ℰ vs |𝔹|
if false
    D = 2.87e9
    E = 0
    S = 1 // 1

    NBs = 250

    # ISO 80000-2:2019 convention, θ polar, ϕ azimuthal
    Bnorms = collect(LinRange(0, 0.25, NBs))
    θ = pi / 3

    energy_levels = Array{Float64}(undef, Int(2 * S + 1), NBs)

    for (i, Bn) in enumerate(Bnorms)
        𝔹 = ThreeHalf.rad2cart(Bn, θ, 0)
        ham = ThreeHalf._ham_sd(𝔹, D, E, S)
        evals, _ = eigen(ham)

        energy_levels[:, i] = real(evals)
        # push!(energy_levels, real(evals))
    end
    plt = plot(Bnorms, energy_levels', seriestype = :line)
    display(plt)
end

# 3 // 2 :G ℰ vs |𝔹|
if false
    D = 1e6 / 2
    E = 0
    S = 3 // 2

    NBs = 250

    # ISO 80000-2:2019 convention, θ polar, ϕ azimuthal
    Bnorms = collect(LinRange(0, 0.05, NBs))
    θ = pi / 4
    # θ = 0

    energy_levels = Array{Float64}(undef, Int(2 * S + 1), NBs)

    for (i, Bn) in enumerate(Bnorms)
        𝔹 = ThreeHalf.rad2cart(Bn, θ, 0)
        ham = ThreeHalf._ham_sd(𝔹, D, E, S)
        evals, _ = eigen(ham)

        energy_levels[:, i] = real(evals)
        # push!(energy_levels, real(evals))
    end
    plt = plot(Bnorms, energy_levels', seriestype = :line)
    display(plt)
end

# NV :G ℰ vs θ
if false
    D = 2.87e9
    E = 0
    S = 1 // 1

    Nθs = 50

    # ISO 80000-2:2019 convention, θ polar, ϕ azimuthal
    θs = collect(LinRange(0, π, Nθs))
    Bnorm = 35e-3
    energy_levels = Array{Float64}(undef, Int(2 * S + 1), Nθs)

    for (i, θ) in enumerate(θs)
        𝔹 = ThreeHalf.rad2cart(Bnorm, θ, 0)
        ham = ThreeHalf._ham_sd(𝔹, D, E, S)
        evals, _ = eigen(ham)

        energy_levels[:, i] = real(evals)
    end
    plt = plot(θs, energy_levels', seriestype = :line)
    display(plt)
end

# 3 // 2 :G ℰ vs θ -> boring
if false
    D = 1e6 / 2
    E = 0
    S = 3 // 2

    Nθs = 50

    # ISO 80000-2:2019 convention, θ polar, ϕ azimuthal
    θs = collect(LinRange(0, π, Nθs))
    Bnorm = 35e-3
    energy_levels = Array{Float64}(undef, Int(2 * S + 1), Nθs)

    for (i, θ) in enumerate(θs)
        𝔹 = ThreeHalf.rad2cart(Bnorm, θ, 0)
        ham = ThreeHalf._ham_sd(𝔹, D, E, S)
        evals, _ = eigen(ham)

        energy_levels[:, i] = real(evals)
    end
    plt = plot(θs, energy_levels', seriestype = :line)
    display(plt)
end

# 3 // 2 ms_charac vs θ
if false
    D = 1e6 # 1e6 / 2
    E = 0
    S = 3 // 2

    Nθs = 50

    # ISO 80000-2:2019 convention, θ polar, ϕ azimuthal
    θs = collect(LinRange(0, π, Nθs))
    Bnorm = 35e-3
    characs = Array{Int64}(undef, Int(2 * S + 1), Nθs)
    # characs = Array{Float64}(undef, Int(2*S+1), Nθs)

    for (i, θ) in enumerate(θs)
        𝔹 = ThreeHalf.rad2cart(Bnorm, θ, 0)
        ham = ThreeHalf._ham_sd(𝔹, D, E, S)
        _, evecs = eigen(ham)
        ids = ThreeHalf._define_ms_charac(S, evecs)

        characs[:, i] = ids
    end
    plt = plot(θs, characs', seriestype = :line, lw = 2)
    display(plt)
end

# NV :G Δ sweep, mixed
if false
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

    # simples for now
    function lineshape(Δ, Ω, Γ)
        return Ω^2 / (Δ^2 + Γ^2)
    end

    NVALS = 250
    𝐶s = Vector{Float64}(undef, NVALS)
    # width = 0.05e9
    width = 0.1e9
    ωs = collect(LinRange(2.87e9 - width / 2, 2.87e9 + width / 2, NVALS))

    Γ = 1.5e6
    Ω = 1e9

    for (i, ω) in pairs(IndexLinear(), ωs)
        PL_on = ThreeHalf.pl_mw(𝔹, D, E, S, 𝕜₀, ℙ, ω, lineshape, (Ω, Γ))
        PL_off = ThreeHalf.pl_mw(𝔹, D, E, S, 𝕜₀, ℙ, ω, lineshape, (0, Γ))

        𝐶s[i] = (PL_on - PL_off) / (PL_off == 0 ? 1 : PL_off)
    end
    plt = plot(ωs, 𝐶s, seriestype = :scatter)
    vline!(plt, [2.87e9], linestyle = :dash)
    display(plt)
end

# 3 // 2 :G Δ sweep, mixed
if false
    𝔹 = ThreeHalf.rad2cart(10e-3, pi / 4, 0)
    # 𝔹 = ThreeHalf.rad2cart(10e-3, 0, 0)
    D = 10e6 / 2
    # D = 1e9 / 2
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
    # ℙ = Dict((1 // 2) => 1, (-1 // 2) => 1, (3 // 2) => 0, (-3 // 2) => 0)

    function lineshape(Δ, Ω, Γ)
        #         return Ω * exp(-1 * Δ^2 / (2 * Γ^2))
        return Ω^2 / (Δ^2 + Γ^2)
        #         return  Ω * exp(-abs(Δ) / Ω)
    end

    esr_freqs, esr_old_ids, _ = ThreeHalf.calc_esr_freqs(𝔹, D, E, S)
    print("esr_freqs: ")
    println(esr_freqs)
    print("esr_old_ids: ")
    println(esr_old_ids)
    print("<esr_freqs>: ")
    println(sum(esr_freqs) / length(esr_freqs))

    NVALS = 1000
    # WIDTH = 2e9
    #         WIDTH = 1e9
    # WIDTH = 0.05e9
    WIDTH = 0.5e9
    ωs = collect(
        LinRange(
            sum(esr_freqs) / length(esr_freqs) - WIDTH,
            sum(esr_freqs) / length(esr_freqs) + WIDTH,
            NVALS,
        ),
    )
    𝐶s = Vector{Float64}(undef, NVALS)

    #         Γ = 0
    #     Γ = 1e3
    #  Γ = 5e6
    Γ = 1e7
    #     Γ = 1.5e7
    Ω = 1e9
    #     Ω = 1e7
    #     Ω = 1e11
    for (i, ω) in pairs(IndexLinear(), ωs)
        PL_on = ThreeHalf.pl_mw(𝔹, D, E, S, 𝕜₀, ℙ, ω, lineshape, (Ω, Γ))
        PL_off = ThreeHalf.pl_mw(𝔹, D, E, S, 𝕜₀, ℙ, ω, lineshape, (0, Γ))

        𝐶s[i] = (PL_on - PL_off) / (PL_off == 0 ? 1 : PL_off)
    end
    plt = plot(ωs, 𝐶s, seriestype = :scatter)
    _, idx_min = findmin(𝐶s)
    display(plt)
end

# NV :G Δ sweep, <θ>, mixed
if false
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

    function lineshape(Δ, Ω, Γ)
        #         return Ω * exp(-1 * Δ^2 / (2 * Γ^2))
        #         return  Ω * exp(-abs(Δ) / Ω)
        return Ω^2 / (Δ^2 + Γ^2)
    end

    𝔹norm = 20e-4

    Nθs = 50
    Nωs = 250
    # WIDTH = 5e9
    # WIDTH = 1e9
    WIDTH = 0.2e9
    # WIDTH = 0.1e9
    ωs = collect(LinRange(D - WIDTH / 2, D + WIDTH / 2, Nωs))
    # ωs = collect(LinRange(2.92e9, 2.93e9, Nωs))
    # ωs = range(start = 0, stop = 6e9, length = Nωs)

    cs = zeros(Nωs)
    sig = zeros(Nωs)
    ref = zeros(Nωs)
    # ISO 80000-2:2019 convention, θ polar, ϕ azimuthal
    θs = LinRange(0, π, Nθs)
    # θs = LinRange(-π/4, π/4, Nθs)
    # θs = [0, pi/2]
    # θs = [0]

    Γ = 1e6
    Ω = 1e10
    # Ω = 5e11

    for θ in θs
        # u = ThreeHalf.rad2cart(1, θ, 0)
        # 𝔹 = sum([0, 0, 𝔹norm] .* u) .* [0, 0, 1]
        # 𝔹 = [0, 0, 𝔹norm * u[3]]
        𝔹 = ThreeHalf.rad2cart(𝔹norm, θ, 0)

        esr_freqs, esr_old_ids, _ = ThreeHalf.calc_esr_freqs(𝔹, D, E, S)

        for (i, ω) in pairs(IndexLinear(), ωs)
            PL_on = ThreeHalf.pl_mw(𝔹, D, E, S, 𝕜₀, ℙ, ω, lineshape, (Ω, Γ))
            PL_off = ThreeHalf.pl_mw(𝔹, D, E, S, 𝕜₀, ℙ, ω, lineshape, (0, Γ))

            sig[i] += PL_on
            ref[i] += PL_off
        end
    end
    cs .+= (sig .- ref) ./ ref
    plt = plot(ωs, cs, seriestype = :scatter, legend_position = :bottomright)
    # vline!(plt, [2.87e9], linestyle = :dash)
    display(plt)
end

# 3 // 2 :G Δ sweep, <θ>, mixed
if true
    E = 0
    S = 3 // 2
    k1 = 1e7
    k2 = 1e9
    𝕜₀ = Dict(
        (1 // 2, 3 // 2) => k1,
        (1 // 2, -3 // 2) => k1,
        (-1 // 2, 3 // 2) => k1,
        (-1 // 2, -3 // 2) => k1,
               
        (3 // 2, 1 // 2) => k2,
        (-3 // 2, 1 // 2) => k2,
        (3 // 2, -1 // 2) => k2,
        (-3 // 2, -1 // 2) => k2,
    )
    # ℙ = Dict((1 // 2) => 0, (-1 // 2) => 0, (3 // 2) => 1, (-3 // 2) => 1)
        ℙ = Dict((1 // 2) => 1, (-1 // 2) => 1, (3 // 2) => 0, (-3 // 2) => 0)

    function lineshape(Δ, Ω, Γ)
        #         return Ω * exp(-1 * Δ^2 / (2 * Γ^2))
        #         return  Ω * exp(-abs(Δ) / Ω)
        return Ω^2 / (Δ^2 + Γ^2)
    end

    D = 10e6 / 2
    # D = 25e6 / 2
    𝔹norm = 35e-3
    mean_freqs, _, _ =
        ThreeHalf.calc_esr_freqs([0, 0, 𝔹norm], D, E, S, selrules = false)

    Nθs = 100
    Nωs = 250
    # WIDTH = 5e9
    # WIDTH = 2.5e9
    # WIDTH = 1e9
    # WIDTH = 0.5e9
    WIDTH = 0.25e9
    # ωs = collect(LinRange(mean_freqs[2] - WIDTH/2, mean_freqs[2]+ WIDTH /2, Nωs))
    ωs = collect(LinRange(1e9 - WIDTH / 2, 1e9 + WIDTH / 2, Nωs))

    # ωs = collect(LinRange(0.5e9, 2.5e9, Nωs))

    cs = zeros(Nωs)
    sig = zeros(Nωs)
    ref = zeros(Nωs)
    # ISO 80000-2:2019 convention, θ polar, ϕ azimuthal
    θs = LinRange(0, π, Nθs)
    # θs = LinRange(-π/4, π/4, Nθs)
    # θs = LinRange(-π/5, π/5, Nθs)
    # θs = [0]
    # θs = [pi / 4]
    # θs = [pi/2]
    # θs = [3*pi/5]

    # Γ = 1e7
    Γ = 5e6
    # Γ = 2e6
    # Γ = 1e6
    Ω = 1e10
    # Ω = 5e11

    for θ in θs
        # u = ThreeHalf.rad2cart(1, θ, 0)
        # 𝔹 = sum([0, 0, 𝔹norm] .* u) .* [0, 0, 1]
        # 𝔹 = [0, 0, 𝔹norm * u[3]]
        𝔹 = ThreeHalf.rad2cart(𝔹norm, θ, 0)

        esr_freqs, esr_old_ids, _ =
            ThreeHalf.calc_esr_freqs(𝔹, D, E, S, selrules = false)

        for (i, ω) in pairs(IndexLinear(), ωs)
            PL_on = ThreeHalf.pl_mw(𝔹, D, E, S, 𝕜₀, ℙ, ω, lineshape, (Ω, Γ))
            PL_off = ThreeHalf.pl_mw(𝔹, D, E, S, 𝕜₀, ℙ, ω, lineshape, (0, Γ))

            sig[i] += PL_on
            ref[i] += PL_off
        end
    end
    cs .+= (sig .- ref) ./ ref
    plt = plot(ωs, cs, seriestype = :scatter, legend_position = :bottomright)
    display(plt)
end

end # module

nothing
