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
    ๐น = (0, 0, 1e-3) # T 
    D = 2.87e9
    E = 0
    S = 1 // 1
    ๐โ = Dict(
        # pump
        (1 // 1, 0 // 1) => 1e9,
        (-1 // 1, 0 // 1) => 1e9,
        # small decay
        (0 // 1, 1 // 1) => 1e3,
        (0 // 1, -1 // 1) => 1e3,
    )
    # 1e20 and 1e9 drive the same transition -> should be in same element?
    ๐แตข = Dict(Set((0 // 1, 1 // 1)) => 1e10)
    โ = Dict((0 // 1) => 1, (1 // 1) => 0, (-1 // 1) => 0)
    PL = ThreeHalf.pl_simple(๐น, D, E, S, ๐โ, ๐แตข, โ)
    pprintln(PL)
    @assert โ(PL, 1.0)
end

# 3 // 2 :G only
if false
    ๐น = (0, 0, 50e-3) # T 
    D = 10e6
    E = 0
    S = 3 // 2
    ๐โ = Dict(
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

    ๐แตข = Dict(
        #         Set((-1 // 2, -3 // 2)) => 1e9,
        Set((1 // 2, 3 // 2)) => 1e9,
        # Set((1 // 2, -3 // 2)) => 1e9,
        #         Set((-1 // 2, 3 // 2)) => 1e9,
    )
    โ = Dict((1 // 2) => 0, (-1 // 2) => 0, (3 // 2) => 1, (-3 // 2) => 1)
    PLon = ThreeHalf.pl_simple(๐น, D, E, S, ๐โ, ๐แตข, โ)
    PLoff = ThreeHalf.pl_simple(๐น, D, E, S, ๐โ, Dict(), โ)
    contrast = 100 * (PLon - PLoff) / PLoff
    printfmtln("๐ถ: {:.2f}% ", contrast)
end

# NV :G ฮ sweep
if false
    ๐น = (0, 1e-3, 1e-3) # T 
    D = 2.87e9
    E = 0
    S = 1 // 1
    ๐โ = Dict(
        # pump
        (1 // 1, 0 // 1) => 1e9,
        (-1 // 1, 0 // 1) => 1e9,
    )
    โ = Dict((0 // 1) => 1, (1 // 1) => 0, (-1 // 1) => 0)

    # simples for now
    function lineshape(ฮ, ฮฉ, ฮ)
        return ฮฉ^2 / (ฮ^2 + ฮ^2)
    end

    NVALS = 250
    ๐ถs = Vector{Float64}(undef, NVALS)
    ฯs = collect(LinRange(2.87e9 - 0.1e9, 2.87e9 + 0.1e9, NVALS))

    ฮ = 1.5e6
    ฮฉ = 1e9

    esr_freqs, esr_old_ids, _ = ThreeHalf.calc_esr_freqs(๐น, D, E, S)

    for (i, ฯ) in pairs(IndexLinear(), ฯs)
        ๐แตข = Dict()
        for (res_id, ฯโ) in zip(esr_old_ids, esr_freqs)
            ๐แตข[res_id] = lineshape(ฯ - ฯโ, ฮฉ, ฮ)
        end

        PL_on = ThreeHalf.pl_simple(๐น, D, E, S, ๐โ, ๐แตข, โ)
        PL_off = ThreeHalf.pl_simple(๐น, D, E, S, ๐โ, Dict(), โ)
        ๐ถs[i] = (PL_on - PL_off) / (PL_off == 0 ? 1 : PL_off)
    end
    plt = plot(ฯs, ๐ถs, seriestype = :scatter)
    vline!(plt, [2.87e9], linestyle = :dash)
    display(plt)
end

# 3 // 2 :G ฮ sweep
if false
    # ๐น = (0, 0, -1e-3) # T 
    ๐น = (0, 0, -10e-3) # T 
    #     defined as half, to get correct splitting
    #     D = 10e6
    D = 10e6 / 2
    #     D = 1e9 / 2
    E = 0
    S = 3 // 2
    ๐โ = Dict(
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
    โ = Dict((1 // 2) => 0, (-1 // 2) => 0, (3 // 2) => 1, (-3 // 2) => 1)
    # โ = Dict((1 // 2) => 1, (-1 // 2) => 1, (3 // 2) => 0, (-3 // 2) => 0)

    function lineshape(ฮ, ฮฉ, ฮ)
        #         return ฮฉ * exp(-1 * ฮ^2 / (2 * ฮ^2))
        return ฮฉ^2 / (ฮ^2 + ฮ^2)
        #         return  ฮฉ * exp(-abs(ฮ) / ฮฉ)
    end

    esr_freqs, esr_old_ids, _ = ThreeHalf.calc_esr_freqs(๐น, D, E, S)
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
    ฯs = collect(
        LinRange(
            sum(esr_freqs) / length(esr_freqs) - WIDTH,
            sum(esr_freqs) / length(esr_freqs) + WIDTH,
            NVALS,
        ),
    )
    ๐ถs = Vector{Float64}(undef, NVALS)

    #         ฮ = 0
    #     ฮ = 1e3
    #  ฮ = 5e6
    ฮ = 1e7
    #     ฮ = 1.5e7
    ฮฉ = 1e9
    #     ฮฉ = 1e7
    #     ฮฉ = 1e11
    for (i, ฯ) in pairs(IndexLinear(), ฯs)
        ๐แตข = Dict()

        for (res_id, ฯโ) in zip(esr_old_ids, esr_freqs)
            ๐แตข[res_id] = lineshape(ฯ - ฯโ, ฮฉ, ฮ)
        end

        PL_on = ThreeHalf.pl_simple(๐น, D, E, S, ๐โ, ๐แตข, โ)
        PL_off = ThreeHalf.pl_simple(๐น, D, E, S, ๐โ, Dict(), โ)
        ๐ถs[i] = (PL_on - PL_off) / (PL_off == 0 ? 1 : PL_off)
    end
    plt = plot(ฯs, ๐ถs, seriestype = :scatter)
    _, idx_min = findmin(๐ถs)
    display(plt)
end

# NV :G ฮ sweep, <ฮธ>
if false
    D = 2.87e9
    E = 0
    S = 1 // 1
    ๐โ = Dict(
        # pump
        (1 // 1, 0 // 1) => 1e9,
        (-1 // 1, 0 // 1) => 1e9,
        # small decay
        (0 // 1, 1 // 1) => 1e3,
        (0 // 1, -1 // 1) => 1e3,
    )
    โ = Dict((0 // 1) => 1, (1 // 1) => 0, (-1 // 1) => 0)

    function lineshape(ฮ, ฮฉ, ฮ)
        #         return ฮฉ * exp(-1 * ฮ^2 / (2 * ฮ^2))
        #         return  ฮฉ * exp(-abs(ฮ) / ฮฉ)
        return ฮฉ^2 / (ฮ^2 + ฮ^2)
    end

    ๐นnorm = 20e-4

    Nฮธs = 50
    Nฯs = 250
    # WIDTH = 5e9
    # WIDTH = 1e9
    WIDTH = 0.2e9
    ฯs = collect(LinRange(D - WIDTH / 2, D + WIDTH / 2, Nฯs))

    cs = zeros(Nฯs)
    sig = zeros(Nฯs)
    ref = zeros(Nฯs)
    # ISO 80000-2:2019 convention, ฮธ polar, ฯ azimuthal
    ฮธs = LinRange(0, ฯ, Nฮธs)
    # ฮธs = LinRange(-ฯ/4, ฯ/4, Nฮธs)
    # ฮธs = [0, pi/2]

    ฮ = 1e6
    ฮฉ = 1e10
    # ฮฉ = 5e11

    for ฮธ in ฮธs
        # u = ThreeHalf.rad2cart(1, ฮธ, 0)
        # ๐น = sum([0, 0, ๐นnorm] .* u) .* u
        # ๐น = ๐นnorm .* u
        # ๐น = sum([0, 0, ๐นnorm] .* u) .* [0, 0, 1]
        ๐น = ThreeHalf.rad2cart(๐นnorm, ฮธ, 0)
        # FIXME cleanup?

        esr_freqs, esr_old_ids, _ = ThreeHalf.calc_esr_freqs(๐น, D, E, S)

        for (i, ฯ) in pairs(IndexLinear(), ฯs)
            ๐แตข = Dict()
            for (res_id, ฯโ) in zip(esr_old_ids, esr_freqs)
                ๐แตข[res_id] = lineshape(ฯ - ฯโ, ฮฉ, ฮ)
            end
            PL_on = ThreeHalf.pl_simple(๐น, D, E, S, ๐โ, ๐แตข, โ)
            PL_off = ThreeHalf.pl_simple(๐น, D, E, S, ๐โ, Dict(), โ)
            sig[i] += PL_on
            ref[i] += PL_off
        end
    end
    cs .+= (sig .- ref) ./ ref
    plt = plot(ฯs, cs, seriestype = :scatter)
    vline!(plt, [2.87e9], linestyle = :dash)
    display(plt)
end

# NV :G ฮ sweep, <๐น>
if false
    D = 2.87e9
    E = 0
    S = 1 // 1
    ๐โ = Dict(
        # pump
        (1 // 1, 0 // 1) => 1e9,
        (-1 // 1, 0 // 1) => 1e9,
        # small decay
        (0 // 1, 1 // 1) => 1e3,
        (0 // 1, -1 // 1) => 1e3,
    )
    โ = Dict((0 // 1) => 1, (1 // 1) => 0, (-1 // 1) => 0)

    function lineshape(ฮ, ฮฉ, ฮ)
        #         return ฮฉ * exp(-1 * ฮ^2 / (2 * ฮ^2))
        #         return  ฮฉ * exp(-abs(ฮ) / ฮฉ)
        return ฮฉ^2 / (ฮ^2 + ฮ^2)
    end

    ๐นnorm = 20e-4

    NBs = 250
    Nฯs = 250
    # WIDTH = 5e9
    # WIDTH = 1e9
    WIDTH = 0.2e9
    ฯs = collect(LinRange(D - WIDTH / 2, D + WIDTH / 2, Nฯs))

    cs = zeros(Nฯs)
    sig = zeros(Nฯs)
    ref = zeros(Nฯs)
    # ISO 80000-2:2019 convention, ฮธ polar, ฯ azimuthal
    Random.seed!(1234)
    d = Distributions.Normal(0, 20e-4 / 3)
    # Bnorms = LinRange(0, 20e-4, NBs)
    Bnorms = Random.rand(d, NBs)
    ฮธ = 0

    ฮ = 1e6
    ฮฉ = 1e10
    # ฮฉ = 5e11

    for Bnorm in Bnorms
        ๐น = ThreeHalf.rad2cart(Bnorm, ฮธ, 0)

        esr_freqs, esr_old_ids, _ = ThreeHalf.calc_esr_freqs(๐น, D, E, S)

        for (i, ฯ) in pairs(IndexLinear(), ฯs)
            ๐แตข = Dict()
            for (res_id, ฯโ) in zip(esr_old_ids, esr_freqs)
                ๐แตข[res_id] = lineshape(ฯ - ฯโ, ฮฉ, ฮ)
            end
            PL_on = ThreeHalf.pl_simple(๐น, D, E, S, ๐โ, ๐แตข, โ)
            PL_off = ThreeHalf.pl_simple(๐น, D, E, S, ๐โ, Dict(), โ)
            sig[i] += PL_on
            ref[i] += PL_off
        end
    end
    cs .+= (sig .- ref) ./ ref
    plt = plot(ฯs, cs, seriestype = :scatter)
    vline!(plt, [2.87e9], linestyle = :dash)
    display(plt)
end

# 3 // 2 :G ฮ sweep, <ฮธ>
if false
    E = 0
    S = 3 // 2
    ๐โ = Dict(
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
    โ = Dict((1 // 2) => 0, (-1 // 2) => 0, (3 // 2) => 1, (-3 // 2) => 1)
    #     โ = Dict((1 // 2) => 1, (-1 // 2) => 1, (3 // 2) => 0, (-3 // 2) => 0)

    function lineshape(ฮ, ฮฉ, ฮ)
        #         return ฮฉ * exp(-1 * ฮ^2 / (2 * ฮ^2))
        #         return  ฮฉ * exp(-abs(ฮ) / ฮฉ)
        return ฮฉ^2 / (ฮ^2 + ฮ^2)
    end

    D = 10e6 / 2
    # D = 25e6 / 2
    ๐นnorm = 35e-3
    mean_freqs, _, _ = ThreeHalf.calc_esr_freqs([0, 0, ๐นnorm], D, E, S)

    Nฮธs = 100
    Nฯs = 500
    # WIDTH = 5e9
    WIDTH = 2.5e9
    # WIDTH = 1e9
    # WIDTH = 0.5e9
    # WIDTH = 0.25e9
    # ฯs = collect(LinRange(mean_freqs[2] - WIDTH/2, mean_freqs[2]+ WIDTH /2, Nฯs))
    # ฯs = collect(LinRange(1e9 - WIDTH / 2, 1e9 + WIDTH / 2, Nฯs))
    ฯs = collect(LinRange(0.5e9, 2.5e9, Nฯs))
    # ฯs = [2e9]

    cs = zeros(Nฯs)
    sig = zeros(Nฯs)
    ref = zeros(Nฯs)
    # ISO 80000-2:2019 convention, ฮธ polar, ฯ azimuthal
    ฮธs = LinRange(0, ฯ, Nฮธs)
    # ฮธs = LinRange(-ฯ/4, ฯ/4, Nฮธs)
    # ฮธs = LinRange(-ฯ/5, ฯ/5, Nฮธs)
    # ฮธs = [0]
    # ฮธs = [pi / 4]
    # ฮธs = [pi/2]
    # ฮธs = [3*pi/5]

    # ฮ = 1e7
    ฮ = 5e6
    ฮฉ = 1e10
    # ฮฉ = 5e11

    for ฮธ in ฮธs
        # u = ThreeHalf.rad2cart(1, ฮธ, 0)
        # ๐น = sum([0, 0, ๐นnorm] .* u) .* [0, 0, 1]
        # ๐น = [0, 0, ๐นnorm * u[3]]
        ๐น = ThreeHalf.rad2cart(๐นnorm, ฮธ, 0)

        esr_freqs, esr_old_ids, _ = ThreeHalf.calc_esr_freqs(๐น, D, E, S)
        # println(esr_freqs)
        # println(esr_old_ids)

        for (i, ฯ) in pairs(IndexLinear(), ฯs)
            ๐แตข = Dict()
            for (res_id, ฯโ) in zip(esr_old_ids, esr_freqs)
                ๐แตข[res_id] = lineshape(ฯ - ฯโ, ฮฉ, ฮ)
            end
            PL_on = ThreeHalf.pl_simple(๐น, D, E, S, ๐โ, ๐แตข, โ)
            PL_off = ThreeHalf.pl_simple(๐น, D, E, S, ๐โ, Dict(), โ)
            sig[i] += PL_on
            ref[i] += PL_off
        end
    end
    cs .+= (sig .- ref) ./ ref
    plt = plot(ฯs, cs, seriestype = :scatter)
    display(plt)
end

# NV :G โฐ vs |๐น|
if false
    D = 2.87e9
    E = 0
    S = 1 // 1

    NBs = 250

    # ISO 80000-2:2019 convention, ฮธ polar, ฯ azimuthal
    Bnorms = collect(LinRange(0, 0.25, NBs))
    ฮธ = pi / 3

    energy_levels = Array{Float64}(undef, Int(2 * S + 1), NBs)

    for (i, Bn) in enumerate(Bnorms)
        ๐น = ThreeHalf.rad2cart(Bn, ฮธ, 0)
        ham = ThreeHalf._ham_sd(๐น, D, E, S)
        evals, _ = eigen(ham)

        energy_levels[:, i] = real(evals)
        # push!(energy_levels, real(evals))
    end
    plt = plot(Bnorms, energy_levels', seriestype = :line)
    display(plt)
end

# 3 // 2 :G โฐ vs |๐น|
if false
    D = 1e6 / 2
    E = 0
    S = 3 // 2

    NBs = 250

    # ISO 80000-2:2019 convention, ฮธ polar, ฯ azimuthal
    Bnorms = collect(LinRange(0, 0.05, NBs))
    ฮธ = pi / 4
    # ฮธ = 0

    energy_levels = Array{Float64}(undef, Int(2 * S + 1), NBs)

    for (i, Bn) in enumerate(Bnorms)
        ๐น = ThreeHalf.rad2cart(Bn, ฮธ, 0)
        ham = ThreeHalf._ham_sd(๐น, D, E, S)
        evals, _ = eigen(ham)

        energy_levels[:, i] = real(evals)
        # push!(energy_levels, real(evals))
    end
    plt = plot(Bnorms, energy_levels', seriestype = :line)
    display(plt)
end

# NV :G โฐ vs ฮธ
if false
    D = 2.87e9
    E = 0
    S = 1 // 1

    Nฮธs = 50

    # ISO 80000-2:2019 convention, ฮธ polar, ฯ azimuthal
    ฮธs = collect(LinRange(0, ฯ, Nฮธs))
    Bnorm = 35e-3
    energy_levels = Array{Float64}(undef, Int(2 * S + 1), Nฮธs)

    for (i, ฮธ) in enumerate(ฮธs)
        ๐น = ThreeHalf.rad2cart(Bnorm, ฮธ, 0)
        ham = ThreeHalf._ham_sd(๐น, D, E, S)
        evals, _ = eigen(ham)

        energy_levels[:, i] = real(evals)
    end
    plt = plot(ฮธs, energy_levels', seriestype = :line)
    display(plt)
end

# 3 // 2 :G โฐ vs ฮธ -> boring
if false
    D = 1e6 / 2
    E = 0
    S = 3 // 2

    Nฮธs = 50

    # ISO 80000-2:2019 convention, ฮธ polar, ฯ azimuthal
    ฮธs = collect(LinRange(0, ฯ, Nฮธs))
    Bnorm = 35e-3
    energy_levels = Array{Float64}(undef, Int(2 * S + 1), Nฮธs)

    for (i, ฮธ) in enumerate(ฮธs)
        ๐น = ThreeHalf.rad2cart(Bnorm, ฮธ, 0)
        ham = ThreeHalf._ham_sd(๐น, D, E, S)
        evals, _ = eigen(ham)

        energy_levels[:, i] = real(evals)
    end
    plt = plot(ฮธs, energy_levels', seriestype = :line)
    display(plt)
end

# 3 // 2 ms_charac vs ฮธ
if false
    D = 1e6 # 1e6 / 2
    E = 0
    S = 3 // 2

    Nฮธs = 50

    # ISO 80000-2:2019 convention, ฮธ polar, ฯ azimuthal
    ฮธs = collect(LinRange(0, ฯ, Nฮธs))
    Bnorm = 35e-3
    characs = Array{Int64}(undef, Int(2 * S + 1), Nฮธs)
    # characs = Array{Float64}(undef, Int(2*S+1), Nฮธs)

    for (i, ฮธ) in enumerate(ฮธs)
        ๐น = ThreeHalf.rad2cart(Bnorm, ฮธ, 0)
        ham = ThreeHalf._ham_sd(๐น, D, E, S)
        _, evecs = eigen(ham)
        ids = ThreeHalf._define_ms_charac(S, evecs)

        characs[:, i] = ids
    end
    plt = plot(ฮธs, characs', seriestype = :line, lw = 2)
    display(plt)
end

# NV :G ฮ sweep, mixed
if false
    ๐น = (0, 0, 1e-3) # T 
    D = 2.87e9
    E = 0
    S = 1 // 1
    ๐โ = Dict(
        # pump
        (1 // 1, 0 // 1) => 1e9,
        (-1 // 1, 0 // 1) => 1e9,
    )
    โ = Dict((0 // 1) => 1, (1 // 1) => 0, (-1 // 1) => 0)

    # simples for now
    function lineshape(ฮ, ฮฉ, ฮ)
        return ฮฉ^2 / (ฮ^2 + ฮ^2)
    end

    NVALS = 250
    ๐ถs = Vector{Float64}(undef, NVALS)
    # width = 0.05e9
    width = 0.1e9
    ฯs = collect(LinRange(2.87e9 - width / 2, 2.87e9 + width / 2, NVALS))

    ฮ = 1.5e6
    ฮฉ = 1e9

    for (i, ฯ) in pairs(IndexLinear(), ฯs)
        PL_on = ThreeHalf.pl_mw(๐น, D, E, S, ๐โ, โ, ฯ, lineshape, (ฮฉ, ฮ))
        PL_off = ThreeHalf.pl_mw(๐น, D, E, S, ๐โ, โ, ฯ, lineshape, (0, ฮ))

        ๐ถs[i] = (PL_on - PL_off) / (PL_off == 0 ? 1 : PL_off)
    end
    plt = plot(ฯs, ๐ถs, seriestype = :scatter)
    vline!(plt, [2.87e9], linestyle = :dash)
    display(plt)
end

# 3 // 2 :G ฮ sweep, mixed
if false
    ๐น = ThreeHalf.rad2cart(10e-3, pi / 4, 0)
    # ๐น = ThreeHalf.rad2cart(10e-3, 0, 0)
    D = 10e6 / 2
    # D = 1e9 / 2
    E = 0
    S = 3 // 2
    ๐โ = Dict(
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
    โ = Dict((1 // 2) => 0, (-1 // 2) => 0, (3 // 2) => 1, (-3 // 2) => 1)
    # โ = Dict((1 // 2) => 1, (-1 // 2) => 1, (3 // 2) => 0, (-3 // 2) => 0)

    function lineshape(ฮ, ฮฉ, ฮ)
        #         return ฮฉ * exp(-1 * ฮ^2 / (2 * ฮ^2))
        return ฮฉ^2 / (ฮ^2 + ฮ^2)
        #         return  ฮฉ * exp(-abs(ฮ) / ฮฉ)
    end

    esr_freqs, esr_old_ids, _ = ThreeHalf.calc_esr_freqs(๐น, D, E, S)
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
    ฯs = collect(
        LinRange(
            sum(esr_freqs) / length(esr_freqs) - WIDTH,
            sum(esr_freqs) / length(esr_freqs) + WIDTH,
            NVALS,
        ),
    )
    ๐ถs = Vector{Float64}(undef, NVALS)

    #         ฮ = 0
    #     ฮ = 1e3
    #  ฮ = 5e6
    ฮ = 1e7
    #     ฮ = 1.5e7
    ฮฉ = 1e9
    #     ฮฉ = 1e7
    #     ฮฉ = 1e11
    for (i, ฯ) in pairs(IndexLinear(), ฯs)
        PL_on = ThreeHalf.pl_mw(๐น, D, E, S, ๐โ, โ, ฯ, lineshape, (ฮฉ, ฮ))
        PL_off = ThreeHalf.pl_mw(๐น, D, E, S, ๐โ, โ, ฯ, lineshape, (0, ฮ))

        ๐ถs[i] = (PL_on - PL_off) / (PL_off == 0 ? 1 : PL_off)
    end
    plt = plot(ฯs, ๐ถs, seriestype = :scatter)
    _, idx_min = findmin(๐ถs)
    display(plt)
end

# NV :G ฮ sweep, <ฮธ>, mixed
if false
    D = 2.87e9
    E = 0
    S = 1 // 1
    ๐โ = Dict(
        # pump
        (1 // 1, 0 // 1) => 1e9,
        (-1 // 1, 0 // 1) => 1e9,
        # small decay
        (0 // 1, 1 // 1) => 1e3,
        (0 // 1, -1 // 1) => 1e3,
    )
    โ = Dict((0 // 1) => 1, (1 // 1) => 0, (-1 // 1) => 0)

    function lineshape(ฮ, ฮฉ, ฮ)
        #         return ฮฉ * exp(-1 * ฮ^2 / (2 * ฮ^2))
        #         return  ฮฉ * exp(-abs(ฮ) / ฮฉ)
        return ฮฉ^2 / (ฮ^2 + ฮ^2)
    end

    ๐นnorm = 20e-4

    Nฮธs = 50
    Nฯs = 250
    # WIDTH = 5e9
    # WIDTH = 1e9
    WIDTH = 0.2e9
    # WIDTH = 0.1e9
    ฯs = collect(LinRange(D - WIDTH / 2, D + WIDTH / 2, Nฯs))
    # ฯs = collect(LinRange(2.92e9, 2.93e9, Nฯs))
    # ฯs = range(start = 0, stop = 6e9, length = Nฯs)

    cs = zeros(Nฯs)
    sig = zeros(Nฯs)
    ref = zeros(Nฯs)
    # ISO 80000-2:2019 convention, ฮธ polar, ฯ azimuthal
    ฮธs = LinRange(0, ฯ, Nฮธs)
    # ฮธs = LinRange(-ฯ/4, ฯ/4, Nฮธs)
    # ฮธs = [0, pi/2]
    # ฮธs = [0]

    ฮ = 1e6
    ฮฉ = 1e10
    # ฮฉ = 5e11

    for ฮธ in ฮธs
        # u = ThreeHalf.rad2cart(1, ฮธ, 0)
        # ๐น = sum([0, 0, ๐นnorm] .* u) .* [0, 0, 1]
        # ๐น = [0, 0, ๐นnorm * u[3]]
        ๐น = ThreeHalf.rad2cart(๐นnorm, ฮธ, 0)

        esr_freqs, esr_old_ids, _ = ThreeHalf.calc_esr_freqs(๐น, D, E, S)

        for (i, ฯ) in pairs(IndexLinear(), ฯs)
            PL_on = ThreeHalf.pl_mw(๐น, D, E, S, ๐โ, โ, ฯ, lineshape, (ฮฉ, ฮ))
            PL_off = ThreeHalf.pl_mw(๐น, D, E, S, ๐โ, โ, ฯ, lineshape, (0, ฮ))

            sig[i] += PL_on
            ref[i] += PL_off
        end
    end
    cs .+= (sig .- ref) ./ ref
    plt = plot(ฯs, cs, seriestype = :scatter, legend_position = :bottomright)
    # vline!(plt, [2.87e9], linestyle = :dash)
    display(plt)
end

# 3 // 2 :G ฮ sweep, <ฮธ>, mixed
if true
    E = 0
    S = 3 // 2
    k1 = 1e7
    k2 = 1e9
    ๐โ = Dict(
        (1 // 2, 3 // 2) => k1,
        (1 // 2, -3 // 2) => k1,
        (-1 // 2, 3 // 2) => k1,
        (-1 // 2, -3 // 2) => k1,
               
        (3 // 2, 1 // 2) => k2,
        (-3 // 2, 1 // 2) => k2,
        (3 // 2, -1 // 2) => k2,
        (-3 // 2, -1 // 2) => k2,
    )
    # โ = Dict((1 // 2) => 0, (-1 // 2) => 0, (3 // 2) => 1, (-3 // 2) => 1)
        โ = Dict((1 // 2) => 1, (-1 // 2) => 1, (3 // 2) => 0, (-3 // 2) => 0)

    function lineshape(ฮ, ฮฉ, ฮ)
        #         return ฮฉ * exp(-1 * ฮ^2 / (2 * ฮ^2))
        #         return  ฮฉ * exp(-abs(ฮ) / ฮฉ)
        return ฮฉ^2 / (ฮ^2 + ฮ^2)
    end

    D = 10e6 / 2
    # D = 25e6 / 2
    ๐นnorm = 35e-3
    mean_freqs, _, _ =
        ThreeHalf.calc_esr_freqs([0, 0, ๐นnorm], D, E, S, selrules = false)

    Nฮธs = 100
    Nฯs = 250
    # WIDTH = 5e9
    # WIDTH = 2.5e9
    # WIDTH = 1e9
    # WIDTH = 0.5e9
    WIDTH = 0.25e9
    # ฯs = collect(LinRange(mean_freqs[2] - WIDTH/2, mean_freqs[2]+ WIDTH /2, Nฯs))
    ฯs = collect(LinRange(1e9 - WIDTH / 2, 1e9 + WIDTH / 2, Nฯs))

    # ฯs = collect(LinRange(0.5e9, 2.5e9, Nฯs))

    cs = zeros(Nฯs)
    sig = zeros(Nฯs)
    ref = zeros(Nฯs)
    # ISO 80000-2:2019 convention, ฮธ polar, ฯ azimuthal
    ฮธs = LinRange(0, ฯ, Nฮธs)
    # ฮธs = LinRange(-ฯ/4, ฯ/4, Nฮธs)
    # ฮธs = LinRange(-ฯ/5, ฯ/5, Nฮธs)
    # ฮธs = [0]
    # ฮธs = [pi / 4]
    # ฮธs = [pi/2]
    # ฮธs = [3*pi/5]

    # ฮ = 1e7
    ฮ = 5e6
    # ฮ = 2e6
    # ฮ = 1e6
    ฮฉ = 1e10
    # ฮฉ = 5e11

    for ฮธ in ฮธs
        # u = ThreeHalf.rad2cart(1, ฮธ, 0)
        # ๐น = sum([0, 0, ๐นnorm] .* u) .* [0, 0, 1]
        # ๐น = [0, 0, ๐นnorm * u[3]]
        ๐น = ThreeHalf.rad2cart(๐นnorm, ฮธ, 0)

        esr_freqs, esr_old_ids, _ =
            ThreeHalf.calc_esr_freqs(๐น, D, E, S, selrules = false)

        for (i, ฯ) in pairs(IndexLinear(), ฯs)
            PL_on = ThreeHalf.pl_mw(๐น, D, E, S, ๐โ, โ, ฯ, lineshape, (ฮฉ, ฮ))
            PL_off = ThreeHalf.pl_mw(๐น, D, E, S, ๐โ, โ, ฯ, lineshape, (0, ฮ))

            sig[i] += PL_on
            ref[i] += PL_off
        end
    end
    cs .+= (sig .- ref) ./ ref
    plt = plot(ฯs, cs, seriestype = :scatter, legend_position = :bottomright)
    display(plt)
end

end # module

nothing
