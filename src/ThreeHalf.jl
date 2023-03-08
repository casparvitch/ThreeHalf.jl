module ThreeHalf

using Plots
using LinearAlgebra

# ---------------------------------- Rabi Stuff -------------------------------
# p1002 in Drake
function rabi(τ, Ω, Δ)
    Ωg = √(Ω^2 + Δ^2)
    return (Ω / Ωg)^2 * sin(Ωg * τ / 2)^2
end

function dt_rabi(τ, Ω, Δ)
    Ωg = √(Ω^2 + Δ^2)
    return (Ω^2 / (2 * Ωg)) * sin(τ * Ωg)
end

function series_plotter(x, ys)
    plt = plot()
    for y in ys
        plot!(plt, x, y)
    end
    return plt
end

function series_plotter(x, ys, labels)
    plt = plot()
    for (label, y) in zip(labels, ys)
        plot!(plt, x, y; label)
    end
    return plt
end

function plot_rabi(x, y)
    return plot(x, y)
end

# ----------------------------------- ESR Stuff -------------------------------
# I use 'ms' where I mean 's' a bunch here
# or we can say 's' is the 'ms' of the B-applied system
# [these two lines don't make sense together? ^ ]
# [I think I mean where I say s_idxs etc.]
const 𝐠 = 2
const 𝛄 = 𝐠 * 14e9 # Hz/T

# https://en.wikipedia.org/wiki/Spin_(physics)#Higher_spins
# basis ordering is mₛ = (1, 0, -1) etc.
function spinop(S::Rational)
    δ(x, y) = ==(x, y)
    σx = Matrix{ComplexF64}(undef, (Int(2 * S + 1), Int(2 * S + 1)))
    σy = Matrix{ComplexF64}(undef, (Int(2 * S + 1), Int(2 * S + 1)))
    σz = Matrix{ComplexF64}(undef, (Int(2 * S + 1), Int(2 * S + 1)))

    @inbounds for i in 1:Int(2 * S + 1)
        for j in 1:Int(2 * S + 1)
            sqrtfac = √(Complex((S + 1) * (i + j - 1) - i * j))
            σx[i, j] = (1 / 2) * (δ(i, j + 1) + δ(i + 1, j)) * sqrtfac
            σy[i, j] = (1im / 2) * (δ(i, j + 1) - δ(i + 1, j)) * sqrtfac
            σz[i, j] = (S + 1 - i) * δ(i, j)
        end
    end
    return [σx, σy, σz]
end

# ISO 80000-2:2019 convention, θ polar, ϕ azimuthal
function rad2cart(r, θ, ϕ)
    return r * sin(θ) * cos(ϕ), r * sin(θ) * sin(ϕ), r * cos(θ)
end

# t-dep hamiltonian method here: 
# https://www.nature.com/articles/s41467-019-09429-x    

function _ham_sd(𝔹, D, E, S::Rational)
    sx, sy, sz = spinop(S)
    # hmmm... any reason to not have a negative before γ there?
    return D * sz * sz + E * (sx * sx - sy * sy) + 𝛄 * sum(𝔹 .* [sx, sy, sz])
end

# new estates are returned ordered by increasing energy
function _esr_calculator(S, evals, evecs; selrules)
    ms_charac = _define_ms_charac(S, evecs)

    new_charac = _define_new_charac()
    s_idxs = collect(1:Int(2 * S + 1))

    esr_freqs = Array{Float64}(undef, 0)
    # vec of pairs identifying old/new estates for above
    esr_old_ids = Array{Set{Rational}}(undef, 0)
    esr_new_ids = Array{Set{Symbol}}(undef, 0)
    for s_idx in s_idxs
        for other_s_idx in s_idxs
            this_energy = evals[s_idx]
            other_energy = evals[other_s_idx]
            this_ms = -1 * (ms_charac[s_idx] - 1 - S)
            other_ms = -1 * (ms_charac[other_s_idx] - 1 - S)
            this_new = new_charac[s_idx]
            other_new = new_charac[other_s_idx]
            if other_s_idx == s_idx || Set((this_new, other_new)) in esr_new_ids
                continue
            else
                if !selrules # simple mode... (no selection rules)
                    push!(esr_freqs, abs(this_energy - other_energy))
                    push!(esr_old_ids, Set((this_ms, other_ms)))
                    push!(esr_new_ids, Set((this_new, other_new)))
                else
                    # check for angular momentum selection rules
                    # on the 'old' identities only
                    # only valid when there's minimal state-mixing 
                    #   (e.g. low off-axis 𝔹)
                    if abs(this_ms - other_ms) != 1
                        continue
                    else
                        push!(esr_freqs, abs(this_energy - other_energy))
                        push!(esr_old_ids, Set((this_ms, other_ms)))
                        push!(esr_new_ids, Set((this_new, other_new)))
                    end
                end
            end
        end
    end
    return esr_freqs, esr_old_ids, esr_new_ids
end

function calc_esr_freqs(𝔹, D, E, S::Rational; selrules::Bool = true)
    ham = _ham_sd(𝔹, D, E, S)
    evals, evecs = eigen(ham)

    return _esr_calculator(S, evals, evecs, selrules = selrules)
end

function _define_new_charac()
    return [
        :α,
        :β,
        :γ,
        :δ,
        :ϵ,
        :ζ,
        :η,
        :θ,
        :ι,
        :κ,
        :λ,
        :μ,
        :ν,
        :ξ,
        :ο,
        :π,
        :ρ,
        :σ,
        :τ,
        :υ,
        :χ,
        :ψ,
        :ω,
    ]
end

function _define_ms_charac(S, evecs)
    s_idxs = collect(1:Int(2 * S + 1))
    # idx in this array represents new_ms, value is idx in old_ms array
    what_old_ms_is_this_new_s_most_like =
        Array{Union{Nothing, Number}}(nothing, Int(2 * S + 1))
    # print("+++ "); println(-1 .* (collect(1:Int(2 * S + 1)) .- S .- 1))

    # very messy algorithm...
    # iterate over new estate indexes, find candidate it is most like (by
    # projection onto evectors). Assume 1:1 old-new estate transform.

    for new_idx in s_idxs
        projs = abs.(evecs[:, new_idx]) # projection onto evec bases
        old_s_candidates = collect(1:Int(2 * S + 1)) # candidates (idx in s_idxs)
        while true
            _, most_similar_ms_idx = findmax(projs)
            old_s = old_s_candidates[most_similar_ms_idx]
            # println((-1 .* (new_idx - S - 1), -1 .* (old_s -S -1), projs, sum(projs.^2)))
            if isnothing(what_old_ms_is_this_new_s_most_like[new_idx]) &&
               !(old_s in what_old_ms_is_this_new_s_most_like)
                what_old_ms_is_this_new_s_most_like[new_idx] = old_s
                break # ok we found 'it', break out *to* for loop
            else
                # println("haha I'm in danger")
                # this new_s basis already used, retry after deleting this optn
                deleteat!(projs, most_similar_ms_idx)
                deleteat!(old_s_candidates, most_similar_ms_idx)
            end
        end
    end
    # print("=== "); println(what_old_ms_is_this_new_s_most_like)
    return what_old_ms_is_this_new_s_most_like
end

# ----------------------------------- PL Stuff -------------------------------

function zero_field_optical_transitions(𝕜₀, S)
    old_charac = collect(S:-1:(-S))
    nstates = Int(2 * S + 1)
    k0 = zeros((nstates, nstates))
    for (tsn, rate) in pairs(𝕜₀)
        first, second = tsn
        idx = findfirst(old_charac .== first)
        jdx = findfirst(old_charac .== second)
        k0[idx, jdx] = rate
    end
    return k0
end

function add_mixed_optical_transitions!(k, k0, 𝔸)
    nstates = size(k, 1)
    for i in 1:nstates
        for j in 1:nstates
            for p in 1:nstates
                for q in 1:nstates
                    # k0[p, q] is in **old basis**
                    # below line fails
                    # k[i, j] += abs(𝔸[i, p])^2 * abs(𝔸[j, q])^2 * k0[p, q]
                    k[i, j] += abs(𝔸[p, i])^2 * abs(𝔸[q, j])^2 * k0[p, q]
                    # k[i, j] += abs(𝔸'[i, p])^2 * abs(𝔸'[j, q])^2 * k0[p, q]
                end
            end
        end
    end
end

function add_unmixed_mw_transtitions!(k, 𝔸, 𝕜ᵢ, S)
    ms_charac = _define_ms_charac(S, 𝔸)
    # set mw transition rates
    for (tsn, rate) in pairs(𝕜ᵢ)
        tsnr = collect(tsn)
        if length(tsnr) == 1
            t = Int(tsnr[1] + S + 1)
            first, second = t, t
        else
            first, second = Int.(tsnr .+ S .+ 1)
        end
        idx = findfirst(ms_charac .== first)
        jdx = findfirst(ms_charac .== second)
        k[idx, jdx] += rate
        k[jdx, idx] += idx != jdx ? rate : 0
    end
end

function add_mixed_mw_transitions!(k, 𝔸, ω, S, evals, evecs, ls, ls_params)
    nstates = size(k, 1)
    selrule(a, b) = abs(a - b) == 1 # Δmₛ = ±1
    esr_freqs, esr_old_ids, esr_new_ids =
        _esr_calculator(S, evals, evecs, selrules = false)
    new_charac = _define_new_charac()
    for (ω₀, res_ids) in zip(esr_freqs, esr_new_ids)
        idar = collect(res_ids)
        from = findfirst(new_charac .== idar[1]) # new basis 'from'
        to = findfirst(new_charac .== idar[2]) # new basis 'to'
        for p in 1:nstates # old basis 'from'
            for q in 1:nstates # old basis 'to'
                old_ms_from = -1 * (p - S - 1)
                old_ms_to = -1 * (q - S - 1)
                # |𝔸|'s: weight by mixing,
                # ls: lineshape (detuning)
                # selrule: selection rule, switch Ω⁰ on/off if allowed/not
                omega =
                    abs(𝔸[p, from])^2 *
                    abs(𝔸[q, to])^2 *
                    ls(ω - ω₀, ls_params...) *
                    selrule(old_ms_from, old_ms_to)
                k[from, to] += omega
                k[to, from] += to != from ? omega : 0
            end
        end
    end
end

function solve_steady_state(k)
    nstates = size(k, 1)
    # define rate equation matrix d(pop)/dt = M*pop
    𝕄 = zeros((nstates, nstates))
    for i in 1:nstates
        for j in 1:nstates
            if i == j
                for q in 1:nstates
                    𝕄[i, i] -= k[i, q]
                end
            else
                𝕄[i, j] = k[j, i]
            end
        end
    end

    # solve rate equation in the steady state (eval ~ 0)
    evals_M, evecs_M = eigen(𝕄)
    _, idx_ss = findmin(abs.(evals_M))

    pops = real(evecs_M[:, idx_ss])
    pops /= sum(pops)
    return pops
end

function pl_from_unmixed_pops(pops, ℙ, S, 𝔸)
    # approximation - take ms this is most-like
    ms_charac = _define_ms_charac(S, 𝔸)
    PL = 0
    for (state, rate) in pairs(ℙ)
        idx = findfirst(ms_charac .== (state + S + 1))
        PL += pops[idx] * rate
    end
    return PL
end

function pl_from_mixed_pops(pops, ℙ, S, 𝔸)
    old_charac = collect(S:-1:(-S))
    nstates = Int(2 * S + 1)
    PL = 0
    for (pidx, pop) in zip(1:nstates, pops)
        for (state, rate) in pairs(ℙ)
            sidx = findfirst(old_charac .== state)
            PL += abs(𝔸[sidx, pidx])^2 * rate * pop
        end
    end
    return PL
end

# only for gs for now...
function pl_simple(𝔹, D, E, S::Rational, 𝕜₀, 𝕜ᵢ, ℙ)
    ham = _ham_sd(𝔹, D, E, S)
    evals, evecs = eigen(ham)
    𝔸 = evecs

    nstates = Int(2 * S + 1)

    k = zeros((nstates, nstates)) # transition matrix
    k0 = zero_field_optical_transitions(𝕜₀, S)

    add_mixed_optical_transitions!(k, k0, 𝔸)

    add_unmixed_mw_transtitions!(k, 𝔸, 𝕜ᵢ, S)
    pops = solve_steady_state(k)

    return pl_from_unmixed_pops(pops, ℙ, S, 𝔸)
end

function pl_mw(𝔹, D, E, S::Rational, 𝕜₀, ℙ, ω, ls, ls_params)
    ham = _ham_sd(𝔹, D, E, S)
    evals, evecs = eigen(ham)
    # 𝔸 = evecs' # unnecessary
    𝔸 = evecs

    nstates = Int(2 * S + 1)

    k = zeros((nstates, nstates)) # transition matrix
    k0 = zero_field_optical_transitions(𝕜₀, S)

    add_mixed_optical_transitions!(k, k0, 𝔸)

    add_mixed_mw_transitions!(k, 𝔸, ω, S, evals, evecs, ls, ls_params)

    pops = solve_steady_state(k)

    return pl_from_mixed_pops(pops, ℙ, S, 𝔸)
end

# NOTE this guy has not been updated for e.g. ±𝔹 bugfixes
# keep it here to track how to do excited state transition stuff
function pl(𝔹, Ds, Es, S::Rational, 𝕜₀, 𝕜ᵢ)
    ham_gs = _ham_sd(𝔹, Ds[1], Es[1], S)
    evals_gs, evecs_gs = eigen(ham_gs)
    ham_es = _ham_sd(𝔹, Ds[2], Es[2], S)
    evals_es, evecs_es = eigen(ham_es)

    nstates = 2 * Int(2 * S + 1) + 1
    ms_charac_gs = _define_ms_charac(S, evecs_gs)
    ms_charac_es = _define_ms_charac(S, evecs_es)
    old_charac = collect(S:-1:(-S)) # true for gs and es

    𝔸 = zeros((nstates, nstates))
    𝔸[1:Int(2 * S + 1), 1:Int(2 * S + 1)] = evecs_gs
    # using lastindex instead of end is nicer on my editor for code folding...
    𝔸[
        Int(2 * S + 2):(lastindex(𝔸, 2) - 1),
        Int(2 * S + 2):(lastindex(𝔸, 2) - 1),
    ] = evecs_es
    𝔸[lastindex(𝔸, 1), lastindex(𝔸, 2)] = 1 # metastable doesn't care about field

    # define transition matrix, [gs..., es..., singlet]
    k = zeros((nstates, nstates))
    # initial rates, e.g. pumping, in old basis
    k0 = zeros((nstates, nstates))
    for (tsn, rate) in pairs(𝕜₀)
        # todo check tsn is in -S:S
        (level1, spin1), (level2, spin2) = tsn
        idx = if level1 == :M
            nstates
        else
            findfirst(old_charac .== spin1) + (level1 == :E ? Int(2 * S + 1) : 0)
        end
        jdx = if level2 == :M
            nstates
        else
            findfirst(old_charac .== spin2) + (level2 == :E ? Int(2 * S + 1) : 0)
        end
        # this k0 needs to be defined in the **old** basis
        k0[idx, jdx] = rate
    end

    # add terms due to mixing
    for i in 1:nstates
        for j in 1:nstates
            for p in 1:nstates
                for q in 1:nstates
                    # NB: k0[p, q] is in **old basis**
                    k[i, j] += abs(𝔸[i, p])^2 * abs(𝔸[j, q])^2 * k0[p, q]
                end
            end
        end
    end

    # set mw transition rates
    for (tsn, rate) in pairs(𝕜ᵢ)
        # todo check tsn is in -S:S
        tsnr = collect(tsn)
        if length(tsnr) == 1
            (level1, spin1), (level2, spin2) = tsnr, tsnr
        else
            (level1, spin1), (level2, spin2) = tsnr
        end
        if level1 == :M
            idx = nstates # prolly shouldn't be using this here anyway...
        elseif level1 == :G
            idx = findfirst(ms_charac_gs .== Int(spin1 + S + 1))
        elseif level1 == :E
            idx = findfirst(ms_charac_es .== Int(spin1 + S + 1))
        else
            throw(ArgumentError("bad 'level' symbol, use :M, :G or :E"))
        end
        if level2 == :M
            jdx = nstates # prolly shouldn't be using this here anyway...
        elseif level2 == :G
            jdx = findfirst(ms_charac_gs .== Int(spin2 + S + 1))
        elseif level2 == :E
            jdx = findfirst(ms_charac_es .== Int(spin2 + S + 1))
        else
            ArgumentError("bad 'level' symbol, use :M, :G or :E")
        end
        # should be a += or an =? (in general a += I believe)
        k[idx, jdx] += rate
        k[jdx, idx] += idx != jdx ? rate : 0 # don't add twice ? unlikely occur.
    end

    # define rate equation matrix d(pop)/dt = M*pop
    𝕄 = zeros((nstates, nstates))
    for i in 1:nstates
        for j in 1:nstates
            if i == j
                for q in 1:nstates
                    𝕄[i, i] -= k[i, q]
                end
            else
                𝕄[i, j] = k[j, i]
            end
        end
    end

    # solve rate equation in the steady state (eval ~ 0)
    evals_M, evecs_M = eigen(𝕄)
    _, idx_ss = findmin(abs.(evals_M))

    pops = real(evecs_M[:, idx_ss])
    pops /= sum(pops)

    k_rad = zeros((nstates))
    for ii in Int(2 * S + 2):(2 * Int(2 * S + 1)) # es
        for jj in 1:Int(2 * S + 1) # gs
            k_rad[ii] += k[ii, jj]
        end
    end

    PL = 0 # sum over emission from excited states
    for ii in Int(2 * S + 2):(2 * Int(2 * S + 1)) # es
        PL += pops[ii] * k_rad[ii]
    end
    return PL
end

end # module
