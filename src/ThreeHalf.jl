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
const 𝐠 = 2
const 𝛄 = 𝐠 * 14e9 # Hz/T

# https://en.wikipedia.org/wiki/Spin_(physics)#Higher_spins
# basis = (1, 0, -1) etc.
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

function _rad2cart(r, θ, ϕ)
    return r * sin(θ) * cos(ϕ), r * sin(θ) * sin(ϕ), r * cos(θ)
end

# t-dep hamiltonian method here: 
# https://www.nature.com/articles/s41467-019-09429-x    

function _ham_sd(𝔹, D, E, S::Rational)
    sx, sy, sz = spinop(S)
    return D * sz * sz + E * (sx * sx - sy * sy) + 𝛄 * sum(𝔹 .* [sx, sy, sz])
end

# new estates are returned ordered by increasing energy
function _esr_calculator(S, evals, evecs)
    ms_charac = _define_ms_charac(S, evecs)
    new_charac = _define_new_charac()
    s_idxs = collect(1:Int(2 * S + 1))

    esr_freqs = Array{Float64}(undef, 0)
    # vec of pairs identifying old/new estates for above
    esr_old_ids = Array{Set{Rational}}(undef, 0)
    esr_new_ids = Array{Set{Symbol}}(undef, 0)
    for s_idx in eachindex(evals)
        for other_s_idx in s_idxs
            this_energy = evals[s_idx]
            other_energy = evals[other_s_idx]
            this_ms = ms_charac[s_idx] - 1 - S
            other_ms = ms_charac[other_s_idx] - 1 - S
            this_new = new_charac[s_idx]
            other_new = new_charac[other_s_idx]
            if other_s_idx == s_idx || Set((this_new, other_new)) in esr_new_ids
                continue
            else
                if false # simple mode... (debugging)
                    push!(esr_freqs, abs(this_energy - other_energy))
                    push!(esr_old_ids, Set((this_ms, other_ms)))
                    push!(esr_new_ids, Set((this_new, other_new)))
                else
                    # check for angular momentum selection rules
                    # on the 'old' identities only
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

function calc_esr_freqs(𝔹, D, E, S::Rational)
    ham = _ham_sd(𝔹, D, E, S)
    evals, evecs = eigen(ham)

    return _esr_calculator(S, evals, evecs)
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

# new states are ordered by increasing energy
function _define_ms_charac(S, evecs)
    s_idxs = collect(1:Int(2 * S + 1))
    what_old_ms_is_this_new_s_most_like =
        Array{Union{Nothing, Number}}(nothing, Int(2 * S + 1))
    for old_s in s_idxs # very messy algorithm...
        new_s_projs = abs.(evecs[:, old_s]) # projections
        new_s_cands = collect(1:Int(2 * S + 1)) # candidates (idx in s_idxs)
        while true
            _, this_idx = findmax(new_s_projs)
            new_s = new_s_cands[this_idx]
            if isnothing(what_old_ms_is_this_new_s_most_like[new_s])
                what_old_ms_is_this_new_s_most_like[new_s] = old_s
                break # ok we found 'it', break out to for loop
            else
                # this new_s vector already used, retry after deleting this optn
                deleteat!(new_s_projs, this_idx)
                deleteat!(new_s_cands, this_idx)
            end
        end
    end
    return what_old_ms_is_this_new_s_most_like
end

# ----------------------------------- PL Stuff -------------------------------

# only for gs for now...
function pl_simple(𝔹, D, E, S::Rational, 𝕜₀, 𝕜ᵢ, ℙ)
    ham = _ham_sd(𝔹, D, E, S)
    evals, 𝔸 = eigen(ham) # 𝔸 = evecs

    nstates = Int(2 * S + 1)
    ms_charac = _define_ms_charac(S, 𝔸)
    old_charac = collect(S:-1:(-S))

    # define transition matrix
    k = zeros((nstates, nstates))
    # pumping rates etc. in old basis
    k0 = zeros((nstates, nstates))

    for (tsn, rate) in pairs(𝕜₀)
        first, second = tsn
        idx = findfirst(old_charac .== first)
        jdx = findfirst(old_charac .== second)
        k0[idx, jdx] = rate
    end

    # first add terms due to mixing
    for i in 1:nstates
        for j in 1:nstates
            for p in 1:nstates
                for q in 1:nstates
                    # k0[p, q] is in **old basis**
                    k[i, j] += abs(𝔸[i, p])^2 * abs(𝔸[j, q])^2 * k0[p, q]
                end
            end
        end
    end

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
    
    PL = 0
    for (state, rate) in pairs(ℙ)
        idx = findfirst(ms_charac .== (state + S + 1))
        PL += pops[idx] * rate
    end
    return PL
end

function pl(𝔹, Ds, Es, S::Rational, 𝕜₀, 𝕜ᵢ)
    ham_gs = _ham_sd(𝔹, Ds[1], Es[1], S)
    evals_gs, evecs_gs = eigen(ham_gs)
    ham_es = _ham_sd(𝔹, Ds[2], Es[2], S)
    evals_es, evecs_es = eigen(ham_es)

    #esr_freqs, esr_old_ids, esr_new_ids = _esr_calculator(S, evals, evecs)

    nstates = 2 * Int(2 * S + 1) + 1
    ms_charac_gs = _define_ms_charac(S, evecs_gs)
    ms_charac_es = _define_ms_charac(S, evecs_es)
    old_charac = collect(S:-1:(-S)) # true for gs and es

    𝔸 = zeros((nstates, nstates))
    𝔸[1:Int(2 * S + 1), 1:Int(2 * S + 1)] = evecs_gs
    𝔸[Int(2 * S + 2):(end - 1), Int(2 * S + 2):(end - 1)] = evecs_es
    𝔸[end, end] = 1 # metastable doesn't care about field

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
            ArgumentError("bad 'level' symbol, use :M, :G or :E")
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

    PL = 0
    for ii in Int(2 * S + 2):(2 * Int(2 * S + 1)) # es
        PL += pops[ii] * k_rad[ii]
    end
    return PL
end

end # module
