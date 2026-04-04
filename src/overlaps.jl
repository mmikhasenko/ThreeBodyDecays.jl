"""
    ChainAmplitudeCache{T,S}

Heavy reusable cache of stripped per-event chain amplitudes.

The `amplitudes` array is stored as `(nchains, nhelicity, nsamples)` and contains the
bare chain amplitudes `a_i(σ_k)` without couplings applied.
"""
struct ChainAmplitudeCache{T,S}
    labels::Vector{String}
    samples::Vector{S}
    amplitudes::Array{Complex{T},3}
end

"""
    OverlapContribution{T}

Labeled additive overlap contribution.

This shared additive container can represent one-event contributions, partial sums, full
sums, and sample means. Derived normalized observables such as fit fractions or
interference fractions should be computed only after aggregation.
"""
struct OverlapContribution{T}
    labels::Vector{String}
    matrix::Hermitian{Complex{T},Matrix{Complex{T}}}
end

const OverlapMatrix = OverlapContribution

_recoupling_label(h::NoRecoupling) = "$(d2(h.two_λa)), $(d2(h.two_λb))"
_recoupling_label(h::ParityRecoupling) = "$(d2(h.two_λa)), $(d2(h.two_λb))"
_recoupling_label(h::RecouplingLS) = "l=$(d2(h.two_ls[1])), s=$(d2(h.two_ls[2]))"
_recoupling_label(h) = sprint(show, MIME"text/plain"(), h)

function chain_label(name, chain::DecayChain)
    string(
        name,
        "_{k=",
        chain.k,
        ", ij=(",
        _recoupling_label(chain.Hij.h),
        "), Rk=(",
        _recoupling_label(chain.HRk.h),
        ")}",
    )
end

chain_label(name, ::AbstractDecayChain) = String(name)

function add_event_contribution!(overlaps, amplitudes, k)
    nchains, nhelicity, _ = size(amplitudes)
    @inbounds for j in 1:nchains
        for i in 1:j
            z = zero(eltype(overlaps))
            for h in 1:nhelicity
                z += conj(amplitudes[i, h, k]) * amplitudes[j, h, k]
            end
            overlaps[i, j] += z
        end
    end
    overlaps
end

"""
    chain_amplitudes(model, samples)

Compute and store stripped per-event chain amplitudes for a `ThreeBodyDecay` model.
"""
function chain_amplitudes(model::ThreeBodyDecay, samples)
    labels = Vector{String}(undef, length(model.chains))
    for i in eachindex(model.chains)
        labels[i] = String(chain_label(model.names[i], model.chains[i]))
    end
    collected_samples = collect(samples)
    nchains = length(model.chains)
    nhelicity = length(amplitude(model.chains[1], collected_samples[1]))
    amplitudes = Array{ComplexF64}(undef, nchains, nhelicity, length(collected_samples))

    for (k, σs) in pairs(collected_samples)
        for i in eachindex(model.chains)
            chain_amplitude = amplitude(model.chains[i], σs)
            @inbounds for h in eachindex(chain_amplitude)
                amplitudes[i, h, k] = chain_amplitude[h]
            end
        end
    end

    ChainAmplitudeCache(labels, collected_samples, amplitudes)
end

"""
    chain_overlap_matrix(cache)

Construct the stripped chain-basis overlap contribution from a `ChainAmplitudeCache`.
"""
function chain_overlap_matrix(cache::ChainAmplitudeCache{T}) where {T}
    nchains, _, nsamples = size(cache.amplitudes)
    overlaps = zeros(Complex{T}, nchains, nchains)

    @inbounds for k in 1:nsamples
        add_event_contribution!(overlaps, cache.amplitudes, k)
    end

    OverlapContribution(cache.labels, Hermitian(overlaps, :U))
end

chain_overlap_contribution(cache::ChainAmplitudeCache) = chain_overlap_matrix(cache)

"""
    event_overlap_contributions(cache)

Construct one stripped chain-basis overlap contribution per event.
"""
function event_overlap_contributions(cache::ChainAmplitudeCache{T}) where {T}
    nchains, _, nsamples = size(cache.amplitudes)
    xs = Vector{OverlapContribution{T}}(undef, nsamples)

    @inbounds for k in 1:nsamples
        overlaps = zeros(Complex{T}, nchains, nchains)
        add_event_contribution!(overlaps, cache.amplitudes, k)
        xs[k] = OverlapContribution(cache.labels, Hermitian(overlaps, :U))
    end

    xs
end

event_overlaps(cache::ChainAmplitudeCache) = event_overlap_contributions(cache)

"""
    chain_overlap_matrix(model, samples)

Convenience constructor for the stripped chain-basis overlap contribution.
"""
chain_overlap_matrix(model::ThreeBodyDecay, samples) =
    chain_overlap_matrix(chain_amplitudes(model, samples))
chain_overlap_contribution(model::ThreeBodyDecay, samples) = chain_overlap_matrix(model, samples)

function _check_overlap_compatibility(xs::AbstractVector{<:OverlapContribution})
    isempty(xs) && throw(ArgumentError("Need a non-empty vector of overlap contributions"))
    labels = xs[1].labels
    n = length(labels)
    for x in xs
        x.labels == labels ||
            throw(ArgumentError("All overlap contributions must have identical labels"))
        size(x.matrix) == (n, n) ||
            throw(ArgumentError("All overlap contributions must have identical matrix size"))
    end
    labels
end

function _sum_overlap_contributions(xs::AbstractVector{<:OverlapContribution})
    labels = _check_overlap_compatibility(xs)
    acc = zeros(ComplexF64, length(labels), length(labels))
    for x in xs
        acc .+= Matrix(x.matrix)
    end
    OverlapContribution(labels, Hermitian(acc, :U))
end

function Statistics.mean(xs::AbstractVector{<:OverlapContribution})
    total = _sum_overlap_contributions(xs)
    OverlapContribution(total.labels, Hermitian(Matrix(total.matrix) ./ length(xs), :U))
end

function Statistics.std(xs::AbstractVector{<:OverlapContribution}; corrected::Bool=true, mean=nothing)
    labels = _check_overlap_compatibility(xs)
    μ = isnothing(mean) ? Statistics.mean(xs) : mean
    n = length(labels)
    sre = zeros(Float64, n, n)
    sim = zeros(Float64, n, n)
    scale = corrected ? length(xs) - 1 : length(xs)
    scale > 0 ||
        throw(ArgumentError("Need at least two overlap contributions for corrected standard deviation"))

    for x in xs
        Δ = Matrix(x.matrix) .- Matrix(μ.matrix)
        sre .+= real.(Δ) .^ 2
        sim .+= imag.(Δ) .^ 2
    end

    (labels=labels, re=sqrt.(sre ./ scale), im=sqrt.(sim ./ scale))
end

"""
    stderr(xs)

Elementwise Monte Carlo standard error for a vector of event-level overlap contributions.
"""
function stderr(xs::AbstractVector{<:OverlapContribution}; corrected::Bool=true)
    s = Statistics.std(xs; corrected)
    n = length(xs)
    (labels=s.labels, re=s.re ./ sqrt(n), im=s.im ./ sqrt(n))
end

"""
    physical_overlap(overlap, coefficients)

Apply physical coefficients to a stripped overlap contribution.
"""
function physical_overlap(overlap::OverlapContribution, coefficients)
    length(coefficients) == length(overlap.labels) ||
        throw(DimensionMismatch("Expected $(length(overlap.labels)) coefficients, got $(length(coefficients))"))
    weighted = Matrix(overlap.matrix)
    @inbounds for j in eachindex(coefficients), i in eachindex(coefficients)
        weighted[i, j] = conj(coefficients[i]) * weighted[i, j] * coefficients[j]
    end
    OverlapContribution(overlap.labels, Hermitian(weighted, :U))
end

"""
    group_overlap(overlap, groups)

Coarse-grain an overlap contribution by summing blocks according to `groups`.
"""
function group_overlap(overlap::OverlapContribution, groups)
    length(groups) == length(overlap.labels) ||
        throw(DimensionMismatch("Expected $(length(overlap.labels)) group labels, got $(length(groups))"))

    labels = String[]
    lookup = Dict{String,Int}()
    component_index = Vector{Int}(undef, length(groups))
    for (i, group) in pairs(groups)
        label = String(group)
        idx = get!(lookup, label) do
            push!(labels, label)
            length(labels)
        end
        component_index[i] = idx
    end

    grouped = zeros(ComplexF64, length(labels), length(labels))
    @inbounds for j in eachindex(component_index), i in eachindex(component_index)
        grouped[component_index[i], component_index[j]] += overlap.matrix[i, j]
    end

    OverlapContribution(labels, Hermitian(grouped, :U))
end

"""
    total_intensity(overlap)

Return the total intensity represented by a labeled overlap contribution.
"""
function total_intensity(overlap::OverlapContribution)
    total = 0.0
    @inbounds for j in eachindex(overlap.labels)
        total += real(overlap.matrix[j, j])
        for i in 1:(j - 1)
            total += 2 * real(overlap.matrix[i, j])
        end
    end
    total
end

"""
    diagonal_integrals(overlap)

Return the diagonal contributions of a labeled overlap contribution.
"""
diagonal_integrals(overlap::OverlapContribution) =
    [real(overlap.matrix[i, i]) for i in eachindex(overlap.labels)]

diagonal_contributions(overlap::OverlapContribution) = diagonal_integrals(overlap)

"""
    fit_fractions(overlap; normalize=true, percent=true, sort=true)

Return labeled diagonal contributions of an aggregated overlap contribution.
"""
function fit_fractions(overlap::OverlapContribution; normalize::Bool=true, percent::Bool=true, sort::Bool=true)
    values = diagonal_integrals(overlap)
    scale = percent ? 100.0 : 1.0
    values = normalize ? values ./ total_intensity(overlap) .* scale : values .* scale
    rows = [(; label, fraction) for (label, fraction) in zip(overlap.labels, values)]
    sort ? Base.sort(rows; by=x -> x.fraction, rev=true) : rows
end

fit_fractions(xs::AbstractVector{<:OverlapContribution}; kwargs...) =
    fit_fractions(mean(xs); kwargs...)

"""
    interference_terms(overlap; normalize=true, percent=true, sort=true)

Return labeled pairwise interference contributions from the off-diagonal part of an
aggregated overlap contribution.
"""
function interference_terms(overlap::OverlapContribution; normalize::Bool=true, percent::Bool=true, sort::Bool=true)
    total = normalize ? total_intensity(overlap) : 1.0
    scale = percent ? 100.0 : 1.0
    rows = NamedTuple[]

    @inbounds for j in 2:length(overlap.labels)
        for i in 1:(j - 1)
            value = 2 * real(overlap.matrix[i, j])
            push!(rows, (; label_i=overlap.labels[i], label_j=overlap.labels[j], fraction=value / total * scale))
        end
    end

    sort ? Base.sort(rows; by=x -> abs(x.fraction), rev=true) : rows
end

interference_terms(xs::AbstractVector{<:OverlapContribution}; kwargs...) =
    interference_terms(mean(xs); kwargs...)
