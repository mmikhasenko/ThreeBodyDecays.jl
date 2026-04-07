using Test
using Random
using Statistics
using LinearAlgebra
using ThreeBodyDecays

function overlap_test_model()
    tbs = ThreeBodySystem(
        ms = ThreeBodyMasses(0.141, 0.142, 0.143; m0 = 3.09),
        two_js = ThreeBodySpins(0, 0, 0; two_h0 = 2),
    )

    ch1 = DecayChain(;
        k = 1,
        two_j = 2,
        Xlineshape = σ -> 1 / (4.1^2 - σ - 0.1im),
        Hij = RecouplingLS((2, 0)) |> VertexFunction,
        HRk = RecouplingLS((2, 2)) |> VertexFunction,
        tbs,
    )
    ch2 = DecayChain(ch1; k = 2)
    ch3 = DecayChain(ch1; k = 3)

    chains = [ch1, ch2, ch3]
    couplings = ComplexF64[4.0, 2.0, 3.0]
    names = ["K(892)", "K(892)", "L(1520)"]
    ThreeBodyDecay(chains, couplings, names)
end

function physical_samples(model, nsamples)
    ms = masses(model)
    samples = typeof(y2σs(rand(2), ms))[]
    while length(samples) < nsamples
        σs = y2σs(rand(2), ms)
        isphysical(σs, ms) && push!(samples, σs)
    end
    samples
end

function grouped_submodels(model)
    labels = unique(model.names)
    labels, [model[model.names .== label] for label in labels]
end

function collapse_chain_matrix(result, model)
    labels = unique(model.names)
    lookup = Dict(label => i for (i, label) in enumerate(labels))
    matrix = zeros(ComplexF64, length(labels), length(labels))
    for i in eachindex(model.chains), j in eachindex(model.chains)
        ii = lookup[model.names[i]]
        jj = lookup[model.names[j]]
        matrix[ii, jj] += result.matrix[i, j]
    end
    labels, matrix
end

Random.seed!(1234)
model = overlap_test_model()
samples = physical_samples(model, 24)
cache = chain_amplitudes(model, samples)
event_chain_overlaps = event_overlap_contributions(cache)

chain_stripped = chain_overlap_contribution(cache)
chain_physical = physical_overlap(chain_stripped, model.couplings)
name_physical = group_overlap(chain_physical, model.names)

@testset "Hermitian structure" begin
    @test chain_physical.matrix ≈ chain_physical.matrix'
    @test chain_stripped.matrix ≈ chain_stripped.matrix'
    @test name_physical.matrix ≈ name_physical.matrix'
    @test maximum(
        abs(imag(chain_physical.matrix[i, i])) for i in eachindex(chain_physical.labels)
    ) ≤ 1e-10
    @test all(
        real(chain_physical.matrix[i, i]) ≥ -1e-10 for i in eachindex(chain_physical.labels)
    )
end

@testset "Event-level overlaps reduce correctly" begin
    mean_chain = mean(event_chain_overlaps)
    @test Matrix(mean_chain.matrix) ≈
          Matrix(chain_stripped.matrix) ./ length(event_chain_overlaps) atol = 1e-10 rtol =
        1e-10

    se_chain = ThreeBodyDecays.stderr(event_chain_overlaps)
    @test se_chain.labels == chain_stripped.labels
    @test size(se_chain.re) == size(chain_stripped.matrix)
    @test size(se_chain.im) == size(chain_stripped.matrix)
    @test all(se_chain.re .≥ 0.0)
    @test all(se_chain.im .≥ 0.0)
end

@testset "Quadratic form reproduces direct intensities" begin
    direct_total = sum(σs -> unpolarized_intensity(model, σs), samples)
    @test total_intensity(chain_physical) ≈ direct_total atol = 1e-9 rtol = 1e-10
    @test real(model.couplings' * Matrix(chain_stripped.matrix) * model.couplings) ≈
          direct_total atol = 1e-9 rtol = 1e-10

    random_couplings =
        complex.(randn(length(model.couplings)), randn(length(model.couplings)))
    random_model =
        ThreeBodyDecay(collect(model.chains), random_couplings, collect(model.names))
    random_total = sum(σs -> unpolarized_intensity(random_model, σs), samples)
    @test real(random_couplings' * Matrix(chain_stripped.matrix) * random_couplings) ≈
          random_total atol = 1e-8 rtol = 1e-10
end

@testset "Grouped overlaps match collapsed chain overlaps" begin
    labels, collapsed = collapse_chain_matrix(chain_physical, model)
    @test labels == name_physical.labels
    @test collapsed ≈ Matrix(name_physical.matrix) atol = 1e-9 rtol = 1e-10
end

@testset "Fit fractions and pair interferences match direct grouped sums" begin
    labels, submodels = grouped_submodels(model)
    direct_diagonal = Dict(
        label => sum(σs -> unpolarized_intensity(submodel, σs), samples) for
        (label, submodel) in zip(labels, submodels)
    )
    direct_pairs = Dict{Tuple{String,String},Float64}()
    for j = 2:length(labels), i = 1:(j-1)
        pair_model = model[(model.names .== labels[i]) .| (model.names .== labels[j])]
        pair_total = sum(σs -> unpolarized_intensity(pair_model, σs), samples)
        direct_pairs[(labels[i], labels[j])] =
            pair_total - direct_diagonal[labels[i]] - direct_diagonal[labels[j]]
    end

    fractions = Dict(
        row.label => row.fraction for
        row in fit_fractions(name_physical; percent = false, sort = false)
    )
    for label in labels
        @test fractions[label] ≈ direct_diagonal[label] / total_intensity(name_physical) atol =
            1e-10 rtol = 1e-10
    end

    pair_rows = Dict(
        (row.label_i, row.label_j) => row.fraction for
        row in interference_terms(name_physical; percent = false, sort = false)
    )
    for (pair, value) in direct_pairs
        @test pair_rows[pair] ≈ value / total_intensity(name_physical) atol = 1e-9 rtol =
            1e-10
    end
end

@testset "Vector APIs aggregate before normalization" begin
    event_name_physical = group_overlap.(
        physical_overlap.(event_chain_overlaps, Ref(model.couplings)),
        Ref(model.names),
    )

    fractions_from_vector = Dict(
        row.label => row.fraction for
        row in fit_fractions(event_name_physical; percent = false, sort = false)
    )
    fractions_from_mean = Dict(
        row.label => row.fraction for
        row in fit_fractions(mean(event_name_physical); percent = false, sort = false)
    )
    @test fractions_from_vector == fractions_from_mean

    interferences_from_vector = Dict(
        (row.label_i, row.label_j) => row.fraction for
        row in interference_terms(event_name_physical; percent = false, sort = false)
    )
    interferences_from_mean = Dict(
        (row.label_i, row.label_j) => row.fraction for row in
        interference_terms(mean(event_name_physical); percent = false, sort = false)
    )
    @test interferences_from_vector == interferences_from_mean
end

@testset "Overlap utilities preserve real element type" begin
    labels = ["a", "b"]
    m1 = Hermitian(ComplexF32[1+0im 0.25+0.5im; 0.25-0.5im 2+0im], :U)
    m2 = Hermitian(ComplexF32[2+0im 0.5+0.25im; 0.5-0.25im 3+0im], :U)
    xs = OverlapContribution{Float32}.([labels, labels], [m1, m2])

    μ = mean(xs)
    σ = std(xs; corrected = false)
    grouped = group_overlap(μ, ["all", "all"])

    @test μ isa OverlapContribution{Float32}
    @test eltype(Matrix(μ.matrix)) === ComplexF32
    @test eltype(σ.re) === Float32
    @test eltype(σ.im) === Float32
    @test grouped isa OverlapContribution{Float32}
    @test eltype(Matrix(grouped.matrix)) === ComplexF32
    @test total_intensity(grouped) isa Float32
    @test fit_fractions(grouped; percent = false, sort = false)[1].fraction isa Float32
    @test interference_terms(grouped; percent = false, sort = false) == NamedTuple[]
end
