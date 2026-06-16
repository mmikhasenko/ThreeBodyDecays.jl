#!/usr/bin/env julia
# Benchmark heterogeneous ThreeBodyDecay chain storage (issue #121).
# Run from repo root: julia --project=benchmark benchmark/bench_tuple_chains.jl

using BenchmarkTools
using Test
using ThreeBodyDecays
using ThreeBodyDecays.Parameters: @unpack

function heterogeneous_model(n_per_topology = 2)
    tbs = ThreeBodySystem(
        ms = ThreeBodyMasses(0.141, 0.142, 0.143; m0 = 3.09),
        two_js = ThreeBodySpins(0, 0, 0; two_h0 = 2),
    )
    base = DecayChain(;
        k = 1,
        two_j = 2,
        Xlineshape = σ -> 1 / (4.1^2 - σ - 0.1im),
        Hij = RecouplingLS((2, 0)) |> VertexFunction,
        HRk = RecouplingLS((2, 2)) |> VertexFunction,
        tbs,
    )
    chain_groups = (
        ntuple(i -> DecayChain(base; k = mod1(i, 3)), n_per_topology)...,
        ntuple(
            i -> DecayChain(;
                k = mod1(i, 3),
                two_j = 2,
                Xlineshape = σ -> 1 / (5.0^2 - σ - 0.5im),
                Hij = RecouplingLS((2, 0)) |> VertexFunction,
                HRk = RecouplingLS((2, 2)) |> VertexFunction,
                tbs,
            ),
            n_per_topology,
        )...,
        ntuple(
            i -> DecayChainLS(;
                k = mod1(i, 3),
                Xlineshape = σ -> 1.0,
                jp = jp"1-",
                Ps = ThreeBodyParities('+', '+', '+'; P0 = '+'),
                tbs,
            ),
            n_per_topology,
        )...,
    )
    couplings = Tuple(ComplexF64.(randn(length(chain_groups))))
    names = Tuple("ch$i" for i in eachindex(chain_groups))
    ThreeBodyDecay(chain_groups, couplings, names)
end

function heterogeneous_model_via_vector(n_per_topology = 2)
    model = heterogeneous_model(n_per_topology)
    ThreeBodyDecay(collect(model.chains), collect(model.couplings), collect(model.names))
end

function homogeneous_model(n = 6)
    tbs = ThreeBodySystem(
        ms = ThreeBodyMasses(0.141, 0.142, 0.143; m0 = 3.09),
        two_js = ThreeBodySpins(0, 0, 0; two_h0 = 2),
    )
    base = DecayChain(;
        k = 1,
        two_j = 2,
        Xlineshape = σ -> 1 / (4.1^2 - σ - 0.1im),
        Hij = RecouplingLS((2, 0)) |> VertexFunction,
        HRk = RecouplingLS((2, 2)) |> VertexFunction,
        tbs,
    )
    chains = [DecayChain(base; k = mod1(i, 3)) for i = 1:n]
    couplings = ComplexF64.(randn(n))
    names = ["ch$i" for i = 1:n]
    ThreeBodyDecay(chains, couplings, names)
end

model_hetero = heterogeneous_model(4)  # 12 chains, 3 distinct chain types
model_hetero_vector = heterogeneous_model_via_vector(4)
model_homo = homogeneous_model(12)

ms = masses(model_hetero)
σs = x2σs([0.1, 0.9], ms; k = 2)

println("=== Type stability ===")
println("heterogeneous chains field: ", fieldtype(typeof(model_hetero), :chains))
println("heterogeneous model type:   ", typeof(model_hetero))
println(
    "heterogeneous (vector path) chains field: ",
    fieldtype(typeof(model_hetero_vector), :chains),
)
println("homogeneous chains field:   ", fieldtype(typeof(model_homo), :chains))

@inferred ComplexF64 amplitude(model_hetero, σs, spins(model_hetero))
@inferred Float64 unpolarized_intensity(model_hetero, σs)

println("\n=== Benchmark: heterogeneous 12-chain model (tuple storage) ===")
t_hetero = @belapsed unpolarized_intensity($model_hetero, $σs)
println("unpolarized_intensity: ", round(t_hetero * 1e6; digits = 2), " µs")

println("\n=== Benchmark: heterogeneous 12-chain model (vector constructor) ===")
t_hetero_vec = @belapsed unpolarized_intensity($model_hetero_vector, $σs)
println("unpolarized_intensity: ", round(t_hetero_vec * 1e6; digits = 2), " µs")

println("\n=== Benchmark: homogeneous 12-chain model ===")
t_homo = @belapsed unpolarized_intensity($model_homo, $σs)
println("unpolarized_intensity: ", round(t_homo * 1e6; digits = 2), " µs")

println("\n=== Benchmark: 20-chain vcat model (test_ls_amplitude setup) ===")
two_js, Ps = ThreeBodySpinParities("1-", "1/2+", "0-"; jp0 = "3/2+")
tbs20 = ThreeBodySystem(1.1, 2.2, 3.3; m0 = 7.7, two_js)
models20 = map([
    (name = "R3_3h-", k = 3, two_jp = jp"3/2-"),
    (name = "R1_3h-", k = 1, two_jp = jp"3/2-"),
    (name = "R3_1h-", k = 3, two_jp = jp"1/2-"),
    (name = "R1_1h-", k = 1, two_jp = jp"1/2-"),
    (name = "R2", k = 2, two_jp = jp"3-"),
]) do (; k, two_jp, name)
    chains = map(possible_lsLS(two_jp, two_js, Ps; k)) do conf
        @unpack two_LS, two_ls = conf
        DecayChain(;
            k,
            two_jp.two_j,
            tbs = tbs20,
            Xlineshape = identity,
            HRk = RecouplingLS(two_LS) |> VertexFunction,
            Hij = RecouplingLS(two_ls) |> VertexFunction,
        )
    end
    ThreeBodyDecay(name .=> zip(ones(length(chains)), chains))
end
model20 = vcat(models20...)
σs20 = x2σs([0.5, 0.3], masses(model20); k = 1)
t20 = @belapsed unpolarized_intensity($model20, $σs20)
println("unpolarized_intensity (20 chains): ", round(t20 * 1e3; digits = 2), " ms")
println("model20 chains field: ", fieldtype(typeof(model20), :chains))
