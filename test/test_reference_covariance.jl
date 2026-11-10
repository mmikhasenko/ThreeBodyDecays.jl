using ThreeBodyDecays, Random, LinearAlgebra, Test
# This checks covariance of the package's fixed chain conventions, not a conversion
# recipe for arbitrary fitted couplings. Each model owns its reference at creation.
@testset "Reference tensor covariance" begin
    Random.seed!(12876)
    for ms in (
            ThreeBodyMasses(1, 2, 3; m0 = 7.9),
            ThreeBodyMasses(0.938, 0.49367, 0.13957; m0 = 2.46867),
        ),
        sample = 1:8

        s=randomPoint(ms)
        function d(a, b, p, j)
            w=wr(a, b, p)
            ThreeBodyDecays.wignerd_doublearg_sign(
                j,
                cosζ(w, s, ms^2),
                ThreeBodyDecays.ispositive(w),
            )
        end
        for j = 1:4, a = 1:3, b = 1:3, c = 1:3, p = 0:3
            # Half-integer parent, three distinct topologies: the composition
            # crosses the 2↔3 rotation, which lacks the SU(2) full-turn sign (#127).
            if p == 0 && isodd(j) && length(Set((a, b, c))) == 3
                @test_broken d(a, c, p, j) ≈ d(a, b, p, j)*d(b, c, p, j) atol=1e-10
            else
                @test d(a, c, p, j) ≈ d(a, b, p, j)*d(b, c, p, j) atol=1e-10
            end
        end
        tbs=ThreeBodySystem(ms, ThreeBodySpins(1, 2, 0; two_h0 = 1))
        c1=DecayChain(
            k = 1,
            Xlineshape = x->1+0im,
            two_j = 4,
            HRk = Vertex(RecouplingLS((4, 3))),
            Hij = Vertex(RecouplingLS((4, 2))),
            tbs = tbs,
        )
        c3=DecayChain(
            k = 3,
            Xlineshape = x->1+0im,
            two_j = 1,
            HRk = Vertex(RecouplingLS((2, 1))),
            Hij = Vertex(RecouplingLS((2, 3))),
            tbs = tbs,
        )
        models=[
            ThreeBodyDecay(
                ["a", "b"] .=> zip([1+0.37im, -0.62+0.81im], [c1, c3]);
                reference_k = k,
            ) for k = 1:3
        ]
        A=[amplitude(m, s) for m in models]
        for r = 2:3
            D=[d(1, r, p, j) for (p, j) in zip((1, 2, 3, 0), (1, 2, 0, 1))]
            transformed=map(CartesianIndices(A[1])) do out
                sum(
                    D[4][out[4], inp[4]]*A[1][inp]*prod(D[p][inp[p], out[p]] for p = 1:3)
                    for inp in CartesianIndices(A[1])
                )
            end
            # Reference 2: the topology-3 chain needs the 3→2 parent rotation,
            # which differs by the full-turn sign from 3→1→2 (#127).
            if r == 2
                @test_broken transformed ≈ A[r] atol=1e-11
            else
                @test transformed ≈ A[r] atol=1e-11
            end
        end
        @test_broken sum(abs2, A[1]) ≈ sum(abs2, A[2]) atol=1e-11
        @test sum(abs2, A[1]) ≈ sum(abs2, A[3]) atol=1e-11
    end
end
