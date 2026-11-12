# cspell:words benchmarkable evals
using BenchmarkTools
using ThreeBodyDecays

ms = ThreeBodyMasses(1.1, 3.3, 5.5; m0 = 20.0)
msq = ms^2
σs = Invariants(ms; σ3 = 4.5^2, σ2 = 10.1^2)

suite = BenchmarkGroup()
suite["cosθ wrapper"] = @benchmarkable cosθ23($(Ref(σs))[], $(Ref(msq))[]) evals = 1
suite["cosθ indexed"] = @benchmarkable cosθij($(Ref(σs))[], $(Ref(msq))[]; k = 1) evals = 1
suite["σ wrapper"] = @benchmarkable σ3of1(0.5, $(Ref(σs.σ1))[], $(Ref(msq))[]) evals = 1
suite["σ indexed"] =
    @benchmarkable σjofk(0.5, $(Ref(σs.σ1))[], $(Ref(msq))[]; k = 1) evals = 1
suite["Kibble"] = @benchmarkable Kibble($(Ref(σs))[], $(Ref(msq))[]) evals = 1
suite["limits wrapper"] = @benchmarkable lims1($(Ref(ms))[]) evals = 1
suite["limits indexed"] = @benchmarkable lims($(Ref(ms))[]; k = 1) evals = 1
suite["Dalitz mapping"] = @benchmarkable x2σs([0.5, 0.5], $(Ref(ms))[]; k = 1) evals = 1

display(run(suite))
