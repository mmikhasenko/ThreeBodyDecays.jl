
#                                                                _|
#    _|_|_|    _|_|    _|_|_|      _|_|    _|  _|_|    _|_|_|  _|_|_|_|    _|_|
#  _|    _|  _|_|_|_|  _|    _|  _|_|_|_|  _|_|      _|    _|    _|      _|_|_|_|
#  _|    _|  _|        _|    _|  _|        _|        _|    _|    _|      _|
#    _|_|_|    _|_|_|  _|    _|    _|_|_|  _|          _|_|_|      _|_|    _|_|_|
#        _|
#    _|_|


function flatDalitzPlotSample31(tbs; Nev::Int=10000)
    density = getbinned1dDensity(σ1->sqrt(λ(σ1,tbs.msq[2],tbs.msq[3])*λ(σ1,tbs.s,tbs.msq[1]))/σ1, (tbs.mthsq[1],tbs.sthsq[1]), 500)
    σ1 = [rand(density) for _ in 1:Nev]
    σ3 = [σ3of1(σ,2*rand()-1,tbs) for σ in σ1]
    return σ3, σ1
end
