wignerd_doublearg_sign(two_j, cosθ, ispositive) =
    wignerd_doublearg_sign.(
        two_j,
        -two_j:2:two_j,
        transpose(-two_j:2:two_j),
        cosθ,
        ispositive,
    )

PartialWaveFunctions.wignerD_doublearg(two_j, α, cosβ, γ) =
    wignerD_doublearg.(two_j, -two_j:2:two_j, transpose(-two_j:2:two_j), α, cosβ, γ)
