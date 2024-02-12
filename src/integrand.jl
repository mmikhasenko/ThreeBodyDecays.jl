"""
    three_body_phase_space_integral(function_σs, ms)
    
Computes the phase-space integral over the Dalitz plot (dsigma1 x dsigma3). Returns complex value
"""
function phase_space_integrand_ki(function_σs, ms; k::Int)
    mssq = ms^2
    i,j = ij_from_k(k)
    misq, mjsq, mksq, m0sq = mssq[i], mssq[j], mssq[k], mssq[4]
    σkmin, σkmax = lims(k, ms)
    # 
    function integrand(x)
        σs = x2σs(x, ms; k)
        σk = σs[k]
        jac_z = sqrt(Kallen(m0sq, σk, mksq) * Kallen(σk, misq, mjsq)) / σk
        jac_σk = (σkmax - σkmin)
        value = function_σs(σs) * jac_σk * jac_z / m0sq
        return value
    end
    return integrand
end




