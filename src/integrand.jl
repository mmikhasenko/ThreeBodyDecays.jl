"""
	three_body_phase_space_integral(function_σs, ms)
	
Computes the phase-space integral over the Dalitz plot (dsigma1 x dsigma3). Returns complex value
"""
function phase_space_integrand(function_σs, ms; k::Int)
	mssq = ms^2
	i, j = ij_from_k(k)
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


function projection_integrand(function_σs, ms, σk; k)
	l, h = lims(k, ms)
	!(l < σk < h) && return x -> 0.0
	σjlims = σjofk.([-1, 1], Ref(σk), Ref(ms^2); k)
	function integrand(x)
		σj = σjlims[1] + x[1] * (σjlims[2] - σjlims[1])
		σi = sum(ms^2) - σk - σj
		σt = circleorigin(-k, (σi, σj, σk))
		σs = MandestamTuple{typeof(ms.m0)}(σt)
		return function_σs(σs) * (σjlims[2] - σjlims[1])
	end
	return integrand
end

