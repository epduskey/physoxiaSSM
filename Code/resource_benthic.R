# Benthic resources only, scaled by oxygen
#	What is scaled?
#		1. Carrying capacity scaled as a function of benthic resource sensitivity to hypoxia
#		2. NO scaling of growth rate, given that it is an intrinsic rate
resource_benthic = function(params, n, n_pp, n_other, rates, t, dt, ...) {
    
	# Calculate scaling of carrying capacity
	cc_pp_scale = 1/(1 + exp(-1*params@resource_params$U_oxy * (params@other_params$benthic_oxygen[t+1] - params@resource_params$k_oxy)))
	
	# Benthic resource mortality with assumption of constant mortality in a given time step
	rr_pp_benthic = params@rr_pp
	cc_pp_benthic = params@cc_pp * cc_pp_scale
	mur = rr_pp_benthic + rates$resource_mort
	benthic_tmp = rr_pp_benthic * cc_pp_benthic / mur
	benthic = benthic_tmp + (n_pp - benthic_tmp) * exp(-mur * dt)
	
	# if growth rate and death rate are zero then the above would give NaN
	# whereas the value should simply not change
	sel = mur == 0
	benthic[sel] = n_pp[sel]

	return(benthic)
}