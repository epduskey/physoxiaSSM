# Pelagic resource dynamics
#	What is scaled?
#		1. Carrying capacity scaled as a function of benthic resource sensitivity to hypoxia
#		2. NO scaling of growth rate, given that it is an intrinsic rate
resource_pelagic = function(params, n, n_pp, n_other, rates, t, dt, ...) {

	# Get pelagic resource mortality
	pelagic_mort = as.vector(params@species_params$interaction_resource %*% rates$pred_rate$pelagic)
	
	# Scale according to pelagic oxygen
	cc_pp_scale = 1/(1 + exp(-params@resource_params$pelagic$U_oxy* (params@other_params$pelagic_oxygen[t+1] - params@resource_params$pelagic$k_oxy)))

	# Benthic resource mortality with assumption of constant mortality in a given time step
	rr_pp_pelagic = params@resource_params$pelagic$rr_pp
	cc_pp_pelagic = params@resource_params$pelagic$cc_pp* cc_pp_scale
	mur = rr_pp_pelagic + pelagic_mort
	pelagic_tmp = rr_pp_pelagic * cc_pp_pelagic / mur
	pelagic = pelagic_tmp - (pelagic_tmp - n_other$n_pp_pelagic) * exp(-mur * dt)
	
	# if growth rate and death rate are zero then the above would give NaN
	# whereas the value should simply not change
	sel = mur == 0
	pelagic[sel] = n_other$n_pp_pelagic[sel]
    
	return(pelagic)
}
