# Scale maintenance costs with oxygen
#	What is scaled?
#		1. Maintenance costs (metabolic rate) have a hump-shaped relationship with falling oxygen
#		2. Scale assimilation efficiency as a representation of increasing costs of digestion, ventilation, movement, and waste management
physoxiaEReproAndGrowth <- function(params, n, n_pp, n_other, t, encounter, feeding_level, ...) {
	alpha = params@species_params$alpha * params@other_params$rs_alpha[t+1,,]
	erag = sweep((1 - feeding_level) * encounter, 1, alpha, "*", check.margin = FALSE) - (params@metab * params@other_params$metab_scale[t+1,,])
	return(erag)
}