# Mortality function with added direct mortality due to a hazard function of hypoxia exposure
#	What is scaled?
#		1. Nothing, this is additional mortality as a function of oxygen
hypoxiaMort = function(params, n, n_pp, n_other, t, f_mort, pred_mort, ...) {
	
	# Mortality at zero oxygen
	z0 = params@species_params$z_mort
	
	# Hypoxia tolerance; higher is greater tolerance
	b = params@species_params$b_mort
	
	# Calculate direct mortality due to hypoxia exposure
	hypoxia_mort = z0*exp(-1*b*params@other_params$oxy_exposure[t+1,,])
    
	# All mortality
	mort = pred_mort + params@mu_b + f_mort + hypoxia_mort
    
	# Add contributions from other components
	for (i in seq_along(params@other_mort)) {
		mort = mort + do.call(params@other_mort[[i]], 
        						list(params = params,
        								n = n, 
        								n_pp = n_pp, 
        								n_other = n_other, 
        								t = t,
        								component = names(params@other_mort)[[i]], ...))
	}
	return(mort)
}