# Create a MizerParams object with benthic and pelagic habitats

# Created: March 3, 2021
# Last modified: September 5, 2021 by EPD

# Contents (ctrl-f):
#	0. Common functions
#	I. Two habitats with mizer


########## 0. Common functions ##########

# Add custom predation kernel
#	pred_w: predator weight in g
#	prey_w: prey weight in g
#	model: log-linear model of ppmr with pred_w
#	sigma: ppmr width
#	returns preferred ppmr
my_pred_kernel = function(pred_w, prey_w, model, sigma) {
	beta = predict(model, newdata = list(weight = pred_w), type = "response")
	ppmr = pred_w / prey_w
	phi = exp(-(log(ppmr) - beta)^2 / (2 * sigma^2))
	return(phi)
}

# Add standard lognormal predation kernel
#	beta: ppmr mean
#	sigma: ppmr width
#	pred_w: predator weight in g
#	prey_w: prey weight in g
#	returns preferred ppmr
st_pred_kernel = function(beta, sigma, pred_w, prey_w) {
	Beta = log(beta)
	ppmr = pred_w / prey_w
	phi = exp(-(log(ppmr) - Beta)^2 / (2 * sigma^2))
	return(phi)
}

########## I. Two habitats with mizer ##########

# Create a MizerParams object with benthic and pelagic habitats with user specified parameters
#	species_params: data frame output from sp_mat()
#	kappa_: resource carrying capacity coefficient
#	lambda_: resource carrying capacity scaling exponent
#	r_pp_: resource intrinsic growth rate coefficient
#	n_: resource allometric growth exponent
#	w_pp_: resource upper size cutoff
#	U_oxy_: resource tolerance of hypoxia; higher => more tolerant; scales carrying capacity
#	k_oxy_: resource oxygen at which carrying capacity is 50% of maximum
#	p: allometric metabolic exponent
#	no_w: number of size bins in the consumer spectrum
#	interaction: baseline interaction matrix, here treated as preference of each species for any resource
#	kernel: "custom" means preferred PPMR increases with body size, "standard" means standard lognormal
#	returns Mizer parameters object with two habitats, custom predation kernel, 
bp_setup = function(species_params,
			kappa_benthic,
			kappa_pelagic,
			lambda_benthic,
			lambda_pelagic,
			r_pp_benthic,
			r_pp_pelagic,
			n_benthic,
			n_pelagic,
			w_pp_cutoff_benthic,
			w_pp_cutoff_pelagic,
			U_oxy_benthic,
			U_oxy_pelagic,
			k_oxy_benthic,
			k_oxy_pelagic,
			p,
			no_w,
			interaction,
			kernel = F) {
	
	# Create params object
	params = newMultispeciesParams(species_params,
						no_w = no_w,
						kappa = kappa_benthic, 
						lambda = lambda_benthic,
						r_pp = r_pp_benthic,
						n = n_benthic,
						p = p,
						w_pp_cutoff = w_pp_cutoff_benthic,
						interaction = interaction)
	
	# All parameters for both benthic and pelagic
	params@resource_params$U_oxy = U_oxy_benthic
	params@resource_params$k_oxy = k_oxy_benthic
	params@resource_params$pelagic = list(kappa = kappa_pelagic, lambda = lambda_pelagic, r_pp = r_pp_pelagic, n = n_pelagic, w_pp_cutoff = w_pp_cutoff_pelagic, U_oxy = U_oxy_pelagic, k_oxy = k_oxy_pelagic)
	
	# Pelagic resource instrinsic growth rate
	rr_pp_pelagic = params@resource_params$pelagic$r_pp* params@w_full^(params@resource_params$pelagic$n-1)
	names(rr_pp_pelagic) = names(params@rr_pp)

	# Pelagic resource intrinsic carrying capacity
	cc_pp_pelagic = params@resource_params$pelagic$kappa*params@w_full^(-params@resource_params$pelagic$lambda)
	cc_pp_pelagic[params@w_full > params@resource_params$pelagic$w_pp_cutoff] = 0
	names(cc_pp_pelagic) = names(params@cc_pp)

	# Store in resource_params
	params@resource_params$pelagic$rr_pp = rr_pp_pelagic
	params@resource_params$pelagic$cc_pp = cc_pp_pelagic
	
	# Create custom predation kernel if kernel = T
	if(kernel) {
	
		# Create custom predation kernel
		my_pred_array = array(NA, dim = c(nrow(species_params),length(params@w),length(params@w_full)))
		my_pred_array[1,,] = t(sapply(params@w, my_pred_kernel, model = codppmr.lm, prey_w = params@w_full, sigma = params@species_params["cod","sigma"]))
		my_pred_array[2,,] = t(sapply(params@w, my_pred_kernel, model = flounderppmr.lm, prey_w = params@w_full, sigma = params@species_params["flounder","sigma"]))
		my_pred_array[3,,] = t(sapply(params@w, my_pred_kernel, model = spratppmr.lm, prey_w = params@w_full, sigma = params@species_params["sprat","sigma"]))
		my_pred_array[4,,] = t(sapply(params@w, my_pred_kernel, model = herringppmr.lm, prey_w = params@w_full, sigma = params@species_params["herring","sigma"]))
		
		params = setPredKernel(params, my_pred_array)
	}
	
	return(params)						
						
}