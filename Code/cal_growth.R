# Fitting growth

# A function to fit growth with h or gamma
#	par: input parameters c(all species' maximum intake, all species' search volume) in that order
#	params: mizer params object
#	horg: fit parameter "h" or "gamma"
#	species: which species to fit
#	yr: years to which to fit
#	t: time steps of model to run
#	effort: a vector containing effort values for cod, flounder, sprat, and herring in that order
#	sim: mizer simulation object, if already run
#	cal: if true, takes mean values over all years
#	returns sum of relative errors in all species' SSB 
fit_growth = function(par, params, horg, species, yr = NULL, t, effort, sim = NULL, cal = T) {
	
	# New params object
	ret = params
	
	# Run simulation if not already provided
	if(is.null(sim)) {
		# Assign parameters
		if(horg == "gamma") {
			gamma = vector(mode = "numeric", length = 4)
			for(i in 1:length(gamma)) {
				gamma[i] = ifelse(rownames(species_params(params))[i] %in% species, 10^par[i], species_params(params)$gamma[i])
			}
			species_params(ret)$gamma = gamma
		} else if(horg == "h") {
			h = vector(mode = "numeric", length = 4)
			for(i in 1:length(gamma)) {
				h[i] = ifelse(rownames(species_params(params))[i] %in% species, par[i], species_params(params)$h[i])
			}
			species_params(ret)$h = h
			sp = species_params(ret)
			species_params(ret)$ks = 0.2 * sp$alpha * sp$h * sp$w_mat^(sp$n - sp$p)
		} else{
			stop("horg must be equal to 'gamma' or 'h'")
		}
	
		# Run to steady state
		sim = project(ret, t_max = t, effort = effort)
	}
	
	# Get model growth
	if(cal) {
		
		# Set up params object
		params_growth = sim@params
		params_growth@initial_n[] <- sim@n[dim(sim@n)[1],, ]
    		params_growth@initial_n_pp[] <- sim@n_pp[dim(sim@n)[1], ]
    		params_growth@initial_n_other[] <- sim@n_other[dim(sim@n)[1], ]
    		params_growth@initial_effort[] <- sim@effort[dim(sim@n)[1], ]
    		params_growth@time_modified <- lubridate::now()
	
		# Model growth
		sp = species_params(params_growth)[rownames(species_params(params_growth)) %in% species, ]
		sp_names = rownames(sp)
		sp_maxage = c(15,26,16,13)[rownames(species_params(params_growth)) %in% species]
		growth_mod = array(NA, dim = c(length(sp_maxage),50))
		dimnames(growth_mod) = list(Species = sp_names, Age = seq(0,1,length.out=50))
		
		# von Bertalanffy growth
		L_inf = (sp$w_inf/sp$a)^(1/sp$b)
		growth_obs = matrix(NA, nrow=length(sp_names), ncol = 50)
		for(i in 1:nrow(growth_obs)) {
			growth_mod[i,] = getGrowthCurves(params_growth, species = sp_names[i], max_age = sp_maxage[i], percentage = T)
			growth_obs[i,] = sp$a[i] * (L_inf[i] * (1 - exp(-sp$k_vb[i] * (as.numeric(dimnames(growth_mod)$Age)*sp_maxage[i])))) ^ sp$b[i]
			growth_obs[i,] = (growth_obs[i,]/sp$w_inf[i])*100
		}
		
		error = sum((growth_mod - growth_obs)^2)
		return(error)
	} else {
		
		# Growth matrices
		growth_mod = array(dim = c(length(yr), length(species), 50))
		species_all = c("cod","flounder","sprat","herring")
		dimnames(growth_mod) = list(Year = yr, Species = species_all[species_all %in% species], Age = seq(0,1,length.out=50))
		sp_maxage = c(15,26,16,13)
		growth_obs = growth_mod
		
		for(i in 1:dim(growth_mod)[1]) {
			
			# Set up params object
			params_growth = sim@params
			params_growth@initial_n[] <- sim@n[dim(sim@n)[1]-length(yr)+i,, ]
    			params_growth@initial_n_pp[] <- sim@n_pp[dim(sim@n)[1]-length(yr)+i, ]
    			params_growth@initial_n_other[] <- sim@n_other[dim(sim@n)[1]-length(yr)+i, ]
    			params_growth@initial_effort[] <- sim@effort[dim(sim@n)[1]-length(yr)+i, ]
    			params_growth@time_modified <- lubridate::now()
			
			# Model growth
		
			# von Bertalanffy growth
			sp = species_params(params_growth)[rownames(species_params(params_growth)) %in% species, ]
			L_inf = (sp$w_inf/sp$a)^(1/sp$b)
			for(j in 1:dim(growth_obs)[2]) {
				growth_mod[i,j,] = getGrowthCurves(params_growth, species = rownames(sp)[j], max_age = sp_maxage[j], percentage = T)
				growth_obs[i,j,] = sp$a[j] * (L_inf[j] * (1 - exp(-sp$k_vb[j] * (as.numeric(dimnames(growth_mod)$Age)*sp_maxage[j])))) ^ sp$b[j]
				growth_obs[i,j,] = (growth_obs[i,j,]/sp$w_inf[j])*100
			}
		}
		
		error = sum((growth_mod - growth_obs)^2)
		return(error)
	}
}

# A function to fit growth with h or gamma
#	par: input parameters c(U_crit, a_crit, U_met, a_met) in that order
#	params: mizer params object
#	species: which species to fit
#	yr: years to which to fit
#	t: time steps of model to run
#	effort: a vector containing effort values for cod, flounder, sprat, and herring in that order
#	sim: mizer simulation object, if already run
#	returns sum of relative errors in all species' SSB 
time_growth = function(par, params, species, yr, t, effort, sim = NULL) {
	
	# New params object
	ret = params
	
	# Other variables
	us_cod = species_params(ret)["cod",]$U_search
	us_flounder = species_params(ret)["flounder",]$U_search
	as_cod = species_params(ret)["cod",]$a_search
	as_flounder = species_params(ret)["flounder",]$a_search
	ui_cod = species_params(ret)["cod",]$U_maxin
	ui_flounder = species_params(ret)["flounder",]$U_maxin
	ai_cod = species_params(ret)["cod",]$a_maxin
	ai_flounder = species_params(ret)["flounder",]$a_maxin
	ue_cod = species_params(ret)["cod",]$U_erepro
	ue_flounder = species_params(ret)["flounder",]$U_erepro
	ae_cod = species_params(ret)["cod",]$a_erepro
	ae_flounder = species_params(ret)["flounder",]$a_erepro
	ua_cod = species_params(ret)["cod",]$U_alpha
	ua_flounder = species_params(ret)["flounder",]$U_alpha
	aa_cod = species_params(ret)["cod",]$a_alpha
	aa_flounder = species_params(ret)["flounder",]$a_alpha
	um_cod = species_params(ret)["cod",]$U_met
	um_flounder = species_params(ret)["flounder",]$U_met
	am_cod = species_params(ret)["cod",]$a_met
	am_flounder = species_params(ret)["flounder",]$a_met
	zm_cod = species_params(ret)["cod",]$z_mort
	zm_flounder = species_params(ret)["flounder",]$z_mort
	bm_cod = species_params(ret)["cod",]$b_mort
	bm_flounder = species_params(ret)["flounder",]$z_mort
	
	# Set up occupancy optim object if absent
	if(!("occ.optim" %in% ls())) { 
		occ.optim = list()
		occ.optim$par = c(10^5, 0)
	}
	
	# Run simulation if not already provided
	if(is.null(sim)) {
		# Assign parameters
		if(species == "cod") {
			ret = oxy_sensitivity(ret, 
				U_hab = c(occ.optim$par[1], 10^5, -10^5, -10^5),
				a_hab = c(occ.optim$par[2], 0, 0, 0),
				U_search = c(par[1], us_flounder, 10^5, 10^5),
				a_search = c(par[2], as_flounder, 0, 0),
				U_maxin = c(par[3], ui_flounder, 10^5, 10^5),
				a_maxin = c(par[4], ai_flounder, 0, 0),
				U_erepro = c(par[5], ue_flounder, 10^5, 10^5),
				a_erepro = c(par[6], ae_flounder, 0, 0),
				U_alpha = c(par[7], ua_flounder, 10^5, 10^5),
				a_alpha = c(par[8], aa_flounder, 0, 0),
				U_met = c(par[9], um_flounder, 10^5, 10^5),
				a_met = c(par[10], am_flounder, -10^5, -10^5),
				z_mort = c(zm_cod, zm_flounder, 0, 0),
				b_mort = c(bm_cod, bm_flounder, 0, 0))
		} else if(species == "flounder") {
			ret = oxy_sensitivity(ret, 
				U_hab = c(occ.optim$par[1], 10^5, -10^5, -10^5),
				a_hab = c(occ.optim$par[2], 0, 0, 0),
				U_search = c(us_cod, par[1], 10^5, 10^5),
				a_search = c(as_cod, par[2], 0, 0),
				U_maxin = c(ui_cod, par[3], 10^5, 10^5),
				a_maxin = c(ai_cod, par[4], 0, 0),
				U_erepro = c(ue_cod, par[5], 10^5, 10^5),
				a_erepro = c(ae_cod, par[6], 0, 0),
				U_alpha = c(ua_cod, par[7], 10^5, 10^5),
				a_alpha = c(aa_cod, par[8], 0, 0),
				U_met = c(um_cod, par[9], 10^5, 10^5),
				a_met = c(am_cod, par[10], -10^5, -10^5),
				z_mort = c(zm_cod, zm_flounder, 0, 0),
				b_mort = c(bm_cod, bm_flounder, 0, 0))
		}
		
		# Create an oxygen vector
		times = 0:t
		benthic_oxygen = vector(mode = "numeric", length = length(times))
		pelagic_oxygen = vector(mode = "numeric", length = length(times))

		# Time series oxygen
		benthic_oxygen[1:length(times)] = c(rep(mean(preds_all$df[preds_all$df$Year %in% yr, ]$Oxygen_benthic),t-length(yr)+1), preds_all$df[preds_all$df$Year %in% yr, ]$Oxygen_benthic)
		pelagic_oxygen[1:length(times)] = c(rep(mean(preds_all$df[preds_all$df$Year %in% yr, ]$Oxygen_pelagic),t-length(yr)+1), preds_all$df[preds_all$df$Year %in% yr, ]$Oxygen_pelagic)

		# Store these in other_params
		ret@other_params$benthic_oxygen = benthic_oxygen
		ret@other_params$pelagic_oxygen = pelagic_oxygen

		# Choose reference temperature
		Tref = 19
		ret@other_params$Tref = Tref

		# Create a temperature vector
		benthic_temp = vector(mode = "numeric", length = length(times))
		pelagic_temp = vector(mode = "numeric", length = length(times))

		# Constant temperature scenario
		benthic_temp[1:length(times)] = 19
		pelagic_temp[1:length(times)] = 19

		# Store these in other_params
		ret@other_params$benthic_temp = benthic_temp
		ret@other_params$pelagic_temp = pelagic_temp

		# Scale rates by oxygen and temperature: see otmscale.R
		ret = rate_scale(ret, t)
		ret = mscale(ret, t)
	
		# Run to steady state
		sim = project(ret, t_max = t, effort = effort)
	}
		
	# Growth matrices
	growth_mod = array(dim = c(length(yr), length(species), 50))
	species_all = c("cod","flounder","sprat","herring")
	dimnames(growth_mod) = list(Year = yr, Species = species_all[species_all %in% species], Age = seq(0,1,length.out=50))
	sp_maxage = c(15,26,16,13)
	growth_obs = growth_mod
	
	# Growth data
	sp_obs = list(cod = cod_ivb_cal, flounder = flounder_ivb_cal, sprat = sprat_ivb_cal, herring = herring_ivb_cal)[species]
	
	for(i in 1:dim(growth_mod)[1]) {
		
		# Model growth
	
		# von Bertalanffy growth
		sp = species_params(ret)[rownames(species_params(ret)) %in% species, ]
		L_inf = (sp$w_inf/sp$a)^(1/sp$b)
		for(j in 1:dim(growth_obs)[2]) {
			growth_mod[i,j,] = myGrowthCurves(sim, (t-length(yr))+i, species = rownames(sp)[j], max_age = sp_maxage[j], percentage = T)
			growth_obs[i,j,] = sp_obs[[j]][i,]$a * (sp_obs[[j]][i,]$L_inf * (1 - exp(-sp_obs[[j]][i,]$k * (as.numeric(dimnames(growth_mod)$Age)*sp_maxage[j])))) ^ sp_obs[[j]][i,]$b
			growth_obs[i,j,] = (growth_obs[i,j,]/sp$w_inf[j])*100
		}
	}
		
	error = sum((growth_mod - growth_obs)^2)
	return(error)
}

# A function to return growth error
#	params: mizer params object
#	species: which species to fit
#	yr: years to which to fit
#	t: time steps of model to run
#	effort: a vector containing effort values for cod, flounder, sprat, and herring in that order
#	sim: mizer simulation object, if already run
#	returns sum of relative errors in all species' SSB 
proj_growth = function(params, yr, t, species, sim = NULL) {
	
	# Run simulation if not already provided
	if(is.null(sim)) {
		# Create effort array
		effort = farray(params, fdat, yr, c(2,2,2,2), t)

		# Run to steady state
		sim = project(params, t_max = t, effort = effort)
	}
		
	# Growth matrices
	growth_mod = array(dim = c(length(yr), length(species), 50))
	dimnames(growth_mod) = list(Year = yr, Species = species, Age = seq(0,1,length.out=50))
	sp_maxage = c(15,26,16,13)
	growth_obs = growth_mod
	
	# Growth data
	sp_obs = list(cod = rbind(cod_ivb_cal,cod_ivb_prj), flounder = rbind(flounder_ivb_cal,flounder_ivb_prj), sprat = rbind(sprat_ivb_cal,sprat_ivb_prj), herring = rbind(herring_ivb_cal,herring_ivb_prj))[species]
	
	for(i in 1:dim(growth_mod)[1]) {
		
		# Model growth
	
		# von Bertalanffy growth
		sp = species_params(params)[rownames(species_params(params)) %in% species, ]
		L_inf = (sp$w_inf/sp$a)^(1/sp$b)
		for(j in 1:dim(growth_obs)[2]) {
			sp_temp = sp_obs[[j]][sp_obs[[j]]$Year %in% yr[i], 2:6]
			growth_mod[i,j,] = myGrowthCurves(sim, (t-length(yr))+i, species = rownames(sp)[j], max_age = sp_maxage[j], percentage = T)
			if(sum(is.na(sp_temp)) == 0) {
				growth_obs[i,j,] = sp_temp$a * (sp_temp$L_inf * (1 - exp(-sp_temp$k * (as.numeric(dimnames(growth_mod)$Age)*sp_maxage[j])))) ^ sp_temp$b
			} else {
				growth_obs[i,j,] = NA
				growth_mod[i,j,] = NA
			}
			growth_obs[i,j,] = (growth_obs[i,j,]/sp$w_inf[j])*100
		}
	}
		
	error = sum((growth_mod - growth_obs)^2, na.rm = T)
	return(error)
}

# A function to fit growth with h or gamma
#	par: benthic oxygen sensitivity input parameters c(U_oxy, k_oxy) in that order
#	params: mizer params object
#	yr: years to which to fit
#	t: time steps of model to run
#	effort: a vector containing effort values for cod, flounder, sprat, and herring in that order
#	sim: mizer simulation object, if already run
#	returns sum of relative errors in all species' growth over time 
benthic_growth = function(par, params, yr, t, effort, sim = NULL) {
	
	# New params object
	ret = params
	
	# Run simulation if not already provided
	if(is.null(sim)) {	
		# Assign parameters
		resource_params(ret)$U_oxy = par[1]
		resource_params(ret)$k_oxy = par[2]
		
		# Create an oxygen vector
		times = 0:t
		benthic_oxygen = vector(mode = "numeric", length = length(times))
		pelagic_oxygen = vector(mode = "numeric", length = length(times))

		# Time series oxygen
		benthic_oxygen[1:length(times)] = c(rep(mean(preds_all$df[preds_all$df$Year %in% yr, ]$Oxygen_benthic),t-length(yr)+1), preds_all$df[preds_all$df$Year %in% yr, ]$Oxygen_benthic)
		pelagic_oxygen[1:length(times)] = c(rep(mean(preds_all$df[preds_all$df$Year %in% yr, ]$Oxygen_pelagic),t-length(yr)+1), preds_all$df[preds_all$df$Year %in% yr, ]$Oxygen_pelagic)

		# Store these in other_params
		ret@other_params$benthic_oxygen = benthic_oxygen
		ret@other_params$pelagic_oxygen = pelagic_oxygen

		# Choose reference temperature
		Tref = 19
		ret@other_params$Tref = Tref

		# Create a temperature vector
		benthic_temp = vector(mode = "numeric", length = length(times))
		pelagic_temp = vector(mode = "numeric", length = length(times))

		# Constant temperature scenario
		benthic_temp[1:length(times)] = 19
		pelagic_temp[1:length(times)] = 19

		# Store these in other_params
		ret@other_params$benthic_temp = benthic_temp
		ret@other_params$pelagic_temp = pelagic_temp

		# Scale rates by oxygen and temperature: see otmscale.R
		ret = rate_scale(ret, t)
		ret = mscale(ret, t)
	
		# Run to steady state
		sim = project(ret, t_max = t, effort = effort)
	}
		
	# Growth matrices
	growth_mod = array(dim = c(length(yr), nrow(params@species_params), 50))
	species_all = c("cod","flounder","sprat","herring")
	dimnames(growth_mod) = list(Year = yr, Species = species_all, Age = seq(0,1,length.out=50))
	sp_maxage = c(15,26,16,13)
	growth_obs = growth_mod
	
	# Growth data
	sp_obs = list(cod = cod_ivb_cal, flounder = flounder_ivb_cal, sprat = sprat_ivb_cal, herring = herring_ivb_cal)
	
	for(i in 1:dim(growth_mod)[1]) {
		
		# Model growth
	
		# von Bertalanffy growth
		sp = species_params(ret)
		L_inf = (sp$w_inf/sp$a)^(1/sp$b)
		for(j in 1:dim(growth_obs)[2]) {
			growth_mod[i,j,] = myGrowthCurves(sim, (t-length(yr))+i, species = rownames(sp)[j], max_age = sp_maxage[j], percentage = T)
			growth_obs[i,j,] = sp_obs[[j]][i,]$a * (sp_obs[[j]][i,]$L_inf * (1 - exp(-sp_obs[[j]][i,]$k * (as.numeric(dimnames(growth_mod)$Age)*sp_maxage[j])))) ^ sp_obs[[j]][i,]$b
			growth_obs[i,j,] = (growth_obs[i,j,]/sp$w_inf[j])*100
		}
	}
		
	error = sum((growth_mod - growth_obs)^2, na.rm = T)
	return(error)
}