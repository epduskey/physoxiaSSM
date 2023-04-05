# Fitting yield

# A function to fit yield
#	par: input parameters c(background mortality for each species)
#	params: mizer params object
#	dat: a list containing data frames with "Year" and "SSB" columns for cod, flounder, sprat, and herring, in that order
#	yr: years to which to fit
#	t: number of time steps for the model
#	effort: a vector containing effort values for cod, flounder, sprat, and herring in that order
#	sim: mizer simulation object, if already run
#	cal: if true, takes mean values over all years
#	returns sum of relative errors in all species' SSB 
fit_yield = function(par, params, dat, yr, t, effort, sim = NULL, cal = T) {
		
	# New params object
	ret = params

	# Run simulation if not already provided
	if(is.null(sim)) {
		
		# Reset parameter values
		species_params(ret)$z0 = par
	
		# Run to steady state
		sim = project(ret, t_max = t, effort = effort)
	}
	
	# Set up yield matrices
	yield_obs = matrix(nrow = ifelse(cal,1,length(yr)), ncol = 4)
	rownames(yield_obs) = if(cal) {"cal"} else{yr}
	colnames(yield_obs) = c("cod", "flounder", "sprat", "herring")
	idx_mod = if(cal) {t} else{(t-length(yr)+1):t}
	yield_mod = getYield(sim)[idx_mod, ]/1000/1000
	
	# Get yield for each species
	for(i in 1:ncol(yield_obs)) {
		
		# Use weights greater than w_mat
		idx_obs = yr
	
		# Observed
		obs = dat[[i]][dat[[i]]$Year %in% yr,]$Catch
		yield_obs[,i] = if(cal) {mean(obs)} else{obs}
	}
	
	# Calculate and return relative error
	error = sum((log(yield_mod) - log(yield_obs))^2)
	return(error)
}

# A function to fit biomass time series with mortality parameters
#	par: input parameters c(U_crit, a_crit, U_met, a_met) in that order
#	params: mizer params object
#	dat: a list containing data frames with "Year" and "Catch" columns for cod, flounder, sprat, and herring, in that order
#	species: which species to fit
#	yr: years to which to fit
#	t: time steps of model to run
#	effort: a vector containing effort values for cod, flounder, sprat, and herring in that order
#	sim: mizer simulation object, if already run
#	returns sum of relative errors in all species' SSB 
time_yield = function(par, params, dat, species, yr, t, effort, sim = NULL) {
	
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
	
	# Set occupancy variables if unavailable
	if(!exists("occ.optim")) {occ.optim = list(par = c(10^5,0))}
	
	# Run simulation if not already provided
	if(is.null(sim)) {
		# Assign parameters
		if(species == "cod") {
			ret = oxy_sensitivity(ret, 
				U_hab = c(occ.optim$par[1], 10^5, -10^5, -10^5),
				a_hab = c(occ.optim$par[2], 0, 0, 0),
				U_search = c(us_cod, us_flounder, 10^5, 10^5),
				a_search = c(as_cod, as_flounder, 0, 0),
				U_maxin = c(ui_cod, ui_flounder, 10^5, 10^5),
				a_maxin = c(ai_cod, ai_flounder, 0, 0),
				U_erepro = c(ue_cod, ue_flounder, 10^5, 10^5),
				a_erepro = c(ae_cod, ae_flounder, 0, 0),
				U_alpha = c(ua_cod, ua_flounder, 10^5, 10^5),
				a_alpha = c(aa_cod, aa_flounder, 0, 0),
				U_met = c(um_cod, um_flounder, 10^5, 10^5),
				a_met = c(am_cod, am_flounder, -10^5, -10^5),
				z_mort = c(par[1], zm_flounder, 0, 0),
				b_mort = c(par[2], bm_flounder, 0, 0))
		} else if(species == "flounder") {
			ret = oxy_sensitivity(ret, 
				U_hab = c(occ.optim$par[1], 10^5, -10^5, -10^5),
				a_hab = c(occ.optim$par[2], 0, 0, 0),
				U_search = c(us_cod, us_flounder, 10^5, 10^5),
				a_search = c(as_cod, as_flounder, 0, 0),
				U_maxin = c(ui_cod, ui_flounder, 10^5, 10^5),
				a_maxin = c(ai_cod, ai_flounder, 0, 0),
				U_erepro = c(ue_cod, ue_flounder, 10^5, 10^5),
				a_erepro = c(ae_cod, ae_flounder, 0, 0),
				U_alpha = c(ua_cod, ua_flounder, 10^5, 10^5),
				a_alpha = c(aa_cod, aa_flounder, 0, 0),
				U_met = c(um_cod, um_flounder, 10^5, 10^5),
				a_met = c(am_cod, am_flounder, -10^5, -10^5),
				z_mort = c(zm_cod, par[1], 0, 0),
				b_mort = c(bm_cod, par[2], 0, 0))
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
	
		# Run to steady state
		sim = project(ret, t_max = t, effort = effort)
	}
	
	# Set up yield matrices
	yield_obs = matrix(nrow = length(yr), ncol = 4)
	rownames(yield_obs) = yr
	colnames(yield_obs) = c("cod", "flounder", "sprat", "herring")
	idx_mod = (t-length(yr)+1):t
	yield_mod = getYield(sim)[idx_mod, ]/1000/1000
	
	# Get yield for each species
	for(i in 1:ncol(yield_obs)) {
		
		# Use weights greater than w_mat
		idx_obs = yr
	
		# Observed
		obs = dat[[i]][dat[[i]]$Year %in% yr,]$Catch
		yield_obs[,i] = obs
	}
	
	# Calculate and return relative error
	error = ifelse(species == "cod", sum((log(yield_mod[,1]) - log(yield_obs[,1]))^2), sum((log(yield_mod[,2]) - log(yield_obs[,2]))^2))
	return(error)

}