# Automated calibration code

# Load benthic occupancy data for cod
cod_benthic = read.table("Parameters/Cod/cod_benthic_occupancy.txt", header = T)

# Check benthic occupancy
#	par: the values c(u,a) of f(x) = 1/(1 + exp(-u*(x-a*p_crit))), in that order
#	params: mizer params object
#	dat: occupancy data against which to fit
#	t: time steps of model to run
#	yr: calibration years
#	effort: a vector containing effort values for cod, flounder, sprat, and herring in that order
#	returns sum of squared error
bocc_optim = function(par, params, dat, t, yr, effort) {
	
	# New params object
	ret = params
	
	# Input parameter
	u = par[1]
	a = par[2]
	
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
	
	# Set up params object
	ret = oxy_sensitivity(ret, 
				U_hab = c(u, 10^5, -10^5, -10^5),
				a_hab = c(a, 0, 0, 0),
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
				z_mort = c(zm_cod, zm_flounder, 0, 0),
				b_mort = c(bm_cod, bm_flounder, 0, 0))
	
	# Get occupancy values
	ret@other_params$benthic_oxygen = c(rep(mean(preds_all$df[preds_all$df$Year %in% yr, ]$Oxygen_benthic),t-length(yr)+1), preds_all$df[preds_all$df$Year %in% yr, ]$Oxygen_benthic)
	ret@other_params$pelagic_oxygen = c(rep(mean(preds_all$df[preds_all$df$Year %in% yr, ]$Oxygen_pelagic),t-length(yr)+1), preds_all$df[preds_all$df$Year %in% yr, ]$Oxygen_pelagic)
	ret = occupancy(ret, t)
	
	# Run to steady state
	sim = project(ret, t_max = t, effort = effort)
	
	# Set up params object
	ret = sim@params
	ret@initial_n[] <- sim@n[dim(sim@n)[1],, ]
    	ret@initial_n_pp[] <- sim@n_pp[dim(sim@n)[1], ]
    	ret@initial_n_other[] <- sim@n_other[dim(sim@n)[1], ]
    	ret@initial_effort[] <- sim@effort[dim(sim@n)[1], ]
    	ret@time_modified <- lubridate::now()
	
	# Calculate occupancy
	bocc = ret@other_params$benthic_occupancy[,"cod",params@w >= cod_lh["a","stan"]*species_params(ret)$l50[1]^cod_lh["b","stan"]]
	bocc_wm = rowSums(sweep(bocc, 2, initialN(ret)["cod",params@w >= cod_lh["a","stan"]*species_params(ret)$l50[1]^cod_lh["b","stan"]], "*")) / sum(initialN(ret)["cod",params@w >= cod_lh["a","stan"]*species_params(ret)$l50[1]^cod_lh["b","stan"]])
	
	# Calculate and return error
	error = sum((tail(bocc_wm,length(yr)) - dat[dat$Year %in% yr, ]$Benthic_occ)^2)
	return(error)
}