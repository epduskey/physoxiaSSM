# Load growth data
cod_smalk = rbind(read.table("Data/cod_ivb.txt"),read.table("Data/cod_ivb_sim.txt"))
flounder_smalk = rbind(read.table("Data/flounder_ivb.txt"),read.table("Data/flounder_ivb_sim.txt"))
sprat_bias = rbind(read.table("Data/sprat_ivb.txt"),read.table("Data/sprat_ivb_sim.txt"))
herring_bias = rbind(read.table("Data/herring_ivb.txt"),read.table("Data/herring_ivb_sim.txt"))

# Plot emergent growth against observations
#	object: mizer sim or params object
#	t: time step for plotting
#	yr: years over which to summarize growth
#	cal: save calibration plot results
#	returns growth, saves if user turns cal on
plotGrowthObservedVsModel = function(object, t, yr, cal = F) {
	
	# Get mizer params object
	if(is(object, "MizerSim")) {
		params_temp = object@params
	} else if(is(object, "MizerParams")) {
		params_temp = validParams(object)
	}
	
	# Get length-weight conversion parameters
	a = params_temp@species_params$a
	b = params_temp@species_params$b
	
	# Subset data by year
	cod_dat = cod_smalk[cod_smalk$Year %in% yr, ]
	flounder_dat = flounder_smalk[flounder_smalk$Year %in% yr, ]
	if(dim(flounder_dat)[1] == 0) {flounder_dat = flounder_smalk}
	sprat_dat = sprat_bias[sprat_bias$Year %in% yr, ]
	herring_dat = herring_bias[herring_bias$Year %in% yr, ]
	
	# Adjust age
	cod_dat$Age = cod_dat$Age + 0.5
	flounder_dat$Age = flounder_dat$Age + 0.5
	sprat_dat$Age = sprat_dat$Age + 0.5
	herring_dat$Age = herring_dat$Age + 0.5
	
	# Get growth curves and ages
	gr_mod = array(NA, dim = c(4, 50))
	sp = species_params(params_temp)
	dimnames(gr_mod) = list(sp = rownames(sp), Age = seq(0,1,length.out=50))
	gr_mod[1,] = getGrowthCurves(object, species = "cod", max_age = max_age[1])
	gr_mod[2,] = getGrowthCurves(object, species = "flounder", max_age = max_age[2])
	gr_mod[3,] = getGrowthCurves(object, species = "sprat", max_age = max_age[3])
	gr_mod[4,] = getGrowthCurves(object, species = "herring", max_age = max_age[4])
	age = as.numeric(colnames(gr_mod))
	Linf = (sp$w_inf/sp$a)^(1/sp$b)
	gr_obs = sweep(sweep(sweep(1 - exp(sweep(matrix(max_age, ncol = 1) %*% age, 1, -sp$k_vb, "*")), 1, Linf, "*"), 1, sp$b, "^"), 1, sp$a, "*")
	dimnames(gr_obs) = dimnames(gr_mod)
	
	# Plot
	if(cal) {jpeg("Plots/growth_calibration.jpeg", width = 6324, height = 6324, units = 'px', res = 600)}
	par(mfrow = c(2,2), oma = c(5,3,3,1))
	plot(IndWgt ~ Age, data = cod_dat, pch = 16, cex = 0.5, cex.axis = 1.5, cex.main = 2, xlab = "", ylab = "", main = "Cod")
	lines(age*max_age[1], gr_mod["cod",], lwd = 2, col = "red")
	lines(age*max_age[1], gr_obs["cod",], lwd = 2, lty = 2, col = "red")
	plot(IndWgt ~ Age, data = flounder_dat, pch = 16, cex = 0.5, cex.axis = 1.5, cex.main = 2, xlab = "", ylab = "", main = "Flounder")
	lines(age*max_age[2], gr_mod["flounder",], lwd = 2, col = "red")
	lines(age*max_age[2], gr_obs["flounder",], lwd = 2, lty = 2, col = "red")
	plot(IndWgt ~ Age, data = sprat_dat, pch = 16, cex = 0.5, cex.axis = 1.5, cex.main = 2, xlab = "", ylab = "", main = "Sprat")
	lines(max_age[3]*age, gr_mod["sprat",], lwd = 2, col = "red")
	lines(age*max_age[3], gr_obs["sprat",], lwd = 2, lty = 2, col = "red")
	plot(IndWgt ~ Age, data = herring_dat, pch = 16, cex = 0.5, cex.axis = 1.5, cex.main = 2, xlab = "", ylab = "", main = "Herring")
	lines(max_age[4]*age, gr_mod["herring",], lwd = 2, col = "red")
	lines(age*max_age[4], gr_obs["herring",], lwd = 2, lty = 2, col = "red")
	mtext("Age", 1, cex = 2, outer = T)
	mtext("Weight (g)", 2, cex = 2, outer = T)
	par(mfrow = c(1,1), c(5.1,4.1,4.1,2.1))
	if(cal) {dev.off()}
	
	return(gr)
}

# Modified growth curves to look at growth in a given year
#	object: a mizer params or sim object
#	t: model time step in which to get growth
#	species: species to get
#	max_age: maximum age to which to calculate growth
#	percentage: if true, returns percentage of maximum size
#	returns estimated size at age at time step t
myGrowthCurves = function(object, t, species = NULL, max_age, percentage = FALSE) {

	# Check for valid object
	if (is(object, "MizerSim")) {
		params = object@params
		params@initial_n[] = object@n[t, , ]
    		params@initial_n_pp[] = object@n_pp[t, ]
    		params@initial_n_other[] = object@n_other[t, ]
    		params@initial_effort[] = object@effort[t, ]
    		params@time_modified = lubridate::now()
	} else if (is(object, "MizerParams")) {
		params = validParams(object)
	} else {
		stop("The first argument to `getGrowthCurves()` must be a MizerParams or a MizerSim object.")
	}
	
	# Get species
	species = valid_species_arg(params, species)
	idx = which(params@species_params$species %in% species)
	species = params@species_params$species[idx]
	
	# Set up growth array
	age = seq(0, max_age, length.out = 50)
	ws = array(dim = c(length(species), length(age)), dimnames = list(Species = species, Age = age))
	g = getEGrowth(params)
	
	# Full up array
	for (j in seq_along(species)) {
		i = idx[j]
		g_fn = stats::approxfun(c(params@w, params@species_params$w_inf[[i]]), c(g[i, ], 0))
		myodefun = function(t, state, parameters) {
			return(list(g_fn(state)))
		}
		ws[j, ] = deSolve::ode(y = params@w[params@w_min_idx[i]], times = age, func = myodefun)[, 2]
		if (percentage) {
			ws[j, ] = ws[j, ]/params@species_params$w_inf[i] * 100
		}
	}
	return(ws)
}