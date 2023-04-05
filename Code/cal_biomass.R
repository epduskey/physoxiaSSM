# Fitting biomass

# A function to fit biomass
#	par: input parameters c(all species' rmax exponents, benthic and pelagic carrying capacity exponents) in that order
#	params: mizer params object
#	dat: a list containing data frames with "Year" and "SSB" columns for cod, flounder, sprat, and herring, in that order
#	yr: years to which to fit
#	t: number of time steps for the model
#	effort: a vector containing effort values for cod, flounder, sprat, and herring in that order
#	cal: if true, takes mean values over all years
#	sim: mizer simulation object, if already run
#	returns sum of relative errors in all species' SSB 
fit_bg = function(par, params, dat, yr, t, effort, sim = NULL, cal = T) {
	
	# New params object
	ret = params
	
	# Run simulation if not already provided
	if(is.null(sim)) {

		# Assign parameters
		rmax = c(10^par[1], 10^par[2], 10^par[3], 10^par[4])
		gamma = c(10^par[5], 10^par[6], 10^par[7], 10^par[8])
	
		# Reset parameter values
		species_params(ret)$R_max = rmax
		species_params(ret)$gamma = gamma

		# Run to steady state
		sim = project(ret, t_max = t, effort = effort)
	}
	
	# Set up SSB matrices
	ssb_obs = matrix(nrow = ifelse(cal,1,length(yr)), ncol = 4)
	rownames(ssb_obs) = if(cal) {"cal"} else{yr}
	colnames(ssb_obs) = c("cod", "flounder", "sprat", "herring")
	ssb_mod = ssb_obs
	n_all = if(cal) {sim@n[dim(sim@n)[1],,]} else{sim@n[dim(sim@n)[1]:(dim(sim@n)[1]-length(yr)+1),,]}
	n_all = n_all/1000/1000/1000
	
	# Get SSB for each species
	for(i in 1:ncol(ssb_obs)) {
		
		# Use weights greater than w_mat
		idx_w = ret@w >= species_params(ret)$w_mat[i]
	
		# Observed
		obs = dat[[i]][dat[[i]]$Year %in% yr,]$SSB/1000/1000/1000
		ssb_obs[,i] = if(cal) {mean(obs)} else{obs}
		
		# Model
		n = if(cal) {n_all[i,idx_w]} else{n_all[,i,idx_w]}
		mod = if(cal) {n * ret@w[idx_w] * ret@dw[idx_w]} else{rowSums(sweep(n, 2, (ret@w*params@dw)[idx_w], "*"))}
		ssb_mod[,i] = if(cal) {sum(mod)} else{mod}
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
		sp = species_params(params_growth)
		sp_names = c("cod", "flounder", "sprat", "herring")
		sp_maxage = c(15,26,16,13)
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
		
		error_growth = sum((growth_mod - growth_obs)^2)
	} else {
		
		# Growth matrices
		growth_mod = array(dim = c(length(yr), 4, 50))
		species_all = c("cod","flounder","sprat","herring")
		dimnames(growth_mod) = list(Year = yr, Species = species_all, Age = seq(0,1,length.out=50))
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
			sp = species_params(params_growth)
			L_inf = (sp$w_inf/sp$a)^(1/sp$b)
			for(j in 1:dim(growth_obs)[2]) {
				growth_mod[i,j,] = getGrowthCurves(params_growth, species = rownames(sp)[j], max_age = sp_maxage[j], percentage = T)
				growth_obs[i,j,] = sp$a[j] * (L_inf[j] * (1 - exp(-sp$k_vb[j] * (as.numeric(dimnames(growth_mod)$Age)*sp_maxage[j])))) ^ sp$b[j]
				growth_obs[i,j,] = (growth_obs[i,j,]/sp$w_inf[j])*100
			}
		}
		
		error_growth = sum((growth_mod - growth_obs)^2)
	}
	
	# Calculate and return relative error
	error_ssb = sum((ssb_mod - ssb_obs)^2)
	return(0.01*error_growth + error_ssb)
}