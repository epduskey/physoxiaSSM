# Load oxygen data
load("Data/oxy_model.rda")

# A function to find ideal vs poor oxygen and run respective simulations
#	params: a Mizer params object
#	t: time steps to run model
#	p: is there physiological scaling
#	returns a list with mizer object, size at max obs. age, and diet
oxy_scenario = function(params, t, oxy_benthic, p = F) {
	
	# Create temporary params object
	ret = params
	
	# Create an oxygen vector
	times = 0:t
	bo = vector(mode = "numeric", length = length(times))
	po = vector(mode = "numeric", length = length(times))
	
	# Extend temperature vectors
	ret@other_params$benthic_temp = rep(ret@other_params$benthic_temp[1], t+1)
	ret@other_params$pelagic_temp = rep(ret@other_params$pelagic_temp[1], t+1)

	# Use 6 mL/L as pelagic oxygen
	oxy_pelagic = 6
	
	# Fill oxygen vectors
	bo[1:length(times)] = oxy_benthic
	po[1:length(times)] = oxy_pelagic
	
	# Store these in other params
	ret@other_params$benthic_oxygen = bo
	ret@other_params$pelagic_oxygen = po
	
	# Shorten temperature object if necessary
	ret@other_params$benthic_temp = ret@other_params$benthic_temp[1:length(bo)]
	ret@other_params$pelagic_temp = ret@other_params$pelagic_temp[1:length(po)]
	
	# Scale occupancy
	ret = occupancy(ret, t)
	
	# Scale objects
	if(p) {
		
		# Check whether params object has MP or not
		if(is.null(species_params(ret)$U_crit)) {
			
			# Scale physiological rates
			ret = rate_scale(ret, t)
			ret = mscale(ret, t)
		} else {
			
			# Scale physiological rates
			ret = rs_cphys(ret, t)
			ret = mscale(ret, t)
		}
	}
	if(!p) {
		
		# Check whether params object has MP or not
		if(is.null(species_params(ret)$U_crit)) {
			
			# Scale physiological rates
			ret = rate_scale(ret, t)
			ret = mscale(ret, t)
		
			# Make sure scaling is set to 1
			ret@other_params$rs_search[,,] = 1
			ret@other_params$rs_maxin[,,] = 1
			ret@other_params$rs_erepro[,,] = 1
			ret@other_params$rs_alpha[,,] = 1
			ret@other_params$metab_scale[,,] = 1
			
		} else {
			
			# Scale physiological rates
			ret = rs_cphys(ret, t)
			ret = mscale(ret, t)
			
			# Make sure scaling is set to 1
			ret@other_params$rate_scale[,,] = 1
			ret@other_params$metab_scale[,,] = 1
		}
	}
	
	# Run scenario
	sim = project(ret, t_max = t, effort = c(1,1,1,1))
	
	# Get weight categories and weight category widths
	w_all = sim@params@w
	dw_all = sim@params@dw
	
	# Use mature fish for diet
	idx_cod = w_all >= species_params(sim@params)["cod", ]$w_mat
	idx_flounder = w_all >= species_params(sim@params)["flounder", ]$w_mat
	idx_sprat = w_all >= species_params(sim@params)["sprat", ]$w_mat
	idx_herring = w_all >= species_params(sim@params)["herring", ]$w_mat
	
	# Get maintenance costs as a proportion for all species
	cod_metab = sim@params@metab["cod",idx_cod] / 
		rowSums(getHypoxiaDiet(sim,t+1,proportion=F)["cod",idx_cod,])
	flounder_metab = sim@params@metab["flounder",idx_flounder] / 
		rowSums(getHypoxiaDiet(sim,t_max+1,proportion=F)["flounder",idx_flounder,])
	sprat_metab = sim@params@metab["sprat",idx_sprat] / 
		rowSums(getHypoxiaDiet(sim,t_max+1,proportion=F)["sprat",idx_sprat,])
	herring_metab = sim@params@metab["herring",idx_herring] / 
		rowSums(getHypoxiaDiet(sim,t_max+1,proportion=F)["herring",idx_herring,])
		
	# Get weighted diet (g/year)
	cod = colSums(sweep(getHypoxiaDiet(sim, t_max+1, proportion = F)["cod",idx_cod,],1,-1*cod_metab+1,"*") * (sim@n[dim(sim@n)[1],,] * w_all * dw_all)["cod",idx_cod], na.rm = T) / sum((sim@n[dim(sim@n)[1],,] * w_all * dw_all)["cod",idx_cod])
	flounder = colSums(sweep(getHypoxiaDiet(sim, t_max+1, proportion = F)["flounder",idx_flounder,],1,-1*flounder_metab+1,"*") * (sim@n[dim(sim@n)[1],,] * w_all * dw_all)["flounder",idx_flounder], na.rm = T) / sum((sim@n[dim(sim@n)[1],,] * w_all * dw_all)["flounder",idx_flounder])
	sprat = colSums(sweep(getHypoxiaDiet(sim, t_max+1, proportion = F)["sprat",idx_sprat,],1,-1*sprat_metab+1,"*") * (sim@n[dim(sim@n)[1],,] * w_all * dw_all)["sprat",idx_sprat], na.rm = T) / sum((sim@n[dim(sim@n)[1],,] * w_all * dw_all)["sprat",idx_sprat])
	herring = colSums(sweep(getHypoxiaDiet(sim, t_max+1, proportion = F)["herring",idx_herring,],1,-1*herring_metab+1,"*") * (sim@n[dim(sim@n)[1],,] * w_all * dw_all)["herring",idx_herring], na.rm = T) / sum((sim@n[dim(sim@n)[1],,] * w_all * dw_all)["herring",idx_herring])
	
	# Get body size at maximum observed age
	size = c(cod = myGrowthCurves(sim, t = t_max, species = "cod", max_age = 15)[1,50], 
		flounder = myGrowthCurves(sim, t = t_max, species = "flounder", max_age = 26)[1,50], 
		sprat = myGrowthCurves(sim, t = t_max, species = "sprat", max_age = 16)[1,50], 
		herring = myGrowthCurves(sim, t = t_max, species = "herring", max_age = 13)[1,50])
	
	# Return scenarios
	return(list(sim = sim, 
		size = size, 
		cod = cod, 
		flounder = flounder, 
		sprat = sprat, 
		herring = herring))
}

# A function to get SSB
#	sim: a Mizer sim object
#	returns ssb for each species at the last time step of the model
get_ssb = function(sim) {

	# Set up SSB output matrix
	ssb = matrix(NA, nrow = 1, ncol = 4)
	colnames(ssb) = c("cod", "flounder", "sprat", "herring")
	n_all = sim@n[dim(sim@n)[1],,]
	
	# Get SSB for each species
	for(i in 1:ncol(ssb)) {
		
		# Use weights greater than w_mat
		idx_w = sim@params@w >= species_params(sim@params)$w_mat[i]
		
		# Model
		n = n_all[i,idx_w]
		mod = n * sim@params@w[idx_w] * sim@params@dw[idx_w]
		ssb[,i] = sum(mod)
	}
	
	return(ssb)
}