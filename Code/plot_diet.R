# Load array binding packages
library(abind)
library(ggplot2)

# My diet retrieval function
#	sim: mizer sim object
#	t: model time step
#	proportion: returns diet as a proportion if true
#	returns consumer x body size x prey array
getHypoxiaDiet = function(sim, t, proportion = T) {
	
	# Get abundance at time t 
	n = sim@n[t,,]
	n_pp = sim@n_pp[t,]
	n_other = list(n_pp_pelagic = sim@n_other[[t]])

	# Set up diet array
	params = validParams(sim@params)
	species = params@species_params$species
	no_sp = length(species)
	no_w = length(params@w)
	no_w_full = length(params@w_full)
	no_other = length(params@initial_n_other)
	assert_that(identical(dim(n), c(no_sp, no_w)), length(n_pp) == no_w_full)
	diet = array(0, dim = c(no_sp, no_w, no_sp + 1 + no_other), 
		dimnames = list(predator = species, w = dimnames(params@initial_n)$w, prey = c(as.character(species), "Benthic", "Pelagic")))
	idx_sp = (no_w_full - no_w + 1):no_w_full
	
	# Get hypothetical diet based on predation kernel
	if (!is.null(comment(params@pred_kernel))) {
		ae = matrix(params@pred_kernel[, , idx_sp, drop = FALSE], ncol = no_w) %*% t(sweep(n, 2, params@w * params@dw, "*"))
		diet[,,1:no_sp] = ae
		diet[,,no_sp+1] = rowSums(sweep(params@pred_kernel, 3, params@dw_full * params@w_full * n_pp, "*"), dims = 2)
		diet[,,no_sp+1+1] = rowSums(sweep(params@pred_kernel, 3, params@dw_full * params@w_full * n_other$n_pp_pelagic, "*"), dims = 2)
	} else {
		prey_benthic = matrix(0, nrow = no_sp+1+1, ncol = no_w_full)
		prey_pelagic = prey_benthic
		prey_benthic[1:no_sp, idx_sp] = sweep(n*params@other_params$benthic_occupancy[t,,], 2, params@w * params@dw, "*")
		prey_pelagic[1:no_sp, idx_sp] = sweep(n*(1-params@other_params$benthic_occupancy[t,,]), 2, params@w * params@dw, "*")
		prey_benthic[no_sp+1, ] = n_pp * params@w_full * params@dw_full
		prey_pelagic[no_sp+1+1, ] = n_other$n_pp_pelagic * params@w_full * params@dw_full
		ft_benthic = array(rep(params@ft_pred_kernel_e, times = no_sp+1+1) * rep(mvfft(t(prey_benthic)), each = no_sp), dim = c(no_sp, no_w_full, no_sp+1+1))
		ft_pelagic = array(rep(params@ft_pred_kernel_e, times = no_sp+1+1) * rep(mvfft(t(prey_pelagic)), each = no_sp), dim = c(no_sp, no_w_full, no_sp+1+1))
		ft_benthic = matrix(aperm(ft_benthic, c(2,1,3)), nrow = no_w_full)
		ft_pelagic = matrix(aperm(ft_pelagic, c(2,1,3)), nrow = no_w_full)
		ae_benthic = array(Re(mvfft(ft_benthic, inverse = TRUE)/no_w_full), dim = c(no_w_full, no_sp, no_sp+1+1))
		ae_pelagic = array(Re(mvfft(ft_pelagic, inverse = TRUE)/no_w_full), dim = c(no_w_full, no_sp, no_sp+1+1))
		ae_benthic = ae_benthic[idx_sp, , , drop = FALSE]
		ae_pelagic = ae_pelagic[idx_sp, , , drop = FALSE]
		ae_benthic = aperm(ae_benthic, c(2,1,3))
		ae_pelagic = aperm(ae_pelagic, c(2,1,3))
		ae_benthic[ae_benthic < 1e-18] = 0
		ae_pelagic[ae_pelagic < 1e-18] = 0
		ae_benthic[,,1:no_sp] = sweep(ae_benthic[,,1:no_sp,drop=FALSE], c(2,3), t(params@other_params$benthic_occupancy[t,,]), "*")
		ae_pelagic[,,1:no_sp] = sweep(ae_pelagic[,,1:no_sp,drop=FALSE], c(2,3), t(1-params@other_params$benthic_occupancy[t,,]), "*")
		diet[, , 1:(no_sp+1+1)] = sweep(ae_benthic, c(1,2), params@other_params$benthic_occupancy[t,,], "*") + sweep(ae_pelagic, c(1,2), (1-params@other_params$benthic_occupancy[t,,]), "*")
	}
	
	# Get search volume scaling
	if(is.null(species_params(params)$U_crit)) {
		U_search = params@other_params$rs_search[t,,]
	} else {
		U_search = params@other_params$rate_scale[t,,]
	}
	
	# Add interactions and other encounters
	inter = cbind(params@interaction, params@species_params$interaction_resource, params@species_params$interaction_resource)
	diet[,,1:(no_sp+1+1)] = sweep(sweep(diet[,,1:(no_sp+1+1),drop = FALSE], c(1,3), inter, "*"), c(1,2), U_search*params@search_vol, "*")
	for (i in seq_along(params@other_encounter)) {
		diet[,,no_sp+1+i] = do.call(params@other_encounter[[i]], 
            					list(params = params, 
							n = n, 
							n_pp = n_pp, 
							n_other = n_other, 
							component = names(params@other_encounter)[[i]], 
							t = t-1))
	}
	
	# Include feeding level
	f = getFeedingLevel(sim, n, n_pp, time_range = t-1)[1,,]
	fish_mask = n > 0
	diet = sweep(diet, c(1, 2), (1 - f) * fish_mask, "*")
	
	# Convert to proportion
	if (proportion) {
		total = rowSums(diet, dims = 2)
		diet = sweep(diet, c(1, 2), total, "/")
		diet[is.nan(diet)] <- 0
	}
	
	return(diet)
}

# Plot diet proportion
#	sim: mizer sim object
#	species: one species to plot
#	t: model time step
#	return_data: returns diet data if true
#	cal: saves plot if users turn on
plotHypoxiaDiet = function(sim, species, t, return_data = F, cal = F) {
	
	# Check input
	assert_that(is.flag(return_data))
 	params = validParams(sim@params)
	if(length(species) != 1) {stop("just choose one species, don't get greedy")}
	lab = species
	
	# Get diet
	species = valid_species_arg(sim, species, return.logical = TRUE)
    	diet = getHypoxiaDiet(sim, t)[species,,]
	colnames(diet) = c("cod", "flounder", "sprat", "herring", "benthos", "plankton")
	
	# Construct diet data frame
    	prey = dimnames(diet)$prey
    	prey = factor(prey, levels = rev(prey))
    	plot_dat = data.frame(w = params@w, Proportion = c(diet), Prey = rep(prey, each = length(params@w)))
	plot_dat = plot_dat[plot_dat$Proportion > 0.001, ]
	if (return_data) {return(plot_dat)}
	
	# Plot colors
	color_benthos = rgb(38, 24, 95, maxColorValue=255)
	color_pelagic = rgb(0, 106, 168, maxColorValue = 255)
	color_cod = rgb(0, 166, 174, maxColorValue = 255)
	color_flounder = rgb(252, 255, 221, maxColorValue = 255)
	color_sprat = rgb(205, 240, 203, maxColorValue = 255)
	color_herring = rgb(119, 209, 181, maxColorValue = 255)
	
	# Plot diet
	if(cal) {jpeg(paste("Plots/diet_calibraiton_", lab, ".jpeg", sep=""), width = 6324, height = 6324, units = 'px', res = 600)}
    	legend_colors = c(benthos=color_benthos, plankton=color_pelagic, cod=color_cod, flounder=color_flounder, sprat=color_sprat, herring=color_herring)
	legend_which = intersect(plot_dat$Prey, names(legend_colors))
    	g = ggplot(plot_dat) + 
		geom_area(aes(x = w, y = Proportion, fill = Prey)) + 
		scale_x_log10() + labs(x = "Size [g]") + scale_fill_manual(values = legend_colors[legend_which])
	plot(g)
	if(cal) {dev.off()}
}