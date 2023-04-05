# Plot emergent size spectra
#	sim: mizer sim object
#	t: time step for plotting
#	returns nothing, saves plot if users turn cal on
plot_spectra = function(sim, t, cal = F) {
	
	# Get params object
	params = sim@params
	
	# Index for weights at 1e-5 or greater
	w_idx = params@w_full >= 1e-5

	# Weight cutoffs
	w_idx_benthic = params@w_full <= params@resource_params$w_pp_cutoff
	w_idx_pelagic = params@w_full <= params@resource_params$pelagic$w_pp_cutoff

	# Pelagic spectra
	w_pp_pelagic = aperm(array(unlist(sim@n_other), dim = rev(dim(sim@n_pp))))
	w_pelagic = params@w_full[w_idx & w_idx_pelagic]
	w_pp_pelagic_stable = w_pp_pelagic[t, w_idx & w_idx_pelagic] * w_pelagic

	# Benthic spectra
	w_benthic = params@w_full[w_idx & w_idx_benthic]
	w_pp_benthic_stable = sim@n_pp[t, w_idx & w_idx_benthic] * w_benthic

	# Species spectra
	w_species = t(t(sim@n[t,,]) * params@w)
	
	# Plot results
	ymax = max(c(w_pp_pelagic_stable, w_pp_benthic_stable))
	if(cal) {jpeg("Plots/yield_spectra.jpeg", width = 6324, height = 6324, units = 'px', res = 600)}
	plot(w_benthic, w_pp_benthic_stable, type = 'l', lwd = 4, col = 'green', 
			log = 'xy', xlim = c(1e-5,1e+5), ylim = c(1e-5,ymax), xlab = "Size (g)", ylab = "Biomass density")
	lines(w_pelagic, w_pp_pelagic_stable, lwd = 4, col = 'darkgreen')
	for(i in 1:nrow(w_species)) {lines(params@w, w_species[i, ], lwd = 4, col = i+1)}
	legend("topright", bty = 'n', c("benthic","pelagic","cod","flounder","sprat","herring"), lwd = 4, col = c("green","darkgreen",2,3,4,5))
	par(mfrow = c(1,1))
	if(cal) {dev.off()}
}

# Plot emergent size spectra in a conceptual style with drawing colors
#	sim: mizer sim object
#	t: time step for plotting
concept_spectra = function(sim, t) {
	
	# Get params object
	params = sim@params
	
	# Index for weights at 1e-5 or greater
	w_idx = params@w_full >= 1e-5

	# Weight cutoffs
	w_idx_benthic = params@w_full <= params@resource_params$w_pp_cutoff
	w_idx_pelagic = params@w_full <= params@resource_params$pelagic$w_pp_cutoff

	# Pelagic spectra
	w_pp_pelagic = aperm(array(unlist(sim@n_other), dim = rev(dim(sim@n_pp))))
	w_pelagic = params@w_full[w_idx & w_idx_pelagic]
	w_pp_pelagic_stable = w_pp_pelagic[t, w_idx & w_idx_pelagic] * w_pelagic

	# Benthic spectra
	w_benthic = params@w_full[w_idx & w_idx_benthic]
	w_pp_benthic_stable = sim@n_pp[t, w_idx & w_idx_benthic] * w_benthic

	# Species spectra
	w_species = t(t(sim@n[t,,]) * params@w)
	
	# Plot results
	color_benthos = rgb(38, 24, 95, maxColorValue=255)
	color_pelagic = rgb(0, 106, 168, maxColorValue = 255)
	color_cod = rgb(0, 166, 174, maxColorValue = 255)
	color_flounder = rgb(252, 255, 221, maxColorValue = 255)
	color_sprat = rgb(205, 240, 203, maxColorValue = 255)
	color_herring = rgb(119, 209, 181, maxColorValue = 255)
	ymax = max(c(w_pp_pelagic_stable, w_pp_benthic_stable))
	ymin = min(w_species)
	par(mar = c(3,3,3,1))
	plot(w_benthic, w_pp_benthic_stable, type = 'l', lwd = 10, col = 1, axes = F, main = "Size Spectra", cex.main = 2,
			log = 'xy', xlim = c(1e-5,1e+5), ylim = c(ymin+1e-1,ymax), xlab = "", ylab = "")
	box()
	mtext("log(Weight)", side = 1, cex = 1.5, line = 0.5)
	mtext("log(Biomass)", side = 2, cex = 1.5, line = 0.5)
	lines(w_benthic, w_pp_benthic_stable, lwd = 5, col = color_benthos)
	rasterImage(sad_drawing, 10^-5.2,10^13,10^-4.5,10^16)
	lines(w_pelagic, w_pp_pelagic_stable, lwd = 10, col = 1)
	lines(w_pelagic, w_pp_pelagic_stable, lwd = 5, col = color_pelagic)
	rasterImage(cop_drawing, 10^-2.5,10^13,10^-1.5,10^16)
	lines(params@w, w_species[1, ], lwd = 10, col = 1)
	lines(params@w, w_species[1, ], lwd = 5, col = color_cod)
	rasterImage(cod_drawing, 10^0,10^2,10^5,10^8)
	lines(params@w, w_species[2, ], lwd = 10, col = 1)
	lines(params@w, w_species[2, ], lwd = 5, col = color_flounder)
	rasterImage(fln_drawing, 10^1,10^-1,10^3,10^3)
	lines(params@w, w_species[3, ], lwd = 10, col = 1)
	lines(params@w, w_species[3, ], lwd = 5, col = color_sprat)
	rasterImage(spt_drawing, 10^1,10^11,10^2,10^12)
	lines(params@w, w_species[4, ], lwd = 10, col = 1)
	lines(params@w, w_species[4, ], lwd = 5, col = color_herring)
	rasterImage(her_drawing, 10^2.2,10^8.6,10^3.75,10^10.5)
}