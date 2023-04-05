# Plotting pcrit

# Plot colors
color_cod = rgb(0, 166, 174, maxColorValue = 255)
color_flounder = rgb(252, 255, 221, maxColorValue = 255)
color_sprat = rgb(205, 240, 203, maxColorValue = 255)
color_herring = rgb(119, 209, 181, maxColorValue = 255)
col_all = c(color_cod, color_flounder, color_sprat, color_herring)

# Plots pcrit as a function of body size
#	params: mizer params object
#	jpg: save to jpeg if T, default is F
#	returns nothing, just plots
plot_pcrit = function(params, jpg = F) {
	if(jpg) jpeg("Plots/pcrit.jpeg", width = 20, height = 13, units = 'cm', res = 600)
	par(mar = c(5,5,3,1))
	plot(0, type = 'n', xlim = c(1e-3,15000), ylim = c(1,3), axes = F, cex.lab = 2, log = 'x', xlab = "Weight (g)", ylab = expression(paste(P[crit], " ", (mL/L))))
	axis(1, at = 10^seq(-3,4), labels = 10^seq(-3,4), cex.axis = 1.5)
	axis(2, cex.axis = 1.5)
	box()
	for(i in 1:nrow(species_params(params))) {
		lines(params@w[params@w < species_params(params)$w_inf[i]], params@other_params$P_crit[i, params@w < species_params(params)$w_inf[i]], lwd = 10)
		lines(params@w[params@w < species_params(params)$w_inf[i]], params@other_params$P_crit[i, params@w < species_params(params)$w_inf[i]], lwd = 5, col = col_all[i])
	}
	legend("topleft", bty = 'n', c("", "", "", ""), lwd = 10, col = "black")
	legend("topleft", bty = 'n', c("Cod", "Flounder", "Sprat", "Herring"), lwd = 5, col = col_all)
	if(jpg) dev.off()
}