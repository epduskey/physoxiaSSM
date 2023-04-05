plotFmsy = function(params, effort_res = 40, max_f = 2, t_max = 100) {
	
	# Set new params object
	params_temp = params
	
	# Set catchability to 1 for all species
	catchability(params_temp) = catchability(params_temp) * diag(c(1,1,1,1))
	species_params(params_temp)$catchability = diag(catchability(params_temp))
	
	# Set up data frame
	f = seq(0, max_f, length.out=effort_res)
	temp = data.frame(F=f)
	temp$cod = rep(NA, nrow(temp))
	temp$flounder = rep(NA, nrow(temp))
	temp$sprat = rep(NA, nrow(temp))
	temp$herring = rep(NA, nrow(temp))
	
	# Loop through all fishing mortalities
	for(i in 1:length(f)) {
		sim_temp = projectToSteady(params_temp, t_max = t_max, effort = f[i], return_sim = T)
		y = getYield(sim_temp)
		temp$cod[i] = y[dim(y)[1],"cod"]
		temp$flounder[i] = y[dim(y)[1],"flounder"]
		temp$sprat[i] = y[dim(y)[1],"sprat"]
		temp$herring[i] = y[dim(y)[1],"herring"]
	}
	
	# Plot results
	par(mfrow = c(2,2), mar = c(5,3,3,1), oma = c(2,3,3,1))
	plot(cod/1000/1000/1000 ~ F, temp, xlab = "", ylab = "", main = "cod")
	plot(flounder/1000/1000/1000 ~ F, temp, xlab = "", ylab = "", main = "flounder")
	plot(sprat/1000/1000/1000 ~ F, temp, xlab = "", ylab = "", main = "sprat")
	plot(herring/1000/1000/1000 ~ F, temp, xlab = "", ylab = "", main = "herring")
	mtext("F", side = 1, cex = 1.5, outer = T)
	mtext("Yield", side = 2, cex = 1.5, outer = T)
	par(mfrow = c(1,1), c(5.1,4.1,4.1,2.1))
}