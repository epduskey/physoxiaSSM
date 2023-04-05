# Load catch data and estimates
cod_catch = read.csv("Data/cod_catch.csv")
flounder_catch = read.csv("Data/flounder_catch.csv")
sprat_catch = read.csv("Data/sprat_catch.csv")
herring_catch = read.csv("Data/herring_catch.csv")

# Plot emergent yield against observations
#	sim: mizer simulation object
#	t: time step for plotting
#	yr: years over which to summarize catch
#	cal: save calibration plot results
#	returns nothing, saves if users turn cal on
plotYieldObservedVsModel = function(sim, t, yr, cal = F) {
	
	# Summarize catch data over yr in tonnes
	cod_dat = mean(cod_catch[cod_catch$Year %in% yr, ]$Catch)
	flounder_dat = mean(flounder_catch[flounder_catch$Year %in% yr, ]$Catch)
	sprat_dat = mean(sprat_catch[sprat_catch$Year %in% yr, ]$Catch)
	herring_dat = mean(herring_catch[herring_catch$Year %in% yr, ]$Catch)
	
	# Create yield data frame
	yield = data.frame(sp = c("cod","flounder","sprat","herring"))
	yield$obs = c(cod_dat, flounder_dat, sprat_dat, herring_dat)
	yield$mod = getYield(sim)[t,]/1000/1000
	
	# Create yield matrix
	ymat = t(as.matrix(yield[,c(2,3)]))
	rownames(ymat) = c("observed","model")
	colnames(ymat) = yield$sp
	
	# Plot results
	if(cal) {jpeg("Plots/yield_calibraiton.jpeg", width = 6324, height = 6324, units = 'px', res = 600)}
	par(mar = c(5.1,4.1,4.1,5.3))
	barplot(ymat/1000, names.arg = yield[,1], legend.text = T, beside = T, ylab = "Yield (kt)", cex.axis = 1.2, cex.names = 1.5, cex.lab = 1.5, args.legend = list(x = "topright", inset = c(-0.25,0), bty = "n", cex = 1.0))
	box()
	if(cal) {dev.off()}
	
	return(yield)
}