# Create yearly effort array

# Load fishing mortality data
cod_fyear = read.table("Parameters/Cod/cod_fyear.txt", header = T)
flounder_fyear = read.table("Parameters/Flounder/flounder_fyear.txt", header = T)
sprat_fyear = read.table("Parameters/Sprat/sprat_fyear.txt", header = T)
herring_fyear = read.table("Parameters/Herring/herring_fyear.txt", header = T)
fdat = list(cod = cod_fyear, flounder = flounder_fyear, sprat = sprat_fyear, herring = herring_fyear)

# Function to return time x gear fishing mortality array
#	params: mizer params object
#	dat: a named list with fishing mortality data frames for each species, in the same order as species appearing in params
#	yr: years for data
#	fcol: a vector containing the column in which fishing mortality appears
#	t: the number of time steps
#	returns a time x gear array with effort values to be used in mizer project
farray = function(params, dat, yr, fcol, t) {
	
	# Check data
	if(sum(rownames(species_params(params)) != names(dat)) != 0) {
		stop("dat must be named list with entries for each species in params and in the same order")
	}
	
	# Create return array
	ret = array(1, dim = c(t, length(dat)))
	dimnames(ret) = list(time = seq(t), gear = rownames(species_params(params)))
	
	# Set last years to data values
	for(i in 1:ncol(ret)) {
		
		ret[(t-nrow(dat[[i]][dat[[i]]$Year %in% yr, ])+1):t,i] = dat[[i]][dat[[i]]$Year %in% yr,fcol[i]] / species_params(params)$catchability[i]
	}
	
	return(ret)
}
