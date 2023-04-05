# Scale occupancy and rates with oxygen and temperature

# Created: March 3, 2021
# Last modified: September 5, 2021 by EPD

# Contents (ctrl-f):
#	I. General scaling
#	II. Scale occupancy
#	III. Scale physiological rates


########## I. General scaling ##########

# Logistic scaling
#	U: determines slope at k; greater => greater slope i.e. more precipitous decline at k
#	k: pcrit
#	a: adjusts scaling at pcrit
#	oxy: one or more values for oxygen
#	returns logistic scaling
scaling = function(U, k, a, oxy) {
	return(1/(1+exp(-U*(oxy - a*k))))
}


########## II. Temperature scaling ##########

# Temperature scaling
#	m: body size
#	Temp: temperature in C
#	Tref: reference temperature in C
#	cA: activation constant
#	cD: deactivation constant
#	EA: activation energy
#	ED: deactivation energy
#	k: activation rate
#	TD: deactivation temperature
#	returns temperature scaling
tscale = function(m, Temp, Tref, cA, cD, EA, ED, k, TD) {
	size_act = m^(cA * (Temp - Tref))
	act = exp(-(EA/k) * ((1/Temp) - (1/Tref)))
	size_deact = m^(-cD * (Temp - TD))
	deact = (1 + exp(-(ED/k) * ((1/Temp) - (1/TD))))^(-1)
	refscale = (m^(cD * (Tref - TD))) * (1 + exp(-(ED/k) * ((1/Tref) - (1/TD))))
	ts = size_act * act * size_deact * deact * refscale
	
	return(ts)
}


########## III. Scale occupancy ##########

# Occupancy function
#	params: a MizerParams object which includes oxygen sensitivity values
#	oxy: a vector of oxygen values, one for each year the model runs
#	returns params object with added benthic occupancy to your params object
occupancy = function(params, t_max) {
	
	# Create times vector
	times = seq(0,t_max)
	
	# Create return object
	ret = params
	
	# Create empty array
	benthic_occupancy = array(NA, dim = c(length(times), nrow(ret@species_params), length(ret@w)))
	dimnames(benthic_occupancy) = list(time = times, sp = dimnames(ret@other_params$P_crit)$sp, w = dimnames(ret@other_params$P_crit)$w)
	
	# Calculate habitat scaling
	for(i in 1:dim(benthic_occupancy)[1]) {
		for(j in 1:dim(benthic_occupancy)[2]) {
			benthic_occupancy[i,j,] = scaling(U = ret@species_params[j,"U_hab"], k = ret@other_params$P_crit[j,], a = ret@species_params[j,"a_hab"], oxy = ret@other_params$benthic_oxygen[i])
		}
	}
	
	# Store oxygen
	ret@other_params$benthic_occupancy = benthic_occupancy
		
	return(ret)
	
}


########## IV. Scale physiological rates ##########

# Rate scaling function
#	params: a MizerParams object which includes oxygen sensitivity values
#	returns params object with added mean oxygen exposure and rate scaling arrays to your params object
rate_scale = function(params, t_max) {
	
	# Create times vector
	times = seq(0,t_max)
	
	# Create return object
	ret = params
	
	# Create empty array
	rs = array(NA, dim = c(length(times), nrow(params@species_params), length(params@w)))
	dimnames(rs) = list(time = times, sp = dimnames(params@other_params$P_crit)$sp, w = dimnames(params@other_params$P_crit)$w)
	
	# Create all physiological scaling arrays
	rs_search = rs
	rs_maxin = rs
	rs_erepro = rs
	rs_alpha = rs
	
	# Calculate oxygen exposure
	oxy_exposure = sweep(params@other_params$benthic_occupancy, 1, params@other_params$benthic_oxygen, '*') + 
				sweep(1 - params@other_params$benthic_occupancy, 1, params@other_params$pelagic_oxygen, '*')
	
	# Calculate temperature exposure
	temp_exposure = sweep(params@other_params$benthic_occupancy, 1, params@other_params$benthic_temp, '*') + 
				sweep(1 - params@other_params$benthic_occupancy, 1, params@other_params$pelagic_temp, '*')
	
	# Calculate rate scaling
	for(i in 1:dim(rs)[1]) {
		for(j in 1:dim(rs)[3]) {
			ss = scaling(U = params@species_params[,"U_search"], k = params@other_params$P_crit[,j], a = params@species_params[,"a_search"], oxy = oxy_exposure[i,,j])
			ms = scaling(U = params@species_params[,"U_maxin"], k = params@other_params$P_crit[,j], a = params@species_params[,"a_maxin"], oxy = oxy_exposure[i,,j])
			es = scaling(U = params@species_params[,"U_erepro"], k = params@other_params$P_crit[,j], a = params@species_params[,"a_erepro"], oxy = oxy_exposure[i,,j])
			as = scaling(U = params@species_params[,"U_alpha"], k = params@other_params$P_crit[,j], a = params@species_params[,"a_alpha"], oxy = oxy_exposure[i,,j])
			ts = tscale(m = params@w[j], 
					Temp = temp_exposure[i,,j], 
					Tref = params@other_params$Tref, 
					cA = params@species_params[,"cA"], 
					cD = params@species_params[,"cD"], 
					EA = params@species_params[,"EA"], 
					ED = params@species_params[,"ED"], 
					k = params@species_params[,"k"], 
					TD = params@species_params[,"TD"])
			rs_search[i,,j] = ss * ts
			rs_maxin[i,,j] = ms * ts
			rs_erepro[i,,j] = es * ts
			rs_alpha[i,,j] = as * ts				
		}
	}

	ret@other_params$oxy_exposure = oxy_exposure
	ret@other_params$temp_exposure = temp_exposure
	ret@other_params$rs_search = rs_search
	ret@other_params$rs_maxin = rs_maxin
	ret@other_params$rs_erepro = rs_erepro
	ret@other_params$rs_alpha = rs_alpha
	
	return(ret)
	
}

# Rate scaling function for constant physiological scaling
#	params: a MizerParams object which includes oxygen sensitivity values
#	returns params object with added mean oxygen exposure and rate scaling arrays to your params object
rs_cphys = function(params, t_max) {
	
	# Create times vector
	times = seq(0,t_max)
	
	# Create return object
	ret = params
	
	# Create empty array
	rs = array(NA, dim = c(length(times), nrow(params@species_params), length(params@w)))
	dimnames(rs) = list(time = times, sp = dimnames(params@other_params$P_crit)$sp, w = dimnames(params@other_params$P_crit)$w)
	
	# Calculate oxygen exposure
	oxy_exposure = sweep(params@other_params$benthic_occupancy, 1, params@other_params$benthic_oxygen, '*') + 
				sweep(1 - params@other_params$benthic_occupancy, 1, params@other_params$pelagic_oxygen, '*')
	
	# Calculate temperature exposure
	temp_exposure = sweep(params@other_params$benthic_occupancy, 1, params@other_params$benthic_temp, '*') + 
				sweep(1 - params@other_params$benthic_occupancy, 1, params@other_params$pelagic_temp, '*')
	
	# Calculate rate scaling
	for(i in 1:dim(rs)[1]) {
		for(j in 1:dim(rs)[3]) {
			os = scaling(U = params@species_params[,"U_crit"], k = params@other_params$P_crit[,j], a = params@species_params[,"a_crit"], oxy = oxy_exposure[i,,j])
			ts = tscale(m = params@w[j], 
					Temp = temp_exposure[i,,j], 
					Tref = params@other_params$Tref, 
					cA = params@species_params[,"cA"], 
					cD = params@species_params[,"cD"], 
					EA = params@species_params[,"EA"], 
					ED = params@species_params[,"ED"], 
					k = params@species_params[,"k"], 
					TD = params@species_params[,"TD"])
			rs[i,,j] = os * ts				
		}
	}

	ret@other_params$oxy_exposure = oxy_exposure
	ret@other_params$temp_exposure = temp_exposure
	ret@other_params$rate_scale = rs
	
	return(ret)
	
}


########## V. Metabolic scaling ##########

# Rate scaling function
metab_scale = function(U, pcrit, a, m, oxy, Temp, Tref, cA, cD, EA, ED, k, TD) {

	# Temperature scaling
	size_act = m^(cA * (Temp - Tref))
	act = exp(-(EA/k) * ((1/Temp) - (1/Tref)))
	size_deact = m^(-cD * (Temp - TD))
	deact = (1 + exp(-(ED/k) * ((1/Temp) - (1/TD))))^(-1)
	refscale = (m^(cD * (Tref - TD))) * (1 + exp(-(ED/k) * ((1/Tref) - (1/TD))))
	ts = size_act * act * size_deact * deact * refscale
	
	# # Oxygen scaling
	# os = 2 * a * ((exp((-U/2) * ((pcrit/oxy) - 1))) * ((1 + exp(-U * ((pcrit/oxy) - 1)))^(-1)))
	
	# Oxygen scaling
	os = 1 + exp(-U*(oxy-a*pcrit))
	
	return(ts*os)
}

mscale = function(params, t_max) {
	
	# Create times vector
	times = seq(0,t_max)
	
	# Create empty array
	rs_metab = array(NA, dim = c(length(times), nrow(params@species_params), length(params@w)))
	dimnames(rs_metab) = list(time = times, sp = dimnames(params@other_params$P_crit)$sp, w = dimnames(params@other_params$P_crit)$w)

	# Calculate metabolic scaling
	for(i in 1:dim(rs_metab)[1]) {
		for(j in 1:dim(rs_metab)[3]) {
			rs_metab[i,,j] = metab_scale(U = params@species_params[,"U_met"], 
						pcrit = params@other_params$P_crit[,j], 
						a = params@species_params[,"a_met"],
						m = params@w[j],
						oxy = params@other_params$oxy_exposure[i,,j], 
						Temp = params@other_params$temp_exposure[i,,j], 
						Tref = params@other_params$Tref, 
						cA = params@species_params[,"cA"], 
						cD = params@species_params[,"cD"], 
						EA = params@species_params[,"EA"], 
						ED = params@species_params[,"ED"], 
						k = params@species_params[,"k"], 
						TD = params@species_params[,"TD"])
		}
	}
	
	ret = params
	ret@other_params$metab_scale = rs_metab
	return(ret)
}
