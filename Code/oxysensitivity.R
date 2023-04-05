# Add oxygen sensitivity for each species to the params object

# Created: March 3, 2021
# Last modified: September 5, 2021 by EPD

# Contents (ctrl-f):
#	I. Calculate P_crit
#	II. Oxygen sensitivity


########## I. Calculate P_crit ##########

# Convert kPa to concentration with Henry's law and van't Hoff equation
#	kpa: oxygen pressure
#	temp: temperature in kelvin
#	returns temperature corrected aqueous concentration of oxygen in M
hvh = function(kpa, temp) {
	H = (1/770) * exp(1700 * ((1/temp) - (1/298.15)))
	M = kpa * 0.00986923
	return(M * H)
}

# Convert mizer metabolism (g/year) to RMR (mg O2/kg/h)
#	params: Mizer parameters object
#	Notes:
#		numerator: 1000 converts g O2 to mg O2; 1.429 converts L O2 to g O2; 1.3 converts g food to L O2 required
#		denominator: (1.076 converts kg body mass to L;) 1000 converts g to kg; 24 converts days to hours; 365.25 converts year to days
get_rmr = function(params) {
	rmr = 1000*1.429*1.3*params@metab / (365.25*24*matrix(rep(params@w,each=4),nrow=4)/1000)
	return(rmr)
}

# Calculate oxygen sensitivity from salinity, temperature, body size, and RMR
#	params: mizer parameter object
#	mod: results of stepwise linear regression
#	psu: salinity in practical salinity units
#	temp: temperature in degrees C
#	w: weight in kg
#	returns pcrit in kPa
pcrit = function(params, mod, psu, temp, w) {
	rmr = get_rmr(params)
	P_crit_kpa = predict(mod, newdata = list(Salinity_psu_edited = rep(psu, prod(dim(w))), Temperature_OK = rep(temp, prod(dim(w))), RMR_mgO2_kg_hr_edited = c(rmr), Mean_Mass_kg_edited = c(w)))
	P_crit_kpa = matrix(P_crit_kpa, nrow = nrow(w), ncol = ncol(w))
	P_crit = 1000 * 0.7 * 31.998 * hvh(P_crit_kpa, 292.15)
	dimnames(P_crit) = dimnames(params@maturity)
	return(P_crit)
}

# Add P_crit array to the MizerParams object
#	mod: a model object describing pcrit as a function of psu, temp, and weight
#	params: a MizerParams object created with bp_setup()
#	psu: salinity in practical salinity units
#	temp temperature in degrees C
#	returns modified params object
set_pcrit = function(mod, params, psu, temp) {
	
	# Create return object
	ret = params
	
	# Get weight bins for each species
	w = matrix(rep(params@w,each=4),nrow=4)/1000
	
	# Get pcrit for each species at each body size
	P_crit = pcrit(params = params, mod = mod, psu = psu, temp = temp, w = w)
	
	# Store in params
	ret@other_params$P_crit = P_crit
	
	return(ret)
}


########## II. Oxygen sensitivity ##########

# Add oxygen sensitivity of all species to the MizerParams object
#	params: a MizerParams object created with bp_setup()
#	U_hab: tolerance of benthic occupancy to hypoxia for cod, flounder, sprat, and herring, in that order
#	a_hab: adjusting habitat P_crit for cod, flounder, sprat, and herring, in that order
#	U_phys: tolerance of physiological processes (search, max intake, egg survival, assimilation IN THAT ORDER) to hypoxia for cod, flounder, sprat, and herring, in that order
#	a_phys: adjusting rate scaling (search, max intake, egg survival, assimilation IN THAT ORDER) P_crit for cod, flounder, sprat, and herring, in that order
#	z_mort: mortality at zero oxygen for cod, flounder, sprat, and herring, in that order
#	b_mort: hypoxia tolerance; higher => greater tolerance i.e. lower mortality at low oxygen
#	returns modified params object
oxy_sensitivity = function(params, U_hab, a_hab, U_search, a_search, U_maxin, a_maxin, U_erepro, a_erepro, U_alpha, a_alpha, U_met, a_met, z_mort, b_mort) {
	
	# Create return object
	ret = params
	
	# Assign habitat sensitivity to the species_params data frame
	ret@species_params$U_hab = U_hab
	ret@species_params$a_hab = a_hab
	
	# Assign search scaling sensitivity to the species_params data frame
	ret@species_params$U_search = U_search
	ret@species_params$a_search = a_search

	# Assign search scaling sensitivity to the species_params data frame
	ret@species_params$U_maxin = U_maxin
	ret@species_params$a_maxin = a_maxin

	# Assign search scaling sensitivity to the species_params data frame
	ret@species_params$U_erepro = U_erepro
	ret@species_params$a_erepro = a_erepro

	# Assign search scaling sensitivity to the species_params data frame
	ret@species_params$U_alpha = U_alpha
	ret@species_params$a_alpha = a_alpha

	# Assign metabolic scaling sensitivity to the species_params data frame
	ret@species_params$U_met = U_met
	ret@species_params$a_met = a_met

	# Assign mortality due to hypoxia exposure to the species_params data frame
	ret@species_params$z_mort = z_mort
	ret@species_params$b_mort = b_mort
	
	return(ret)
	
}