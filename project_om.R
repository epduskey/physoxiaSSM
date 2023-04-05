# Hypoxia physiological analysis for OM

# Created: September 3, 2022
# Last modified: February 19, 2023 by EPD

# Set working directory
setwd("/Users/epdus/OneDrive/Breathless/Code/Packets/physoxiaMizer")

# Load mizer
# devtools::install_github("sizespectrum/mizer") 
library(mizer)
library(rje)

# Load optimization
library(optimParallel)

# Contents (ctrl-f):
#	I. Source scripts
#	II. Create MizerParams object
#	III. Add oxygen sensitivity
#	IV. Create oxygen scenario
#	V. Create temperature scenario
#	VI. Add and check scaling
#	VII. Set rate functions
#	VIII. Run model
#	IX. Fit model
#	X. Analysis


########## I. Source scripts ##########

# Source scripts
#	spmat: creates species_params data frame
#	bpsetup: organizes parameters for benthic and pelagic habitats
#	metabscale: scaling maintenance costs i.e. metabolic rate with oxygen and temperature
#	oxysensitivity: organized parameters for sensitivity to hypoxia
#	otscale: scaling occupancy and physiological rates with oxygen and temperature
#	rate_encounter: scaled encounter rate
#	rate_erag: scaled availability of energy for reproduction and growth
#	rate_feeding: scaled feeding level
#	rate_mortality: mortality including direct mortality due to hypoxia exposure
#	rate_predation: scaled predation rate
#	rate_RDI: scaled reproduction
#	rate_resourcemort: scaled resource mortality rate
#	resource_benthic: benthic resource function
#	resource_pelagic: pelagic resource function
#	plot_funcs: assorted functions to process and plot mizer sim output
#	plot_spectra: plot spectra with both benthic and pelagic resources
#	plot_bh: plot Beverton-Holt of params object
#	plot_growth: plot growth at a given time step
#	plot_yield: plot yield at a given time step
#	plot_diet: custom plot diet functions
#	plot_fmsy: plot terminal yield for each species at a given F
#	cal_occupancy: calibrate occupancy sensitivity parameters
#	cal_biomass: calibrate for biomass on average and dynamically
#	cal_growth: calibrate for growth on average and dynamically
#	cal_yield: calibrate for yield on average and dynamically
#	mizer_replace: altered get and set rate functions to operate with custom functions
#	fsetup: get annual effort arrays
source("Code/spmat.R")
source("Code/bpsetup.R")
source("Code/oxysensitivity.R")
source("Code/otmscale.R")
source("Code/rate_encounter.R")
source("Code/rate_erag.R")
source("Code/rate_feeding.R")
source("Code/rate_mortality.R")
source("Code/rate_predation.R")
source("Code/rate_predmort.R")
source("Code/rate_RDI.R")
source("Code/rate_resourcemort.R")
source("Code/resource_benthic.R")
source("Code/resource_pelagic.R")
source("Code/plot_funcs.R")
source("Code/plot_spectra.R")
source("Code/plot_bh.R")
source("Code/plot_growth.R")
source("Code/plot_yield.R")
source("Code/plot_diet.R")
source("Code/plot_fmsy.R")
source("Code/cal_occupancy.R")
source("Code/cal_biomass.R")
source("Code/cal_growth.R")
source("Code/cal_yield.R")
source("Code/mizer_replace.R")
source("Code/fsetup.R")


########## II. Create MizerParams object ##########

# Create species_params object: see spmat.R
species_mat = sp_mat(rmax = c(1.9e+08,7.0e+08,2.1e+11,1.9e+10))

# Create MizerParams object: see bpsetup.R
params = bp_setup(species_params = species_mat$species_params,
					kappa_benthic = 10^11.0,
					kappa_pelagic = 10^11.6,
					lambda_benthic = 2.05,
					lambda_pelagic = 2.05,
					r_pp_benthic = 10,
					r_pp_pelagic = 10,
					n_benthic = 0.89,
					n_pelagic = 0.89,
					w_pp_cutoff_benthic = 100,
					w_pp_cutoff_pelagic = 0.01,
					U_oxy_benthic = 10^5,
					U_oxy_pelagic = 10^5,
					k_oxy_benthic = 0,
					k_oxy_pelagic = 0,
					p = 0.89,
					no_w = 100,
					interaction = species_mat$interaction)


########## III. Add oxygen sensitivity ##########

# Load pcrit database model
load("Parameters/Pcrit/pcrit_modstep.rda")

# Calculate P_crit from salinity, temperature, and rmr: see oxysensitivity.R
psu = 7
temp = 19 + 273.15
params = set_pcrit(mod.step, params, psu, temp)

# Add oxygen sensitivity parameters for each species: see oxysensitivity.R
#	subscript hab is occupancy sensitivity
#	subscript crit is phsyiological sensitivity
#	subscript met is metabolic sensitivity
#	subscript mort is mortality sensitivity
params = oxy_sensitivity(params, 
				U_hab = c(0.8, 10^5, -10^5, -10^5),
				a_hab = c(0.1, 0, 0, 0),
				U_search = c(10^5, 10^5, 10^5, 10^5),
				a_search = c(0, 0, 0, 0),
				U_maxin = c(10^5, 10^5, 10^5, 10^5),
				a_maxin = c(0, 0, 0, 0),
				U_erepro = c(10^5, 10^5, 10^5, 10^5),
				a_erepro = c(0, 0, 0, 0),
				U_alpha = c(10^5, 10^5, 10^5, 10^5),
				a_alpha = c(0, 0, 0, 0),
				U_met = c(10^5, 10^5, 10^5, 10^5),
				a_met = c(-10^5, -10^5, -10^5, -10^5),
				z_mort = c(0.2, 0.1, 0, 0),
				b_mort = c(1.0, 2.0, 0, 0))


########## IV. Create oxygen scenario ##########

# Load oxygen data
load("Data/oxy_model.rda")

# Create an oxygen vector
t_max = 100
times = 0:t_max
benthic_oxygen = vector(mode = "numeric", length = length(times))
pelagic_oxygen = vector(mode = "numeric", length = length(times))

# Constant oxygen scenario (kPa)
benthic_oxygen[1:length(times)] = mean(preds_all$df[preds_all$df$Year %in% seq(1991,2000), ]$Oxygen_benthic)
pelagic_oxygen[1:length(times)] = mean(preds_all$df[preds_all$df$Year %in% seq(1991,2000), ]$Oxygen_pelagic)

# Store these in other_params
params@other_params$benthic_oxygen = benthic_oxygen
params@other_params$pelagic_oxygen = pelagic_oxygen


########## V. Create temperature scenario ##########

# Choose reference temperature
Tref = 19
params@other_params$Tref = Tref

# Create a temperature vector
benthic_temp = vector(mode = "numeric", length = length(times))
pelagic_temp = vector(mode = "numeric", length = length(times))

# Constant temperature scenario
benthic_temp[1:length(times)] = 19
pelagic_temp[1:length(times)] = 19

# Store these in other_params
params@other_params$benthic_temp = benthic_temp
params@other_params$pelagic_temp = pelagic_temp


########## VI. Add and check scaling ##########4

# Scale rates by oxygen and temperature: see otmscale.R
params = occupancy(params, t_max)
params = rate_scale(params, t_max)
params = mscale(params, t_max)

# Make sure physiological scaling is set to one
params@other_params$rs_search[,,] = 1
params@other_params$rs_maxin[,,] = 1
params@other_params$rs_erepro[,,] = 1
params@other_params$rs_alpha[,,] = 1
params@other_params$metab_scale[,,] = 1


########## VII. Set rate functions ##########

# Create new params object
params_hypoxia = params

# Set benthic resource function
params_hypoxia@resource_dynamics = "resource_benthic"

# Set rate functions
params_hypoxia = setRateFunction(params_hypoxia, "Encounter", "physoxiaEncounter")
params_hypoxia = setRateFunction(params_hypoxia, "FeedingLevel", "physoxiaFeedingLevel")
params_hypoxia = setRateFunction(params_hypoxia, "PredRate", "physoxiaPredRate")
params_hypoxia = setRateFunction(params_hypoxia, "PredMort", "hypoxiaPredMort")
params_hypoxia = setRateFunction(params_hypoxia, "Mort", "hypoxiaMort")
params_hypoxia = setRateFunction(params_hypoxia, "ResourceMort", "benthicResourceMort")
params_hypoxia = setRateFunction(params_hypoxia, "EReproAndGrowth", "physoxiaEReproAndGrowth")
params_hypoxia = setRateFunction(params_hypoxia, "RDI", "physoxiaRDI")

# Set a dynamic pelagic resource component
params_hypoxia = setComponent(params = params_hypoxia,
					component = "n_pp_pelagic",
					initial_value = params@resource_params$pelagic$cc_pp,
					dynamics_fun = "resource_pelagic")


########## VIII. Run model ##########

# Fishing effort
effort_hypoxia = c(1,1,1,1)

# # Test convergence
# sim_converge = projectToSteady(params_hypoxia, effort = effort_hypoxia)

# Run model
sim_hypoxia = project(params_hypoxia, t_max = t_max, effort = effort_hypoxia)
plotBiomass(sim_hypoxia) 
plotFeedingLevel(sim_hypoxia)
plotPredMort(sim_hypoxia) 
plotFMort(sim_hypoxia) 
plotSpectra(sim_hypoxia)

# Growth curves
plotGrowthCurves(sim_hypoxia, species = "cod", max_age = 15)
plotGrowthCurves(sim_hypoxia, species = "flounder", max_age = 26)
plotGrowthCurves(sim_hypoxia, species = "sprat", max_age = 16)
plotGrowthCurves(sim_hypoxia, species = "herring", max_age = 13)

# Diet
plotHypoxiaDiet(sim_hypoxia, "cod", t = t_max+1)
plotHypoxiaDiet(sim_hypoxia, "flounder", t_max+1)
plotHypoxiaDiet(sim_hypoxia, "sprat", t_max+1)
plotHypoxiaDiet(sim_hypoxia, "herring", t_max+1)

# Biomass
plotBiomassObservedVsModel(sim_hypoxia)

# Yield
plotYieldObservedVsModel(sim_hypoxia, t_max+1, seq(1991,2000))

# # Check Fmsy
# plotFmsy(params_hypoxia, effort_res = 40, t_max = t_max)


########## IX. Fit model ##########

# Set up new params tuning object
params_tune = params_hypoxia

# Set reproduction level
params_tune = setBevertonHolt(setInitialValues(params_tune, sim_hypoxia), reproduction_level = c(0.90,0.90,0.32,0.86))

# # Check Fmsy
# plotFmsy(params_tune, effort_res = 20, max_f = 1.5, t_max = 100)

# Choose tune runtime
t_tune = 100

# Use test scenarios to choose starting values for mortality, clearance rate, and max consumption
params_test = params_tune
species_params(params_test)$z0 = c(0.00, 0.21, 0.34, 0.26)
species_params(params_test)$gamma = c(4.8e-11, 3.2e-11, 2.5e-11, 3.0e-11)
sim_test = project(params_test, t_max = t_tune, effort = effort_hypoxia)
plotBiomass(sim_test)
plotFeedingLevel(sim_test)
plotGrowthCurves(sim_test, species = "cod", max_age = 15)
plotGrowthCurves(sim_test, species = "flounder", max_age = 26)
plotGrowthCurves(sim_test, species = "sprat", max_age = 16)
plotGrowthCurves(sim_test, species = "herring", max_age = 13)
plotBiomassObservedVsModel(sim_test)
plotYieldObservedVsModel(sim_test, 101, seq(1991,2000))

# Fit occupancy
species_params(params_tune)$z0 = c(0.00, 0.21, 0.34, 0.26)
species_params(params_tune)$gamma = c(4.8e-11, 3.2e-11, 2.5e-11, 3.0e-11)
par_occ = c(0.8, 0.1)
yr_occ = seq(1991,2000)
t_occ = 100

# Set up cores and run
noCores = detectCores() - 1
cl = makeCluster(noCores, setup_timeout = 0.5)
setDefaultCluster(cl = cl)
clusterExport(cl, as.list(ls()))
clusterEvalQ(cl, {
	library(mizer)
	library(optimParallel)
})
occ.optim = optimParallel(par = par_occ, fn = bocc_optim, params = params_tune, dat = cod_benthic, yr = yr_occ, t = t_occ, effort = effort_hypoxia, lower = c(0,0), method = "L-BFGS-B", parallel=list(cl=cl,loginfo=T))
stopCluster(cl)

# Set up params object with fitted habitat values
params_tune = oxy_sensitivity(params_tune, 
				U_hab = c(occ.optim$par[1], 10^5, -10^5, -10^5),
				a_hab = c(occ.optim$par[2], 0, 0, 0),
				U_search = c(10^5, 10^5, 10^5, 10^5),
				a_search = c(0, 0, 0, 0),
				U_maxin = c(10^5, 10^5, 10^5, 10^5),
				a_maxin = c(0, 0, 0, 0),
				U_erepro = c(10^5, 10^5, 10^5, 10^5),
				a_erepro = c(0, 0, 0, 0),
				U_alpha = c(10^5, 10^5, 10^5, 10^5),
				a_alpha = c(0, 0, 0, 0),
				U_met = c(10^5, 10^5, 10^5, 10^5),
				a_met = c(-10^5, -10^5, -10^5, -10^5),
				z_mort = c(0.2, 0.1, 0, 0),
				b_mort = c(1.0, 2.0, 0, 0))

# Scale rates by oxygen and temperature: see otmscale.R
params_time = occupancy(params_tune, t_tune)

# Run to steady state
sim_tune = project(params_tune, t_max = t_tune, effort = effort_hypoxia)
plotBiomass(sim_tune)

# Check growth
plotGrowthCurves(sim_test, species = "cod", max_age = 15)
plotGrowthCurves(sim_test, species = "flounder", max_age = 26)
plotGrowthCurves(sim_test, species = "sprat", max_age = 16)
plotGrowthCurves(sim_test, species = "herring", max_age = 13)

# Check biomass
plotBiomassObservedVsModel(sim_test)

# Check yield
plotYieldObservedVsModel(sim_test, 101, seq(1991,2000))

# Fit yield
par_yield = species_params(params_tune)$z0
dat_yield = list(cod_catch, flounder_catch, sprat_catch, herring_catch)
yr_yield = seq(1991,2000)
t_yield = 100

# Set up cores and run
noCores = detectCores() - 1
cl = makeCluster(noCores, setup_timeout = 0.5)
setDefaultCluster(cl = cl)
clusterExport(cl, as.list(ls()))
clusterEvalQ(cl, {
	library(mizer)
	library(optimParallel)
})
yield.optim = optimParallel(par = par_yield, fn = fit_yield, params = params_tune, dat = dat_yield, yr = yr_yield, t = t_yield, effort = effort_hypoxia, lower = c(0,0,0,0), upper = c(1,1,1,1), method = "L-BFGS-B", parallel=list(cl=cl,loginfo=T))
stopCluster(cl)

# Set new values
species_params(params_tune)$z0 = yield.optim$par

# Run to steady state
sim_tune = project(params_tune, t_max = t_tune, effort = effort_hypoxia)
plotBiomass(sim_tune)

# Check growth and biomass
plotGrowthCurves(sim_tune, species = "cod", max_age = 15)
plotGrowthCurves(sim_tune, species = "flounder", max_age = 26)
plotGrowthCurves(sim_tune, species = "sprat", max_age = 16)
plotGrowthCurves(sim_tune, species = "herring", max_age = 13)
plotBiomassObservedVsModel(sim_tune)
plotYieldObservedVsModel(sim_tune, t_yield, yr_yield)

# Check diet
plotHypoxiaDiet(sim_tune, "cod", t_tune+1)
plotHypoxiaDiet(sim_tune, "flounder", t_tune+1)
plotHypoxiaDiet(sim_tune, "sprat", t_tune+1)
plotHypoxiaDiet(sim_tune, "herring", t_tune+1)

# Use test scenarios to find good starting values for rmax
params_test = params_tune
species_params(params_test)$R_max = c(7.3e+07,3.7e+08,6.0e+11,2.4e+10)
sim_test = project(params_test, t_max = t_tune, effort = effort_hypoxia)
plotBiomassObservedVsModel(sim_test)
plotBiomass(sim_test)

# Fit biomass and growth
par_bg = c(log10(c(7.3e+07,3.7e+08,6.0e+11,2.4e+10)), log10(species_params(params_tune)$gamma))
dat_bg = list(cod_ssb, flounder_ssb, sprat_ssb, herring_ssb)
yr_bg = seq(1991,2000)
t_bg = 100

# Set up cores and run
noCores = detectCores() - 1
cl = makeCluster(noCores, setup_timeout = 0.5)
setDefaultCluster(cl = cl)
clusterExport(cl, as.list(ls()))
clusterEvalQ(cl, {
	library(mizer)
	library(optimParallel)
})
bg.optim = optimParallel(par = par_bg, params = params_tune, fit_bg, dat = dat_bg, yr = yr_bg, t = t_bg, effort = effort_hypoxia, lower = c(5,5,5,5,-15,-15,-15,-15), upper = c(15,15,15,15,-5,-5,-5,-5), method = "L-BFGS-B", control = list(factr=1e10), parallel = list(cl=cl,loginfo=T))
stopCluster(cl)

# Set new values
species_params(params_tune)$R_max = 10^bg.optim$par[1:4]
species_params(params_tune)$gamma = 10^bg.optim$par[5:8]

# Run to steady state
sim_tune = project(params_tune, t_max = t_tune, effort = effort_hypoxia)
plotBiomass(sim_tune)

# Check biomass
plotBiomassObservedVsModel(sim_tune)

# Check growth
plotGrowthCurves(sim_tune, species = "cod", max_age = 15)
plotGrowthCurves(sim_tune, species = "flounder", max_age = 26)
plotGrowthCurves(sim_tune, species = "sprat", max_age = 16)
plotGrowthCurves(sim_tune, species = "herring", max_age = 13)

# Check yield
plotYieldObservedVsModel(sim_tune, t_tune+1, seq(1991,2000))

# Set up new objects for time params
cal_ts = seq(1991,2000)
params_time = params_tune

# Get effort array
effort_time = farray(params_tune, fdat, seq(1991,2000), c(2,2,2,2), t_tune)

# Create an oxygen vector
times = 0:t_tune
benthic_oxygen = vector(mode = "numeric", length = length(times))
pelagic_oxygen = vector(mode = "numeric", length = length(times))

# Constant oxygen scenario (kPa)
benthic_oxygen[1:length(times)] = c(rep(mean(preds_all$df[preds_all$df$Year %in% cal_ts, ]$Oxygen_benthic),t_tune-length(cal_ts)+1), preds_all$df[preds_all$df$Year %in% cal_ts, ]$Oxygen_benthic)
pelagic_oxygen[1:length(times)] = c(rep(mean(preds_all$df[preds_all$df$Year %in% cal_ts, ]$Oxygen_pelagic),t_tune-length(cal_ts)+1), preds_all$df[preds_all$df$Year %in% cal_ts, ]$Oxygen_pelagic)

# Store these in other_params
params_time@other_params$benthic_oxygen = benthic_oxygen
params_time@other_params$pelagic_oxygen = pelagic_oxygen

# Choose reference temperature
Tref = 19
params_time@other_params$Tref = Tref

# Create a temperature vector
benthic_temp = vector(mode = "numeric", length = length(times))
pelagic_temp = vector(mode = "numeric", length = length(times))

# Constant temperature scenario
benthic_temp[1:length(times)] = 19
pelagic_temp[1:length(times)] = 19

# Store these in other_params
params_time@other_params$benthic_temp = benthic_temp
params_time@other_params$pelagic_temp = pelagic_temp

# Scale rates by oxygen and temperature: see otmscale.R
params_time = occupancy(params_time, t_tune)

# Run time series model
sim_time = project(params_time, t_max = t_tune, effort = effort_time)
plotBiomass(sim_time)

# Fit mortality parameters for cod
par_cod_mort = c(species_params(params_time)$z_mort[1], species_params(params_time)$b_mort[1])
yr_mort = seq(1991,2000)
t_mort = 100

# Set up cores and run
noCores = detectCores() - 1
cl = makeCluster(noCores, setup_timeout = 0.5)
setDefaultCluster(cl = cl)
clusterExport(cl, as.list(ls()))
clusterEvalQ(cl, {
	library(mizer)
	library(optimParallel)
	library(rje)
})
cod.mort.optim = optimParallel(par = par_cod_mort, params = params_time, dat = dat_yield, species = "cod", time_yield, yr = yr_mort, t = t_mort, effort = effort_time, lower = c(0,0), method = "L-BFGS-B", parallel = list(cl=cl,loginfo=T))
stopCluster(cl)

# Reset parameters
species_params(params_time)$z_mort[1] = cod.mort.optim$par[1]
species_params(params_time)$b_mort[1] = cod.mort.optim$par[2]

# Run time series model
sim_time = project(params_time, t_max = t_tune, effort = effort_time)
plotBiomass(sim_time)

# Fit mortality parameters for flounder
par_flounder_mort = c(species_params(params_time)$z_mort[2], species_params(params_time)$b_mort[2])
yr_mort = seq(1991,2000)
t_mort = 100

# Set up cores and run
noCores = detectCores() - 1
cl = makeCluster(noCores, setup_timeout = 0.5)
setDefaultCluster(cl = cl)
clusterExport(cl, as.list(ls()))
clusterEvalQ(cl, {
	library(mizer)
	library(optimParallel)
	library(rje)
})
flounder.mort.optim = optimParallel(par = par_flounder_mort, params = params_time, dat = dat_yield, species = "flounder", time_yield, yr = yr_mort, t = t_mort, effort = effort_time, lower = c(0,0), method = "L-BFGS-B", parallel = list(cl=cl,loginfo=T))
stopCluster(cl)

# Reset parameters
species_params(params_time)$z_mort[2] = flounder.mort.optim$par[1]
species_params(params_time)$b_mort[2] = flounder.mort.optim$par[2]

# Run time series model
sim_time = project(params_time, t_max = t_tune, effort = effort_time)
plotBiomass(sim_time)

# Plot growth time series
par(mfrow = c(4,10), mar = c(0,0,3,0), oma = c(5,5,3,1))
for(i in 1:length(cal_ts)) {
	tempmod = myGrowthCurves(sim_time, (t_tune-length(cal_ts))+i, "cod", max_age = 15)
	tempdat = cod_ivb[cod_ivb$Year == cal_ts[i], ]
	tempdat$Prop = tempdat$IndWgt/species_params(params_time)["cod",]$w_inf
	tempdat$Age = tempdat$Age + 0.5
	plot(Prop ~ Age, tempdat, axes = F, xlim = c(-2,17)); axis(1,at=c(0,5,10,15)); box()
	if(i == 1) {axis(2,at=c(0,1),labels=c(0,round(species_params(params_time)["cod",]$w_inf)))}
	tempmod = tempmod/species_params(params_time)["cod",]$w_inf
	lines(dimnames(tempmod)$Age, tempmod, lwd = 2, lty = 1, col = "red")
	tempobs = cod_ivb_cal[i,]$a*(cod_ivb_cal[i,]$L_inf*(1-exp(-cod_ivb_cal[i,]$k*(seq(0,15,length.out=50)))))^cod_ivb_cal[i,]$b
	tempobs = tempobs/species_params(params_time)["cod",]$w_inf
	lines(seq(0,15,length.out=50), tempobs, lwd = 2, col = "red", lty = 2)
}
for(i in 1:length(cal_ts)) {
	tempmod = myGrowthCurves(sim_time, (t_tune-length(cal_ts))+i, "flounder", max_age = 26)
	tempdat = flounder_ivb
	tempdat$Prop = tempdat$IndWgt/species_params(params_time)["flounder",]$w_inf
	tempdat$Age = tempdat$Age + 0.5
	plot(Prop ~ Age, tempdat, axes = F, xlim = c(-2,28), ylim = c(0,1)); axis(1,at=c(0,10,20)); box()
	if(i == 1) {axis(2,at=c(0,1),labels=c(0,round(species_params(params_time)["flounder",]$w_inf)))}
	tempmod = tempmod/species_params(params_time)["flounder",]$w_inf
	lines(dimnames(tempmod)$Age, tempmod, lwd = 2, lty = 1, col = "red")
	tempobs = flounder_lh["a","stan"]*(flounder_lh["L_inf","stan"]*(1-exp(-flounder_lh["k","stan"]*(seq(0,26,length.out=50)))))^flounder_lh["b","stan"]
	tempobs = tempobs/species_params(params_time)["flounder",]$w_inf
	lines(seq(0,26,length.out=50), tempobs, lwd = 2, col = "red", lty = 2)
}
for(i in 1:length(cal_ts)) {
	tempmod = myGrowthCurves(sim_time, (t_tune-length(cal_ts))+i, "sprat", max_age = 16)
	tempmod = tempmod/species_params(params_time)["sprat",]$w_inf
	tempobs = sprat_ivb_cal[i,]$a*(sprat_ivb_cal[i,]$L_inf*(1-exp(-sprat_ivb_cal[i,]$k*(seq(0,16,length.out=50)))))^sprat_ivb_cal[i,]$b
	tempobs = tempobs/species_params(params_time)["sprat",]$w_inf
	tempdat = sprat_ivb[sprat_ivb$Year == cal_ts[i], ]
	tempdat$Age = tempdat$Age + 0.5
	if(nrow(tempdat)>0) {
		tempdat$Prop = tempdat$IndWgt/species_params(params_time)["sprat",]$w_inf
		plot(Prop ~ Age, tempdat, axes = F, xlim = c(-2,18)); axis(1,at=c(0,5,10,15)); box()
		lines(dimnames(tempmod)$Age, tempmod, lwd = 2, lty = 1, col = "red")
		lines(seq(0,16,length.out=50), tempobs, lwd = 2, col = "red", lty = 2)
	} else{
		plot(seq(0,16,length.out=50), tempobs, type = 'l', lwd = 2, lty = 1, col = "red", axes = F, xlim = c(-2,18), ylim = c(0,1)); axis(1,at=c(0,5,10,15)); box()
		lines(dimnames(tempmod)$Age, tempmod, lwd = 2, col = "red", lty = 2)
	}
	if(i == 1) {axis(2,at=c(0,1),labels=c(0,round(species_params(params_time)["sprat",]$w_inf)))}
}
for(i in 1:length(cal_ts)) {
	tempmod = myGrowthCurves(sim_time, (t_tune-length(cal_ts))+i, "herring", max_age = 13)
	tempmod = tempmod/species_params(params_time)["herring",]$w_inf
	tempobs = herring_ivb_cal[i,]$a*(herring_ivb_cal[i,]$L_inf*(1-exp(-herring_ivb_cal[i,]$k*(seq(0,13,length.out=50)))))^herring_ivb_cal[i,]$b
	tempobs = tempobs/species_params(params_time)["herring",]$w_inf
	tempdat = herring_ivb[herring_ivb$Year == cal_ts[i], ]
	tempdat$Age = tempdat$Age + 0.5
	if(nrow(tempdat)>0) {
		tempdat$Prop = tempdat$IndWgt/species_params(params_time)["herring",]$w_inf
		plot(Prop ~ Age, tempdat, axes = F, xlim = c(-2,15)); axis(1,at=c(0,5,10)); box()
		lines(dimnames(tempmod)$Age, tempmod, lwd = 2, lty = 1, col = "red")
		lines(seq(0,13,length.out=50), tempobs, lwd = 2, col = "red", lty = 2)
	} else{
		plot(seq(0,13,length.out=50), tempobs, type = 'l', lwd = 2, lty = 1, col = "red", axes = F, xlim = c(-2,18), ylim = c(0,1)); axis(1,at=c(0,5,10,15)); box()
		lines(dimnames(tempmod)$Age, tempmod, lwd = 2, col = "red", lty = 2)
	}
	if(i == 1) {axis(2,at=c(0,1),labels=c(0,round(species_params(params_time)["herring",]$w_inf)))}
}
mtext("Age", side = 1, line = 3, outer = T, cex = 1.5)
mtext("Weight (g)", side = 2, line = 2, outer = T, cex = 1.5)
par(mfrow = c(1,1))

# Plot biomass time series
par(mfrow = c(2,2))
plot(cal_ts, rowSums((sweep(sim_time@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params@w*params@dw, "*"))[,params@w >= species_params(params_time)$w_mat[1]]), ylim = c(3e10,3e11), log = 'y', type = 'l', lwd = 2, xlab = "Year", ylab = "SSB (g)", main = "Cod")
points(cod_ssb$Year, cod_ssb$SSB, pch = 16, col = "red")
plot(cal_ts, rowSums((sweep(sim_time@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params@w*params@dw, "*"))[,params@w >= species_params(params_time)$w_mat[2]]), ylim = c(8e08,8e10), log = 'y', type = 'l', lwd = 2, xlab = "Year", ylab = "SSB (g)", main = "Flounder")
points(flounder_ssb$Year, flounder_ssb$SSB, pch = 16, col = "red")
plot(cal_ts, rowSums((sweep(sim_time@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params@w*params@dw, "*"))[,params@w >= species_params(params_time)$w_mat[3]]), ylim = c(5e11,3e12), log = 'y', type = 'l', lwd = 2, xlab = "Year", ylab = "SSB (g)", main = "Sprat")
points(sprat_ssb$Year, sprat_ssb$SSB, pch = 16, col = "red")
plot(cal_ts, rowSums((sweep(sim_time@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params@w*params@dw, "*"))[,params@w >= species_params(params_time)$w_mat[4]]), ylim = c(2e11,1e12), log = 'y', type = 'l', lwd = 2, xlab = "Year", ylab = "SSB (g)", main = "Herring")
points(herring_ssb$Year, herring_ssb$SSB, pch = 16, col = "red")
par(mfrow = c(1,1))

# Plot yield time series
yield_time = getYield(sim_time)
par(mfrow = c(2,2))
plot(cal_ts, yield_time[(t_tune-length(cal_ts)+1):t_tune, "cod"], log = 'y', ylim = c(3e10,3e11), type = 'l', lwd = 2, xlab = "Year", ylab = "Yield (g)", main = "Cod")
points(cal_ts, cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000, pch = 16, col = "red")
plot(cal_ts, yield_time[(t_tune-length(cal_ts)+1):t_tune, "flounder"], log = 'y', ylim = c(3e9,3e11), type = 'l', lwd = 2, xlab = "Year", ylab = "Yield (g)", main = "Flounder")
points(cal_ts, flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000, pch = 16, col = "red")
plot(cal_ts, yield_time[(t_tune-length(cal_ts)+1):t_tune, "sprat"], log = 'y', ylim = c(1e11,1e12), type = 'l', lwd = 2, xlab = "Year", ylab = "Yield (g)", main = "Sprat")
points(cal_ts, sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000, pch = 16, col = "red")
plot(cal_ts, yield_time[(t_tune-length(cal_ts)+1):t_tune, "herring"], log = 'y', ylim = c(5e10,5e11), type = 'l', lwd = 2, xlab = "Year", ylab = "Yield (g)", main = "Herring")
points(cal_ts, herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000, pch = 16, col = "red")
par(mfrow = c(1,1))

# Check mean again
sim_tune = project(params_time, t_max = 90, effort = effort_time[1:90, ])
plotBiomass(sim_tune) 
plotFeedingLevel(sim_tune)
plotPredMort(sim_tune) 
plotFMort(sim_tune) 
plotSpectra(sim_tune)
plot_spectra(sim_tune, 90)

# Check growth and biomass
plotGrowthCurves(sim_tune, species = "cod", max_age = 15)
plotGrowthCurves(sim_tune, species = "flounder", max_age = 26)
plotGrowthCurves(sim_tune, species = "sprat", max_age = 16)
plotGrowthCurves(sim_tune, species = "herring", max_age = 13)
plotBiomassObservedVsModel(sim_tune)
plotYieldObservedVsModel(sim_tune, t_tune/2, yr_yield)

# Check diet
plotHypoxiaDiet(sim_tune, "cod", t_tune/2)
plotHypoxiaDiet(sim_tune, "flounder", t_tune/2)
plotHypoxiaDiet(sim_tune, "sprat", t_tune/2)
plotHypoxiaDiet(sim_tune, "herring", t_tune/2)

# Check yield versus F
plotFmsy(params_time)

# Save the fitted objects
om_fit = list(
	occ = occ.optim,
	yield = yield.optim,
	bg = bg.optim,
	codmort = cod.mort.optim,
	floundermort = flounder.mort.optim
)
save(object = om_fit, file = "Calibration/om_fit.rda")

# Save the final objects
save(object = params_tune, file = "Calibration/params_tune_om.rda")
save(object = params_time, file = "Calibration/params_time_om.rda")


########## X. Analysis ##########

# Load the final objects
load("Calibration/params_tune_om.rda")
load("Calibration/params_time_om.rda")

# Create object for analysis
params_full = params_time

# Create an oxygen vector
times = 0:t_tune
benthic_oxygen = vector(mode = "numeric", length = length(times))
pelagic_oxygen = vector(mode = "numeric", length = length(times))

# Varying oxygen scenario
benthic_oxygen[1:length(times)] = c(rep(mean(preds_all$df[preds_all$df$Year %in% seq(1991,2000), ]$Oxygen_benthic),t_tune-dim(preds_all$df[preds_all$df$Year %in% seq(1991,2019), ])[1]+1), preds_all$df[preds_all$df$Year %in% seq(1991,2019), ]$Oxygen_benthic)
pelagic_oxygen[1:length(times)] = c(rep(mean(preds_all$df[preds_all$df$Year %in% seq(1991,2000), ]$Oxygen_pelagic),t_tune-dim(preds_all$df[preds_all$df$Year %in% seq(1991,2019), ])[1]+1), preds_all$df[preds_all$df$Year %in% seq(1991,2019), ]$Oxygen_pelagic)

# Store these in other_params
params_full@other_params$benthic_oxygen = benthic_oxygen
params_full@other_params$pelagic_oxygen = pelagic_oxygen

# Scale rates by oxygen and temperature: see otmscale.R
params_full = occupancy(params_full, t_tune)

# Get effort array
effort_full = farray(params_tune, fdat, seq(1991,2019), c(2,2,2,2), t_tune)

# Run model
sim_full = project(params_full, t_max = t_tune, effort = effort_full)
plotBiomass(sim_full)

# Plot biomass time series
cal_ts = seq(1991,2019)
par(mfrow = c(2,2))
plot(cal_ts, rowSums((sweep(sim_full@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params@w*params@dw, "*"))[,params@w >= species_params(params_full)$w_mat[1]]), ylim = c(3e10,3e11), log = 'y', type = 'l', lwd = 2, xlab = "Year", ylab = "SSB (g)", main = "Cod")
points(cod_ssb$Year, cod_ssb$SSB, pch = 16, col = "red")
plot(cal_ts, rowSums((sweep(sim_full@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params@w*params@dw, "*"))[,params@w >= species_params(params_full)$w_mat[2]]), ylim = c(8e8,8e10), log = 'y', type = 'l', lwd = 2, xlab = "Year", ylab = "SSB (g)", main = "Flounder")
points(flounder_ssb$Year, flounder_ssb$SSB, pch = 16, col = "red")
plot(cal_ts, rowSums((sweep(sim_full@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params@w*params@dw, "*"))[,params@w >= species_params(params_full)$w_mat[3]]), ylim = c(5e11,3e12), log = 'y', type = 'l', lwd = 2, xlab = "Year", ylab = "SSB (g)", main = "Sprat")
points(sprat_ssb$Year, sprat_ssb$SSB, pch = 16, col = "red")
plot(cal_ts, rowSums((sweep(sim_full@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params@w*params@dw, "*"))[,params@w >= species_params(params_full)$w_mat[4]]), ylim = c(2e11,1e12), log = 'y', type = 'l', lwd = 2, xlab = "Year", ylab = "SSB (g)", main = "Herring")
points(herring_ssb$Year, herring_ssb$SSB, pch = 16, col = "red")
par(mfrow = c(1,1))

# Plot yield time series
yield_full = getYield(sim_full)
par(mfrow = c(2,2))
plot(cal_ts, yield_full[(t_tune-length(cal_ts)+1):t_tune, "cod"], log = 'y', ylim = c(2e10,2e11), type = 'l', lwd = 2, xlab = "Year", ylab = "Yield (g)", main = "Cod")
points(cal_ts, cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000, pch = 16, col = "red")
plot(cal_ts, yield_full[(t_tune-length(cal_ts)+1):t_tune, "flounder"], log = 'y', ylim = c(3e9,3e11), type = 'l', lwd = 2, xlab = "Year", ylab = "Yield (g)", main = "Flounder")
points(cal_ts, flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000, pch = 16, col = "red")
plot(cal_ts, yield_full[(t_tune-length(cal_ts)+1):t_tune, "sprat"], log = 'y', ylim = c(1e11,1e12), type = 'l', lwd = 2, xlab = "Year", ylab = "Yield (g)", main = "Sprat")
points(cal_ts, sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000, pch = 16, col = "red")
plot(cal_ts, yield_full[(t_tune-length(cal_ts)+1):t_tune, "herring"], log = 'y', ylim = c(5e10,5e11), type = 'l', lwd = 2, xlab = "Year", ylab = "Yield (g)", main = "Herring")
points(cal_ts, herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000, pch = 16, col = "red")
par(mfrow = c(1,1))

save(object = params_full, file = "Calibration/params_full_om.rda")