# Load and plot calibration results

# Created: March 8, 2022
# Last modified: March 13, 2023 by EPD

# Set working directory
setwd("/Users/epdus/OneDrive/Breathless/Code/Packets/physoxiaMizer")

# Load packages
library(mizer)
library(RColorBrewer)
library(png)
library(xtable)

# Contents (ctrl-f):
#	I. Source scripts
#	II. Load results
#	III. Prepare and run calibrations
#	IV. Prepare and run simulations
#	V. Calculate calibration error
#	VI. Calculate projection error


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


########## II. Load results ##########

# Benthos, occupancy, mortality, and physiological scaling
load("Calibration/params_full_bomp.rda")
params_bomp = params_full

# Benthos, occupancy, and mortality scaling
load("Calibration/params_full_bom.rda")
params_bom = params_full

# Benthos, occupancy, and physiological scaling
load("Calibration/params_full_bop.rda")
params_bop = params_full

# Benthos, mortality, and physiological scaling
load("Calibration/params_full_bmp.rda")
params_bmp = params_full

# Occupancy, mortality, and physiological scaling
load("Calibration/params_full_omp.rda")
params_omp = params_full

# Benthos and occupancy scaling
load("Calibration/params_full_bo.rda")
params_bo = params_full

# Benthos and mortality scaling
load("Calibration/params_full_bm.rda")
params_bm = params_full

# Benthos and physiological scaling
load("Calibration/params_full_bp.rda")
params_bp = params_full

# Occupancy and mortality scaling
load("Calibration/params_full_om.rda")
params_om = params_full

# Occupancy and physiological scaling
load("Calibration/params_full_op.rda")
params_op = params_full

# Occupancy and physiological scaling
load("Calibration/params_full_mp.rda")
params_mp = params_full

# Benthos scaling
load("Calibration/params_full_b.rda")
params_b = params_full

# Occupancy scaling
load("Calibration/params_full_o.rda")
params_o = params_full

# Mortality scaling
load("Calibration/params_full_m.rda")
params_m = params_full

# Physiological scaling
load("Calibration/params_full_p.rda")
params_p = params_full

# None scaling (with left fish)
load("Calibration/params_full_none.rda")
params_none = params_full


########## III. Prepare and run calibrations ##########

# Load oxygen data
load("Data/oxy_model.rda")

# Create an oxygen vector
t_tune = 100
times = 0:t_tune
benthic_oxygen_cal = vector(mode = "numeric", length = length(times))
pelagic_oxygen_cal = vector(mode = "numeric", length = length(times))

# Varying oxygen scenario
benthic_oxygen_cal[1:length(times)] = c(rep(mean(preds_all$df[preds_all$df$Year %in% seq(1991,2000), ]$Oxygen_benthic),t_tune-dim(preds_all$df[preds_all$df$Year %in% seq(1991,2000), ])[1]+1), preds_all$df[preds_all$df$Year %in% seq(1991,2000), ]$Oxygen_benthic)
pelagic_oxygen_cal[1:length(times)] = c(rep(mean(preds_all$df[preds_all$df$Year %in% seq(1991,2000), ]$Oxygen_pelagic),t_tune-dim(preds_all$df[preds_all$df$Year %in% seq(1991,2000), ])[1]+1), preds_all$df[preds_all$df$Year %in% seq(1991,2000), ]$Oxygen_pelagic)

# Store these in other_params
params_bomp_cal = params_bomp
params_bomp_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_bomp_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_bom_cal = params_bom
params_bom_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_bom_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_bop_cal = params_bop
params_bop_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_bop_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_bmp_cal = params_bmp
params_bmp_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_bmp_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_omp_cal = params_omp
params_omp_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_omp_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_bo_cal = params_bo
params_bo_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_bo_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_bm_cal = params_bm
params_bm_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_bm_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_bp_cal = params_bp
params_bp_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_bp_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_om_cal = params_om
params_om_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_om_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_op_cal = params_op
params_op_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_op_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_mp_cal = params_mp
params_mp_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_mp_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_b_cal = params_b
params_b_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_b_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_o_cal = params_bomp
params_o_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_o_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_m_cal = params_m
params_m_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_m_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_p_cal = params_p
params_p_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_p_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

params_none_cal = params_none
params_none_cal@other_params$benthic_oxygen = benthic_oxygen_cal
params_none_cal@other_params$pelagic_oxygen = pelagic_oxygen_cal

# Scale rates by oxygen and temperature: see otmscale.R
params_bomp_cal = occupancy(params_bomp_cal, t_tune)
params_bomp_cal = rate_scale(params_bomp_cal, t_tune)
params_bomp_cal = mscale(params_bomp_cal, t_tune)

params_bom_cal = occupancy(params_bom_cal, t_tune)

params_bop_cal = occupancy(params_bop_cal, t_tune)
params_bop_cal = rate_scale(params_bop_cal, t_tune)
params_bop_cal = mscale(params_bop_cal, t_tune)

params_bmp_cal = rate_scale(params_bmp_cal, t_tune)
params_bmp_cal = mscale(params_bmp_cal, t_tune)

params_omp_cal = occupancy(params_omp_cal, t_tune)
params_omp_cal = rate_scale(params_omp_cal, t_tune)
params_omp_cal = mscale(params_omp_cal, t_tune)

params_bo_cal = occupancy(params_bo_cal, t_tune)

params_bp_cal = rate_scale(params_bp_cal, t_tune)
params_bp_cal = mscale(params_bp_cal, t_tune)

params_om_cal = occupancy(params_om_cal, t_tune)

params_op_cal = occupancy(params_op_cal, t_tune)
params_op_cal = rate_scale(params_op_cal, t_tune)
params_op_cal = mscale(params_op_cal, t_tune)

params_mp_cal = rate_scale(params_mp_cal, t_tune)
params_mp_cal = mscale(params_mp_cal, t_tune)

params_o_cal = occupancy(params_o_cal, t_tune)

params_p_cal = rate_scale(params_p_cal, t_tune)
params_p_cal = mscale(params_p_cal, t_tune)

# Get effort array
effort_cal =  farray(params_bomp, fdat, seq(1991,2000), c(2,2,2,2), t_tune)

# Run calibration model
cal_bomp = project(params_bomp_cal, t_max = 100, effort = effort_cal)
cal_bom = project(params_bom_cal, t_max = 100, effort = effort_cal)
cal_bop = project(params_bop_cal, t_max = 100, effort = effort_cal)
cal_bmp = project(params_bmp_cal, t_max = 100, effort = effort_cal)
cal_omp = project(params_omp_cal, t_max = 100, effort = effort_cal)
cal_bo = project(params_bo_cal, t_max = 100, effort = effort_cal)
cal_bm = project(params_bm_cal, t_max = 100, effort = effort_cal)
cal_bp = project(params_bp_cal, t_max = 100, effort = effort_cal)
cal_om = project(params_om_cal, t_max = 100, effort = effort_cal)
cal_op = project(params_op_cal, t_max = 100, effort = effort_cal)
cal_mp = project(params_mp_cal, t_max = 100, effort = effort_cal)
cal_b = project(params_b_cal, t_max = 100, effort = effort_cal)
cal_o = project(params_o_cal, t_max = 100, effort = effort_cal)
cal_m = project(params_m_cal, t_max = 100, effort = effort_cal)
cal_p = project(params_p_cal, t_max = 100, effort = effort_cal)
cal_none = project(params_none_cal, t_max = 100, effort = effort_cal)


########## IV. Prepare and run simulations ##########

# Create an oxygen vector
t_tune = 100
times = 0:t_tune
benthic_oxygen_sim = vector(mode = "numeric", length = length(times))
pelagic_oxygen_sim = vector(mode = "numeric", length = length(times))

# Varying oxygen scenario
benthic_oxygen_sim[1:length(times)] = c(rep(mean(preds_all$df[preds_all$df$Year %in% seq(1991,2000), ]$Oxygen_benthic),t_tune-dim(preds_all$df[preds_all$df$Year %in% seq(1991,2019), ])[1]+1), preds_all$df[preds_all$df$Year %in% seq(1991,2019), ]$Oxygen_benthic)
pelagic_oxygen_sim[1:length(times)] = c(rep(mean(preds_all$df[preds_all$df$Year %in% seq(1991,2000), ]$Oxygen_pelagic),t_tune-dim(preds_all$df[preds_all$df$Year %in% seq(1991,2019), ])[1]+1), preds_all$df[preds_all$df$Year %in% seq(1991,2019), ]$Oxygen_pelagic)

# Store these in other_params
params_bomp_sim = params_bomp
params_bomp_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_bomp_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_bom_sim = params_bom
params_bom_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_bom_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_bop_sim = params_bop
params_bop_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_bop_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_bmp_sim = params_bmp
params_bmp_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_bmp_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_omp_sim = params_omp
params_omp_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_omp_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_bo_sim = params_bo
params_bo_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_bo_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_bm_sim = params_bm
params_bm_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_bm_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_bp_sim = params_bp
params_bp_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_bp_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_om_sim = params_om
params_om_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_om_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_op_sim = params_op
params_op_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_op_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_mp_sim = params_mp
params_mp_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_mp_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_b_sim = params_b
params_b_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_b_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_o_sim = params_o
params_o_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_o_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_m_sim = params_m
params_m_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_m_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_p_sim = params_p
params_p_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_p_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_none_sim = params_none
params_none_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_none_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

# Scale rates by oxygen and temperature: see otmscale.R
params_bomp_sim = occupancy(params_bomp_sim, t_tune)
params_bomp_sim = rate_scale(params_bomp_sim, t_tune)
params_bomp_sim = mscale(params_bomp_sim, t_tune)

params_bom_sim = occupancy(params_bom_sim, t_tune)

params_bop_sim = occupancy(params_bop_sim, t_tune)
params_bop_sim = rate_scale(params_bop_sim, t_tune)
params_bop_sim = mscale(params_bop_sim, t_tune)

params_bmp_sim = rate_scale(params_bmp_sim, t_tune)
params_bmp_sim = mscale(params_bmp_sim, t_tune)

params_omp_sim = occupancy(params_omp_sim, t_tune)
params_omp_sim = rate_scale(params_omp_sim, t_tune)
params_omp_sim = mscale(params_omp_sim, t_tune)

params_bo_sim = occupancy(params_bo_sim, t_tune)

params_bp_sim = rate_scale(params_bp_sim, t_tune)
params_bp_sim = mscale(params_bp_sim, t_tune)

params_om_sim = occupancy(params_om_sim, t_tune)

params_op_sim = occupancy(params_op_sim, t_tune)
params_op_sim = rate_scale(params_op_sim, t_tune)
params_op_sim = mscale(params_op_sim, t_tune)

params_mp_sim = rate_scale(params_mp_sim, t_tune)
params_mp_sim = mscale(params_mp_sim, t_tune)

params_o_sim = occupancy(params_o_sim, t_tune)

params_p_sim = rate_scale(params_p_sim, t_tune)
params_p_sim = mscale(params_p_sim, t_tune)

# Get effort array
effort_full = farray(params_bomp, fdat, seq(1991,2019), c(2,2,2,2), t_tune)

# Run full model
sim_bomp = project(params_bomp_sim, t_max = 100, effort = effort_full)
sim_bom = project(params_bom_sim, t_max = 100, effort = effort_full)
sim_bop = project(params_bop_sim, t_max = 100, effort = effort_full)
sim_bmp = project(params_bmp_sim, t_max = 100, effort = effort_full)
sim_omp = project(params_omp_sim, t_max = 100, effort = effort_full)
sim_bo = project(params_bo_sim, t_max = 100, effort = effort_full)
sim_bm = project(params_bm_sim, t_max = 100, effort = effort_full)
sim_bp = project(params_bp_sim, t_max = 100, effort = effort_full)
sim_om = project(params_om_sim, t_max = 100, effort = effort_full)
sim_op = project(params_op_sim, t_max = 100, effort = effort_full)
sim_mp = project(params_mp_sim, t_max = 100, effort = effort_full)
sim_b = project(params_b_sim, t_max = 100, effort = effort_full)
sim_o = project(params_o_sim, t_max = 100, effort = effort_full)
sim_m = project(params_m_sim, t_max = 100, effort = effort_full)
sim_p = project(params_p_sim, t_max = 100, effort = effort_full)
sim_none = project(params_none_sim, t_max = 100, effort = effort_full)


########## V. Calculate calibration error ##########

# Calibration years
cal_ts = seq(1991,2000)

# Projection years
prj_ts = seq(2001,2019)

# Model labels
model_labels = c(
	"BOMP (Full)",
	"BOM",
	"BOP",
	"BMP",
	"OMP",
	"BO",
	"BM",
	"BP",
	"OM",
	"OP",
	"MP",
	"B",
	"O",
	"M",
	"P",
	"None"
)

# SSB calibration error
cbe_cod_bomp = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bomp@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_bomp@w*params_bomp@dw, "*"))[,params_bomp@w >= species_params(params_bomp)$w_mat[1]])))^2)
cbe_cod_bom = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bom@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_bom@w*params_bom@dw, "*"))[,params_bom@w >= species_params(params_bom)$w_mat[1]])))^2)
cbe_cod_bop = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bop@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_bop@w*params_bop@dw, "*"))[,params_bop@w >= species_params(params_bop)$w_mat[1]])))^2)
cbe_cod_bmp = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bmp@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_bmp@w*params_bmp@dw, "*"))[,params_bmp@w >= species_params(params_bmp)$w_mat[1]])))^2)
cbe_cod_omp = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_omp@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_omp@w*params_omp@dw, "*"))[,params_omp@w >= species_params(params_omp)$w_mat[1]])))^2)
cbe_cod_bo = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bo@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_bo@w*params_bo@dw, "*"))[,params_bo@w >= species_params(params_bo)$w_mat[1]])))^2)
cbe_cod_bm = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bm@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_bm@w*params_bm@dw, "*"))[,params_bm@w >= species_params(params_bm)$w_mat[1]])))^2)
cbe_cod_bp = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bp@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_bp@w*params_bp@dw, "*"))[,params_bp@w >= species_params(params_bp)$w_mat[1]])))^2)
cbe_cod_om = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_om@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_om@w*params_om@dw, "*"))[,params_om@w >= species_params(params_om)$w_mat[1]])))^2)
cbe_cod_op = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_op@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_op@w*params_op@dw, "*"))[,params_op@w >= species_params(params_op)$w_mat[1]])))^2)
cbe_cod_mp = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_mp@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_mp@w*params_mp@dw, "*"))[,params_mp@w >= species_params(params_mp)$w_mat[1]])))^2)
cbe_cod_b = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_b@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_b@w*params_b@dw, "*"))[,params_b@w >= species_params(params_b)$w_mat[1]])))^2)
cbe_cod_o = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_o@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_o@w*params_o@dw, "*"))[,params_o@w >= species_params(params_o)$w_mat[1]])))^2)
cbe_cod_m = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_m@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_m@w*params_m@dw, "*"))[,params_m@w >= species_params(params_m)$w_mat[1]])))^2)
cbe_cod_p = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_p@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_p@w*params_p@dw, "*"))[,params_p@w >= species_params(params_p)$w_mat[1]])))^2)
cbe_cod_none = sum((log(cod_ssb[cod_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_none@n[(t_tune-length(cal_ts)+1):t_tune,1,], 2, params_none@w*params_none@dw, "*"))[,params_none@w >= species_params(params_none)$w_mat[1]])))^2)

cbe_flounder_bomp = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bomp@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_bomp@w*params_bomp@dw, "*"))[,params_bomp@w >= species_params(params_bomp)$w_mat[2]])))^2)
cbe_flounder_bom = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bom@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_bom@w*params_bom@dw, "*"))[,params_bom@w >= species_params(params_bom)$w_mat[2]])))^2)
cbe_flounder_bop = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bop@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_bop@w*params_bop@dw, "*"))[,params_bop@w >= species_params(params_bop)$w_mat[2]])))^2)
cbe_flounder_bmp = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bmp@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_bmp@w*params_bmp@dw, "*"))[,params_bmp@w >= species_params(params_bmp)$w_mat[2]])))^2)
cbe_flounder_omp = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_omp@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_omp@w*params_omp@dw, "*"))[,params_omp@w >= species_params(params_omp)$w_mat[2]])))^2)
cbe_flounder_bo = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bo@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_bo@w*params_bo@dw, "*"))[,params_bo@w >= species_params(params_bo)$w_mat[2]])))^2)
cbe_flounder_bm = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bm@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_bm@w*params_bm@dw, "*"))[,params_bm@w >= species_params(params_bm)$w_mat[2]])))^2)
cbe_flounder_bp = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bp@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_bp@w*params_bp@dw, "*"))[,params_bp@w >= species_params(params_bp)$w_mat[2]])))^2)
cbe_flounder_om = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_om@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_om@w*params_om@dw, "*"))[,params_om@w >= species_params(params_om)$w_mat[2]])))^2)
cbe_flounder_op = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_op@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_op@w*params_op@dw, "*"))[,params_op@w >= species_params(params_op)$w_mat[2]])))^2)
cbe_flounder_mp = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_mp@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_mp@w*params_mp@dw, "*"))[,params_mp@w >= species_params(params_mp)$w_mat[2]])))^2)
cbe_flounder_b = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_b@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_b@w*params_b@dw, "*"))[,params_b@w >= species_params(params_b)$w_mat[2]])))^2)
cbe_flounder_o = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_o@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_o@w*params_o@dw, "*"))[,params_o@w >= species_params(params_o)$w_mat[2]])))^2)
cbe_flounder_m = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_m@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_m@w*params_m@dw, "*"))[,params_m@w >= species_params(params_m)$w_mat[2]])))^2)
cbe_flounder_p = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_p@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_p@w*params_p@dw, "*"))[,params_p@w >= species_params(params_p)$w_mat[2]])))^2)
cbe_flounder_none = sum((log(flounder_ssb[flounder_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_none@n[(t_tune-length(cal_ts)+1):t_tune,2,], 2, params_none@w*params_none@dw, "*"))[,params_none@w >= species_params(params_none)$w_mat[2]])))^2)

cbe_sprat_bomp = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bomp@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_bomp@w*params_bomp@dw, "*"))[,params_bomp@w >= species_params(params_bomp)$w_mat[3]])))^2)
cbe_sprat_bom = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bom@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_bom@w*params_bom@dw, "*"))[,params_bom@w >= species_params(params_bom)$w_mat[3]])))^2)
cbe_sprat_bop = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bop@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_bop@w*params_bop@dw, "*"))[,params_bop@w >= species_params(params_bop)$w_mat[3]])))^2)
cbe_sprat_bmp = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bmp@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_bmp@w*params_bmp@dw, "*"))[,params_bmp@w >= species_params(params_bmp)$w_mat[3]])))^2)
cbe_sprat_omp = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_omp@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_omp@w*params_omp@dw, "*"))[,params_omp@w >= species_params(params_omp)$w_mat[3]])))^2)
cbe_sprat_bo = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bo@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_bo@w*params_bo@dw, "*"))[,params_bo@w >= species_params(params_bo)$w_mat[3]])))^2)
cbe_sprat_bm = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bm@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_bm@w*params_bm@dw, "*"))[,params_bm@w >= species_params(params_bm)$w_mat[3]])))^2)
cbe_sprat_bp = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bp@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_bp@w*params_bp@dw, "*"))[,params_bp@w >= species_params(params_bp)$w_mat[3]])))^2)
cbe_sprat_om = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_om@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_om@w*params_om@dw, "*"))[,params_om@w >= species_params(params_om)$w_mat[3]])))^2)
cbe_sprat_op = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_op@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_op@w*params_op@dw, "*"))[,params_op@w >= species_params(params_op)$w_mat[3]])))^2)
cbe_sprat_mp = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_mp@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_mp@w*params_mp@dw, "*"))[,params_mp@w >= species_params(params_mp)$w_mat[3]])))^2)
cbe_sprat_b = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_b@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_b@w*params_b@dw, "*"))[,params_b@w >= species_params(params_b)$w_mat[3]])))^2)
cbe_sprat_o = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_o@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_o@w*params_o@dw, "*"))[,params_o@w >= species_params(params_o)$w_mat[3]])))^2)
cbe_sprat_m = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_m@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_m@w*params_m@dw, "*"))[,params_m@w >= species_params(params_m)$w_mat[3]])))^2)
cbe_sprat_p = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_p@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_p@w*params_p@dw, "*"))[,params_p@w >= species_params(params_p)$w_mat[3]])))^2)
cbe_sprat_none = sum((log(sprat_ssb[sprat_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_none@n[(t_tune-length(cal_ts)+1):t_tune,3,], 2, params_none@w*params_none@dw, "*"))[,params_none@w >= species_params(params_none)$w_mat[3]])))^2)

cbe_herring_bomp = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bomp@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_bomp@w*params_bomp@dw, "*"))[,params_bomp@w >= species_params(params_bomp)$w_mat[4]])))^2)
cbe_herring_bom = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bom@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_bom@w*params_bom@dw, "*"))[,params_bom@w >= species_params(params_bom)$w_mat[4]])))^2)
cbe_herring_bop = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bop@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_bop@w*params_bop@dw, "*"))[,params_bop@w >= species_params(params_bop)$w_mat[4]])))^2)
cbe_herring_bmp = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bmp@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_bmp@w*params_bmp@dw, "*"))[,params_bmp@w >= species_params(params_bmp)$w_mat[4]])))^2)
cbe_herring_omp = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_omp@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_omp@w*params_omp@dw, "*"))[,params_omp@w >= species_params(params_omp)$w_mat[4]])))^2)
cbe_herring_bo = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bo@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_bo@w*params_bo@dw, "*"))[,params_bo@w >= species_params(params_bo)$w_mat[4]])))^2)
cbe_herring_bm = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bm@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_bm@w*params_bm@dw, "*"))[,params_bm@w >= species_params(params_bm)$w_mat[4]])))^2)
cbe_herring_bp = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_bp@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_bp@w*params_bp@dw, "*"))[,params_bp@w >= species_params(params_bp)$w_mat[4]])))^2)
cbe_herring_om = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_om@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_om@w*params_om@dw, "*"))[,params_om@w >= species_params(params_om)$w_mat[4]])))^2)
cbe_herring_op = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_op@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_op@w*params_op@dw, "*"))[,params_op@w >= species_params(params_op)$w_mat[4]])))^2)
cbe_herring_mp = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_mp@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_mp@w*params_mp@dw, "*"))[,params_mp@w >= species_params(params_mp)$w_mat[4]])))^2)
cbe_herring_b = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_b@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_b@w*params_b@dw, "*"))[,params_b@w >= species_params(params_b)$w_mat[4]])))^2)
cbe_herring_o = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_o@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_o@w*params_o@dw, "*"))[,params_o@w >= species_params(params_o)$w_mat[4]])))^2)
cbe_herring_m = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_m@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_m@w*params_m@dw, "*"))[,params_m@w >= species_params(params_m)$w_mat[4]])))^2)
cbe_herring_p = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_p@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_p@w*params_p@dw, "*"))[,params_p@w >= species_params(params_p)$w_mat[4]])))^2)
cbe_herring_none = sum((log(herring_ssb[herring_ssb$Year %in% cal_ts, ]$SSB) - log(rowSums((sweep(cal_none@n[(t_tune-length(cal_ts)+1):t_tune,4,], 2, params_none@w*params_none@dw, "*"))[,params_none@w >= species_params(params_none)$w_mat[4]])))^2)

# Sum across species
cbe_cod = c(cbe_cod_bomp, cbe_cod_bom, cbe_cod_bop, cbe_cod_bmp, cbe_cod_omp, cbe_cod_bo, cbe_cod_bm, cbe_cod_bp, cbe_cod_om, cbe_cod_op, cbe_cod_mp, cbe_cod_b, cbe_cod_o, cbe_cod_m, cbe_cod_p, cbe_cod_none)
cbe_flounder = c(cbe_flounder_bomp, cbe_flounder_bom, cbe_flounder_bop, cbe_flounder_bmp, cbe_flounder_omp, cbe_flounder_bo, cbe_flounder_bm, cbe_flounder_bp, cbe_flounder_om, cbe_flounder_op, cbe_flounder_mp, cbe_flounder_b, cbe_flounder_o, cbe_flounder_m, cbe_flounder_p, cbe_flounder_none)
cbe_sprat = c(cbe_sprat_bomp, cbe_sprat_bom, cbe_sprat_bop, cbe_sprat_bmp, cbe_sprat_omp, cbe_sprat_bo, cbe_sprat_bm, cbe_sprat_bp, cbe_sprat_om, cbe_sprat_op, cbe_sprat_mp, cbe_sprat_b, cbe_sprat_o, cbe_sprat_m, cbe_sprat_p, cbe_sprat_none)
cbe_herring = c(cbe_herring_bomp, cbe_herring_bom, cbe_herring_bop, cbe_herring_bmp, cbe_herring_omp, cbe_herring_bo, cbe_herring_bm, cbe_herring_bp, cbe_herring_om, cbe_herring_op, cbe_herring_mp, cbe_herring_b, cbe_herring_o, cbe_herring_m, cbe_herring_p, cbe_herring_none)

# Total and ranked error
cbe_cod_weight = 0.3
cbe_flounder_weight = 0.1 
cbe_sprat_weight = 0.3
cbe_herring_weight = 0.3
cbe_error = cbe_cod_weight*cbe_cod + cbe_flounder_weight*cbe_flounder + cbe_sprat_weight*cbe_sprat + cbe_herring_weight*cbe_herring
cbe_rank = (cbe_error - min(cbe_error)) / (max(cbe_error) - min(cbe_error))

# Yield calibration error
yield_bomp = getYield(cal_bomp)
yield_bom = getYield(cal_bom)
yield_bop = getYield(cal_bop)
yield_bmp = getYield(cal_bmp)
yield_omp = getYield(cal_omp)
yield_bo = getYield(cal_bo)
yield_bm = getYield(cal_bm)
yield_bp = getYield(cal_bp)
yield_om = getYield(cal_om)
yield_op = getYield(cal_op)
yield_mp = getYield(cal_mp)
yield_b = getYield(cal_b)
yield_o = getYield(cal_o)
yield_m = getYield(cal_m)
yield_p = getYield(cal_p)
yield_none = getYield(cal_none)

cye_cod_bomp = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bomp[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_bom = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bom[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_bop = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bop[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_bmp = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bmp[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_omp = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_omp[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_bo = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bo[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_bm = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bm[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_bp = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bp[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_om = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_om[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_op = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_op[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_mp = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_mp[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_b = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_b[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_o = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_o[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_m = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_m[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_p = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_p[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)
cye_cod_none = sum((log(cod_catch[cod_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_none[(t_tune-length(cal_ts)+1):t_tune, "cod"]))^2)

cye_flounder_bomp = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bomp[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_bom = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bom[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_bop = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bop[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_bmp = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bmp[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_omp = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_omp[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_bo = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bo[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_bm = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bm[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_bp = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bp[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_om = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_om[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_op = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_op[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_mp = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_mp[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_b = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_b[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_o = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_o[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_m = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_m[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_p = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_p[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)
cye_flounder_none = sum((log(flounder_catch[flounder_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_none[(t_tune-length(cal_ts)+1):t_tune, "flounder"]))^2)

cye_sprat_bomp = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bomp[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_bom = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bom[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_bop = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bop[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_bmp = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bmp[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_omp = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_omp[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_bo = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bo[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_bm = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bm[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_bp = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bp[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_om = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_om[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_op = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_op[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_mp = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_mp[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_b = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_b[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_o = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_o[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_m = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_m[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_p = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_p[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)
cye_sprat_none = sum((log(sprat_catch[sprat_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_none[(t_tune-length(cal_ts)+1):t_tune, "sprat"]))^2)

cye_herring_bomp = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bomp[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_bom = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bom[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_bop = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bop[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_bmp = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bmp[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_omp = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_omp[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_bo = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bo[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_bm = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bm[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_bp = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_bp[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_om = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_om[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_op = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_op[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_mp = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_mp[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_b = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_b[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_o = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_o[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_m = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_m[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_p = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_p[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)
cye_herring_none = sum((log(herring_catch[herring_catch$Year %in% cal_ts, ]$Catch*1000*1000) - log(yield_none[(t_tune-length(cal_ts)+1):t_tune, "herring"]))^2)

# Sum across species
cye_cod = c(cye_cod_bomp, cye_cod_bom, cye_cod_bop, cye_cod_bmp, cye_cod_omp, cye_cod_bo, cye_cod_bm, cye_cod_bp, cye_cod_om, cye_cod_op, cye_cod_mp, cye_cod_b, cye_cod_o, cye_cod_m, cye_cod_p, cye_cod_none)
cye_flounder = c(cye_flounder_bomp, cye_flounder_bom, cye_flounder_bop, cye_flounder_bmp, cye_flounder_omp, cye_flounder_bo, cye_flounder_bm, cye_flounder_bp, cye_flounder_om, cye_flounder_op, cye_flounder_mp, cye_flounder_b, cye_flounder_o, cye_flounder_m, cye_flounder_p, cye_flounder_none)
cye_sprat = c(cye_sprat_bomp, cye_sprat_bom, cye_sprat_bop, cye_sprat_bmp, cye_sprat_omp, cye_sprat_bo, cye_sprat_bm, cye_sprat_bp, cye_sprat_om, cye_sprat_op, cye_sprat_mp, cye_sprat_b, cye_sprat_o, cye_sprat_m, cye_sprat_p, cye_sprat_none)
cye_herring = c(cye_herring_bomp, cye_herring_bom, cye_herring_bop, cye_herring_bmp, cye_herring_omp, cye_herring_bo, cye_herring_bm, cye_herring_bp, cye_herring_om, cye_herring_op, cye_herring_mp, cye_herring_b, cye_herring_o, cye_herring_m, cye_herring_p, cye_herring_none)

# Total and ranked error
cye_cod_weight = 0.25
cye_flounder_weight = 0.25 
cye_sprat_weight = 0.25
cye_herring_weight = 0.25
cye_error = cye_cod_weight*cye_cod + cye_flounder_weight*cye_flounder + cye_sprat_weight*cye_sprat + cye_herring_weight*cye_herring
cye_rank = (cye_error - min(cye_error)) / (max(cye_error) - min(cye_error))

# Growth calibration error
cge_cod_bomp = proj_growth(params_bomp, cal_ts, 100, "cod")
cge_cod_bom = proj_growth(params_bom, cal_ts, 100, "cod")
cge_cod_bop = proj_growth(params_bop, cal_ts, 100, "cod")
cge_cod_bmp = proj_growth(params_bmp, cal_ts, 100, "cod")
cge_cod_omp = proj_growth(params_omp, cal_ts, 100, "cod")
cge_cod_bo = proj_growth(params_bo, cal_ts, 100, "cod")
cge_cod_bm = proj_growth(params_bm, cal_ts, 100, "cod")
cge_cod_bp = proj_growth(params_bp, cal_ts, 100, "cod")
cge_cod_om = proj_growth(params_om, cal_ts, 100, "cod")
cge_cod_op = proj_growth(params_op, cal_ts, 100, "cod")
cge_cod_mp = proj_growth(params_mp, cal_ts, 100, "cod")
cge_cod_b = proj_growth(params_b, cal_ts, 100, "cod")
cge_cod_o = proj_growth(params_o, cal_ts, 100, "cod")
cge_cod_m = proj_growth(params_m, cal_ts, 100, "cod")
cge_cod_p = proj_growth(params_p, cal_ts, 100, "cod")
cge_cod_none = proj_growth(params_none, cal_ts, 100, "cod")

cge_flounder_bomp = proj_growth(params_bomp, cal_ts, 100, "flounder")
cge_flounder_bom = proj_growth(params_bom, cal_ts, 100, "flounder")
cge_flounder_bop = proj_growth(params_bop, cal_ts, 100, "flounder")
cge_flounder_bmp = proj_growth(params_bmp, cal_ts, 100, "flounder")
cge_flounder_omp = proj_growth(params_omp, cal_ts, 100, "flounder")
cge_flounder_bo = proj_growth(params_bo, cal_ts, 100, "flounder")
cge_flounder_bm = proj_growth(params_bm, cal_ts, 100, "flounder")
cge_flounder_bp = proj_growth(params_bp, cal_ts, 100, "flounder")
cge_flounder_om = proj_growth(params_om, cal_ts, 100, "flounder")
cge_flounder_op = proj_growth(params_op, cal_ts, 100, "flounder")
cge_flounder_mp = proj_growth(params_mp, cal_ts, 100, "flounder")
cge_flounder_b = proj_growth(params_b, cal_ts, 100, "flounder")
cge_flounder_o = proj_growth(params_o, cal_ts, 100, "flounder")
cge_flounder_m = proj_growth(params_m, cal_ts, 100, "flounder")
cge_flounder_p = proj_growth(params_p, cal_ts, 100, "flounder")
cge_flounder_none = proj_growth(params_none, cal_ts, 100, "flounder")

cge_sprat_bomp = proj_growth(params_bomp, cal_ts, 100, "sprat")
cge_sprat_bom = proj_growth(params_bom, cal_ts, 100, "sprat")
cge_sprat_bop = proj_growth(params_bop, cal_ts, 100, "sprat")
cge_sprat_bmp = proj_growth(params_bmp, cal_ts, 100, "sprat")
cge_sprat_omp = proj_growth(params_omp, cal_ts, 100, "sprat")
cge_sprat_bo = proj_growth(params_bo, cal_ts, 100, "sprat")
cge_sprat_bm = proj_growth(params_bm, cal_ts, 100, "sprat")
cge_sprat_bp = proj_growth(params_bp, cal_ts, 100, "sprat")
cge_sprat_om = proj_growth(params_om, cal_ts, 100, "sprat")
cge_sprat_op = proj_growth(params_op, cal_ts, 100, "sprat")
cge_sprat_mp = proj_growth(params_mp, cal_ts, 100, "sprat")
cge_sprat_b = proj_growth(params_b, cal_ts, 100, "sprat")
cge_sprat_o = proj_growth(params_o, cal_ts, 100, "sprat")
cge_sprat_m = proj_growth(params_m, cal_ts, 100, "sprat")
cge_sprat_p = proj_growth(params_p, cal_ts, 100, "sprat")
cge_sprat_none = proj_growth(params_none, cal_ts, 100, "sprat")

cge_herring_bomp = proj_growth(params_bomp, cal_ts, 100, "herring")
cge_herring_bom = proj_growth(params_bom, cal_ts, 100, "herring")
cge_herring_bop = proj_growth(params_bop, cal_ts, 100, "herring")
cge_herring_bmp = proj_growth(params_bmp, cal_ts, 100, "herring")
cge_herring_omp = proj_growth(params_omp, cal_ts, 100, "herring")
cge_herring_bo = proj_growth(params_bo, cal_ts, 100, "herring")
cge_herring_bm = proj_growth(params_bm, cal_ts, 100, "herring")
cge_herring_bp = proj_growth(params_bp, cal_ts, 100, "herring")
cge_herring_om = proj_growth(params_om, cal_ts, 100, "herring")
cge_herring_op = proj_growth(params_op, cal_ts, 100, "herring")
cge_herring_mp = proj_growth(params_mp, cal_ts, 100, "herring")
cge_herring_b = proj_growth(params_b, cal_ts, 100, "herring")
cge_herring_o = proj_growth(params_o, cal_ts, 100, "herring")
cge_herring_m = proj_growth(params_m, cal_ts, 100, "herring")
cge_herring_p = proj_growth(params_p, cal_ts, 100, "herring")
cge_herring_none = proj_growth(params_none, cal_ts, 100, "herring")

# Sum across species
cge_cod = c(cge_cod_bomp, cge_cod_bom, cge_cod_bop, cge_cod_bmp, cge_cod_omp, cge_cod_bo, cge_cod_bm, cge_cod_bp, cge_cod_om, cge_cod_op, cge_cod_mp, cge_cod_b, cge_cod_o, cge_cod_m, cge_cod_p, cge_cod_none)
cge_flounder = c(cge_flounder_bomp, cge_flounder_bom, cge_flounder_bop, cge_flounder_bmp, cge_flounder_omp, cge_flounder_bo, cge_flounder_bm, cge_flounder_bp, cge_flounder_om, cge_flounder_op, cge_flounder_mp, cge_flounder_b, cge_flounder_o, cge_flounder_m, cge_flounder_p, cge_flounder_none)
cge_sprat = c(cge_sprat_bomp, cge_sprat_bom, cge_sprat_bop, cge_sprat_bmp, cge_sprat_omp, cge_sprat_bo, cge_sprat_bm, cge_sprat_bp, cge_sprat_om, cge_sprat_op, cge_sprat_mp, cge_sprat_b, cge_sprat_o, cge_sprat_m, cge_sprat_p, cge_sprat_none)
cge_herring = c(cge_herring_bomp, cge_herring_bom, cge_herring_bop, cge_herring_bmp, cge_herring_omp, cge_herring_bo, cge_herring_bm, cge_herring_bp, cge_herring_om, cge_herring_op, cge_herring_mp, cge_herring_b, cge_herring_o, cge_herring_m, cge_herring_p, cge_herring_none)

# Total and ranked error
cge_cod_weight = 0.25
cge_flounder_weight = 0.25
cge_sprat_weight = 0.25
cge_herring_weight = 0.25
cge_error = cge_cod_weight*cge_cod + cge_flounder_weight*cge_flounder + cge_sprat_weight*cge_sprat + cge_herring_weight*cge_herring
cge_rank = (cge_error - min(cge_error)) / (max(cge_error) - min(cge_error))

# Calculate final weighted error
ssb_weight = 0.3
yield_weight = 0.1
growth_weight = 0.6
final_error_cal = ssb_weight*cbe_rank + yield_weight*cye_rank + growth_weight*cge_rank
model_labels[order(final_error_cal)]

# Create and save error table
cal_error = data.frame(model = model_labels)
cal_error$SSB = ssb_weight*cbe_rank
cal_error$yield = yield_weight*cye_rank
cal_error$growth = ssb_weight*cge_rank
cal_error$rank = cal_error$SSB + cal_error$yield + cal_error$growth
cal_error$SSB_error = cbe_error
cal_error$yield_error = cye_error
cal_error$growth_error = cge_error
write.table(cal_error, "Calibration/cal_error.txt")


########## VI. Calculation projection error ##########

# SSB projection error
pbe_cod_bomp = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bomp@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_bomp@w*params_bomp@dw, "*"))[,params_bomp@w >= species_params(params_bomp)$w_mat[1]])))^2)
pbe_cod_bom = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bom@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_bom@w*params_bom@dw, "*"))[,params_bom@w >= species_params(params_bom)$w_mat[1]])))^2)
pbe_cod_bop = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bop@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_bop@w*params_bop@dw, "*"))[,params_bop@w >= species_params(params_bop)$w_mat[1]])))^2)
pbe_cod_bmp = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bmp@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_bmp@w*params_bmp@dw, "*"))[,params_bmp@w >= species_params(params_bmp)$w_mat[1]])))^2)
pbe_cod_omp = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_omp@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_omp@w*params_omp@dw, "*"))[,params_omp@w >= species_params(params_omp)$w_mat[1]])))^2)
pbe_cod_bo = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bo@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_bo@w*params_bo@dw, "*"))[,params_bo@w >= species_params(params_bo)$w_mat[1]])))^2)
pbe_cod_bm = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bm@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_bm@w*params_bm@dw, "*"))[,params_bm@w >= species_params(params_bm)$w_mat[1]])))^2)
pbe_cod_bp = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bp@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_bp@w*params_bp@dw, "*"))[,params_bp@w >= species_params(params_bp)$w_mat[1]])))^2)
pbe_cod_om = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_om@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_om@w*params_om@dw, "*"))[,params_om@w >= species_params(params_om)$w_mat[1]])))^2)
pbe_cod_op = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_op@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_op@w*params_op@dw, "*"))[,params_op@w >= species_params(params_op)$w_mat[1]])))^2)
pbe_cod_mp = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_mp@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_mp@w*params_mp@dw, "*"))[,params_mp@w >= species_params(params_mp)$w_mat[1]])))^2)
pbe_cod_b = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_b@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_b@w*params_b@dw, "*"))[,params_b@w >= species_params(params_b)$w_mat[1]])))^2)
pbe_cod_o = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_o@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_o@w*params_o@dw, "*"))[,params_o@w >= species_params(params_o)$w_mat[1]])))^2)
pbe_cod_m = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_m@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_m@w*params_m@dw, "*"))[,params_m@w >= species_params(params_m)$w_mat[1]])))^2)
pbe_cod_p = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_p@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_p@w*params_p@dw, "*"))[,params_p@w >= species_params(params_p)$w_mat[1]])))^2)
pbe_cod_none = sum((log(cod_ssb[cod_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_none@n[(t_tune-length(prj_ts)+1):t_tune,1,], 2, params_none@w*params_none@dw, "*"))[,params_none@w >= species_params(params_none)$w_mat[1]])))^2)

pbe_flounder_bomp = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bomp@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_bomp@w*params_bomp@dw, "*"))[,params_bomp@w >= species_params(params_bomp)$w_mat[2]])))^2)
pbe_flounder_bom = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bom@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_bom@w*params_bom@dw, "*"))[,params_bom@w >= species_params(params_bom)$w_mat[2]])))^2)
pbe_flounder_bop = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bop@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_bop@w*params_bop@dw, "*"))[,params_bop@w >= species_params(params_bop)$w_mat[2]])))^2)
pbe_flounder_bmp = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bmp@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_bmp@w*params_bmp@dw, "*"))[,params_bmp@w >= species_params(params_bmp)$w_mat[2]])))^2)
pbe_flounder_omp = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_omp@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_omp@w*params_omp@dw, "*"))[,params_omp@w >= species_params(params_omp)$w_mat[2]])))^2)
pbe_flounder_bo = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bo@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_bo@w*params_bo@dw, "*"))[,params_bo@w >= species_params(params_bo)$w_mat[2]])))^2)
pbe_flounder_bm = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bm@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_bm@w*params_bm@dw, "*"))[,params_bm@w >= species_params(params_bm)$w_mat[2]])))^2)
pbe_flounder_bp = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bp@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_bp@w*params_bp@dw, "*"))[,params_bp@w >= species_params(params_bp)$w_mat[2]])))^2)
pbe_flounder_om = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_om@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_om@w*params_om@dw, "*"))[,params_om@w >= species_params(params_om)$w_mat[2]])))^2)
pbe_flounder_op = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_op@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_op@w*params_op@dw, "*"))[,params_op@w >= species_params(params_op)$w_mat[2]])))^2)
pbe_flounder_mp = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_mp@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_mp@w*params_mp@dw, "*"))[,params_mp@w >= species_params(params_mp)$w_mat[2]])))^2)
pbe_flounder_b = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_b@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_b@w*params_b@dw, "*"))[,params_b@w >= species_params(params_b)$w_mat[2]])))^2)
pbe_flounder_o = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_o@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_o@w*params_o@dw, "*"))[,params_o@w >= species_params(params_o)$w_mat[2]])))^2)
pbe_flounder_m = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_m@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_m@w*params_m@dw, "*"))[,params_m@w >= species_params(params_m)$w_mat[2]])))^2)
pbe_flounder_p = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_p@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_p@w*params_p@dw, "*"))[,params_p@w >= species_params(params_p)$w_mat[2]])))^2)
pbe_flounder_none = sum((log(flounder_ssb[flounder_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_none@n[(t_tune-length(prj_ts)+1):t_tune,2,], 2, params_none@w*params_none@dw, "*"))[,params_none@w >= species_params(params_none)$w_mat[2]])))^2)

pbe_sprat_bomp = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bomp@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_bomp@w*params_bomp@dw, "*"))[,params_bomp@w >= species_params(params_bomp)$w_mat[3]])))^2)
pbe_sprat_bom = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bom@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_bom@w*params_bom@dw, "*"))[,params_bom@w >= species_params(params_bom)$w_mat[3]])))^2)
pbe_sprat_bop = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bop@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_bop@w*params_bop@dw, "*"))[,params_bop@w >= species_params(params_bop)$w_mat[3]])))^2)
pbe_sprat_bmp = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bmp@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_bmp@w*params_bmp@dw, "*"))[,params_bmp@w >= species_params(params_bmp)$w_mat[3]])))^2)
pbe_sprat_omp = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_omp@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_omp@w*params_omp@dw, "*"))[,params_omp@w >= species_params(params_omp)$w_mat[3]])))^2)
pbe_sprat_bo = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bo@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_bo@w*params_bo@dw, "*"))[,params_bo@w >= species_params(params_bo)$w_mat[3]])))^2)
pbe_sprat_bm = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bm@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_bm@w*params_bm@dw, "*"))[,params_bm@w >= species_params(params_bm)$w_mat[3]])))^2)
pbe_sprat_bp = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bp@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_bp@w*params_bp@dw, "*"))[,params_bp@w >= species_params(params_bp)$w_mat[3]])))^2)
pbe_sprat_om = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_om@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_om@w*params_om@dw, "*"))[,params_om@w >= species_params(params_om)$w_mat[3]])))^2)
pbe_sprat_op = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_op@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_op@w*params_op@dw, "*"))[,params_op@w >= species_params(params_op)$w_mat[3]])))^2)
pbe_sprat_mp = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_mp@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_mp@w*params_mp@dw, "*"))[,params_mp@w >= species_params(params_mp)$w_mat[3]])))^2)
pbe_sprat_b = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_b@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_b@w*params_b@dw, "*"))[,params_b@w >= species_params(params_b)$w_mat[3]])))^2)
pbe_sprat_o = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_o@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_o@w*params_o@dw, "*"))[,params_o@w >= species_params(params_o)$w_mat[3]])))^2)
pbe_sprat_m = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_m@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_m@w*params_m@dw, "*"))[,params_m@w >= species_params(params_m)$w_mat[3]])))^2)
pbe_sprat_p = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_p@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_p@w*params_p@dw, "*"))[,params_p@w >= species_params(params_p)$w_mat[3]])))^2)
pbe_sprat_none = sum((log(sprat_ssb[sprat_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_none@n[(t_tune-length(prj_ts)+1):t_tune,3,], 2, params_none@w*params_none@dw, "*"))[,params_none@w >= species_params(params_none)$w_mat[3]])))^2)

pbe_herring_bomp = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bomp@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_bomp@w*params_bomp@dw, "*"))[,params_bomp@w >= species_params(params_bomp)$w_mat[4]])))^2)
pbe_herring_bom = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bom@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_bom@w*params_bom@dw, "*"))[,params_bom@w >= species_params(params_bom)$w_mat[4]])))^2)
pbe_herring_bop = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bop@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_bop@w*params_bop@dw, "*"))[,params_bop@w >= species_params(params_bop)$w_mat[4]])))^2)
pbe_herring_bmp = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bmp@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_bmp@w*params_bmp@dw, "*"))[,params_bmp@w >= species_params(params_bmp)$w_mat[4]])))^2)
pbe_herring_omp = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_omp@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_omp@w*params_omp@dw, "*"))[,params_omp@w >= species_params(params_omp)$w_mat[4]])))^2)
pbe_herring_bo = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bo@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_bo@w*params_bo@dw, "*"))[,params_bo@w >= species_params(params_bo)$w_mat[4]])))^2)
pbe_herring_bm = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bm@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_bm@w*params_bm@dw, "*"))[,params_bm@w >= species_params(params_bm)$w_mat[4]])))^2)
pbe_herring_bp = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_bp@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_bp@w*params_bp@dw, "*"))[,params_bp@w >= species_params(params_bp)$w_mat[4]])))^2)
pbe_herring_om = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_om@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_om@w*params_om@dw, "*"))[,params_om@w >= species_params(params_om)$w_mat[4]])))^2)
pbe_herring_op = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_op@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_op@w*params_op@dw, "*"))[,params_op@w >= species_params(params_op)$w_mat[4]])))^2)
pbe_herring_mp = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_mp@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_mp@w*params_mp@dw, "*"))[,params_mp@w >= species_params(params_mp)$w_mat[4]])))^2)
pbe_herring_b = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_b@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_b@w*params_b@dw, "*"))[,params_b@w >= species_params(params_b)$w_mat[4]])))^2)
pbe_herring_o = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_o@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_o@w*params_o@dw, "*"))[,params_o@w >= species_params(params_o)$w_mat[4]])))^2)
pbe_herring_m = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_m@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_m@w*params_m@dw, "*"))[,params_m@w >= species_params(params_m)$w_mat[4]])))^2)
pbe_herring_p = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_p@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_p@w*params_p@dw, "*"))[,params_p@w >= species_params(params_p)$w_mat[4]])))^2)
pbe_herring_none = sum((log(herring_ssb[herring_ssb$Year %in% prj_ts, ]$SSB) - log(rowSums((sweep(sim_none@n[(t_tune-length(prj_ts)+1):t_tune,4,], 2, params_none@w*params_none@dw, "*"))[,params_none@w >= species_params(params_none)$w_mat[4]])))^2)

# Sum across species
pbe_cod = c(pbe_cod_bomp, pbe_cod_bom, pbe_cod_bop, pbe_cod_bmp, pbe_cod_omp, pbe_cod_bo, pbe_cod_bm, pbe_cod_bp, pbe_cod_om, pbe_cod_op, pbe_cod_mp, pbe_cod_b, pbe_cod_o, pbe_cod_m, pbe_cod_p, pbe_cod_none)
pbe_flounder = c(pbe_flounder_bomp, pbe_flounder_bom, pbe_flounder_bop, pbe_flounder_bmp, pbe_flounder_omp, pbe_flounder_bo, pbe_flounder_bm, pbe_flounder_bp, pbe_flounder_om, pbe_flounder_op, pbe_flounder_mp, pbe_flounder_b, pbe_flounder_o, pbe_flounder_m, pbe_flounder_p, pbe_flounder_none)
pbe_sprat = c(pbe_sprat_bomp, pbe_sprat_bom, pbe_sprat_bop, pbe_sprat_bmp, pbe_sprat_omp, pbe_sprat_bo, pbe_sprat_bm, pbe_sprat_bp, pbe_sprat_om, pbe_sprat_op, pbe_sprat_mp, pbe_sprat_b, pbe_sprat_o, pbe_sprat_m, pbe_sprat_p, pbe_sprat_none)
pbe_herring = c(pbe_herring_bomp, pbe_herring_bom, pbe_herring_bop, pbe_herring_bmp, pbe_herring_omp, pbe_herring_bo, pbe_herring_bm, pbe_herring_bp, pbe_herring_om, pbe_herring_op, pbe_herring_mp, pbe_herring_b, pbe_herring_o, pbe_herring_m, pbe_herring_p, pbe_herring_none)

# Total and ranked error
pbe_cod_weight = 0.3
pbe_flounder_weight = 0.1 
pbe_sprat_weight = 0.3
pbe_herring_weight = 0.3
pbe_error = pbe_cod_weight*pbe_cod + pbe_flounder_weight*pbe_flounder + pbe_sprat_weight*pbe_sprat + pbe_herring_weight*pbe_herring
pbe_rank = (pbe_error - min(pbe_error)) / (max(pbe_error) - min(pbe_error))

# Yield projection error
yield_bomp = getYield(sim_bomp)
yield_bom = getYield(sim_bom)
yield_bop = getYield(sim_bop)
yield_bmp = getYield(sim_bmp)
yield_omp = getYield(sim_omp)
yield_bo = getYield(sim_bo)
yield_bm = getYield(sim_bm)
yield_bp = getYield(sim_bp)
yield_om = getYield(sim_om)
yield_op = getYield(sim_op)
yield_mp = getYield(sim_mp)
yield_b = getYield(sim_b)
yield_o = getYield(sim_o)
yield_m = getYield(sim_m)
yield_p = getYield(sim_p)
yield_none = getYield(sim_none)

pye_cod_bomp = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bomp[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_bom = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bom[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_bop = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bop[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_bmp = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bmp[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_omp = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_omp[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_bo = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bo[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_bm = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bm[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_bp = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bp[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_om = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_om[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_op = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_op[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_mp = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_mp[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_b = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_b[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_o = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_o[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_m = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_m[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_p = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_p[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)
pye_cod_none = sum((log(cod_catch[cod_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_none[(t_tune-length(prj_ts)+1):t_tune, "cod"]))^2)

pye_flounder_bomp = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bomp[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_bom = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bom[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_bop = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bop[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_bmp = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bmp[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_omp = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_omp[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_bo = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bo[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_bm = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bm[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_bp = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bp[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_om = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_om[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_op = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_op[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_mp = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_mp[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_b = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_b[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_o = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_o[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_m = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_m[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_p = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_p[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)
pye_flounder_none = sum((log(flounder_catch[flounder_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_none[(t_tune-length(prj_ts)+1):t_tune, "flounder"]))^2)

pye_sprat_bomp = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bomp[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_bom = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bom[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_bop = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bop[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_bmp = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bmp[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_omp = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_omp[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_bo = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bo[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_bm = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bm[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_bp = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bp[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_om = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_om[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_op = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_op[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_mp = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_mp[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_b = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_b[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_o = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_o[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_m = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_m[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_p = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_p[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)
pye_sprat_none = sum((log(sprat_catch[sprat_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_none[(t_tune-length(prj_ts)+1):t_tune, "sprat"]))^2)

pye_herring_bomp = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bomp[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_bom = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bom[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_bop = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bop[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_bmp = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bmp[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_omp = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_omp[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_bo = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bo[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_bm = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bm[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_bp = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_bp[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_om = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_om[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_op = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_op[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_mp = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_mp[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_b = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_b[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_o = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_o[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_m = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_m[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_p = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_p[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)
pye_herring_none = sum((log(herring_catch[herring_catch$Year %in% prj_ts, ]$Catch*1000*1000) - log(yield_none[(t_tune-length(prj_ts)+1):t_tune, "herring"]))^2)

# Sum across species
pye_cod = c(pye_cod_bomp, pye_cod_bom, pye_cod_bop, pye_cod_bmp, pye_cod_omp, pye_cod_bo, pye_cod_bm, pye_cod_bp, pye_cod_om, pye_cod_op, pye_cod_mp, pye_cod_b, pye_cod_o, pye_cod_m, pye_cod_p, pye_cod_none)
pye_flounder = c(pye_flounder_bomp, pye_flounder_bom, pye_flounder_bop, pye_flounder_bmp, pye_flounder_omp, pye_flounder_bo, pye_flounder_bm, pye_flounder_bp, pye_flounder_om, pye_flounder_op, pye_flounder_mp, pye_flounder_b, pye_flounder_o, pye_flounder_m, pye_flounder_p, pye_flounder_none)
pye_sprat = c(pye_sprat_bomp, pye_sprat_bom, pye_sprat_bop, pye_sprat_bmp, pye_sprat_omp, pye_sprat_bo, pye_sprat_bm, pye_sprat_bp, pye_sprat_om, pye_sprat_op, pye_sprat_mp, pye_sprat_b, pye_sprat_o, pye_sprat_m, pye_sprat_p, pye_sprat_none)
pye_herring = c(pye_herring_bomp, pye_herring_bom, pye_herring_bop, pye_herring_bmp, pye_herring_omp, pye_herring_bo, pye_herring_bm, pye_herring_bp, pye_herring_om, pye_herring_op, pye_herring_mp, pye_herring_b, pye_herring_o, pye_herring_m, pye_herring_p, pye_herring_none)

# Total and ranked error
pye_cod_weight = 0.25
pye_flounder_weight = 0.25 
pye_sprat_weight = 0.25
pye_herring_weight = 0.25
pye_error = pye_cod_weight*pye_cod + pye_flounder_weight*pye_flounder + pye_sprat_weight*pye_sprat + pye_herring_weight*pye_herring
pye_rank = (pye_error - min(pye_error)) / (max(pye_error) - min(pye_error))

# Growth projection error
pge_cod_bomp = proj_growth(params_bomp, prj_ts, 100, "cod")
pge_cod_bom = proj_growth(params_bom, prj_ts, 100, "cod")
pge_cod_bop = proj_growth(params_bop, prj_ts, 100, "cod")
pge_cod_bmp = proj_growth(params_bmp, prj_ts, 100, "cod")
pge_cod_omp = proj_growth(params_omp, prj_ts, 100, "cod")
pge_cod_bo = proj_growth(params_bo, prj_ts, 100, "cod")
pge_cod_bm = proj_growth(params_bm, prj_ts, 100, "cod")
pge_cod_bp = proj_growth(params_bp, prj_ts, 100, "cod")
pge_cod_om = proj_growth(params_om, prj_ts, 100, "cod")
pge_cod_op = proj_growth(params_op, prj_ts, 100, "cod")
pge_cod_mp = proj_growth(params_mp, prj_ts, 100, "cod")
pge_cod_b = proj_growth(params_b, prj_ts, 100, "cod")
pge_cod_o = proj_growth(params_o, prj_ts, 100, "cod")
pge_cod_m = proj_growth(params_m, prj_ts, 100, "cod")
pge_cod_p = proj_growth(params_p, prj_ts, 100, "cod")
pge_cod_none = proj_growth(params_none, prj_ts, 100, "cod")

pge_flounder_bomp = proj_growth(params_bomp, prj_ts, 100, "flounder")
pge_flounder_bom = proj_growth(params_bom, prj_ts, 100, "flounder")
pge_flounder_bop = proj_growth(params_bop, prj_ts, 100, "flounder")
pge_flounder_bmp = proj_growth(params_bmp, prj_ts, 100, "flounder")
pge_flounder_omp = proj_growth(params_omp, prj_ts, 100, "flounder")
pge_flounder_bo = proj_growth(params_bo, prj_ts, 100, "flounder")
pge_flounder_bm = proj_growth(params_bm, prj_ts, 100, "flounder")
pge_flounder_bp = proj_growth(params_bp, prj_ts, 100, "flounder")
pge_flounder_om = proj_growth(params_om, prj_ts, 100, "flounder")
pge_flounder_op = proj_growth(params_op, prj_ts, 100, "flounder")
pge_flounder_mp = proj_growth(params_mp, prj_ts, 100, "flounder")
pge_flounder_b = proj_growth(params_b, prj_ts, 100, "flounder")
pge_flounder_o = proj_growth(params_o, prj_ts, 100, "flounder")
pge_flounder_m = proj_growth(params_m, prj_ts, 100, "flounder")
pge_flounder_p = proj_growth(params_p, prj_ts, 100, "flounder")
pge_flounder_none = proj_growth(params_none, prj_ts, 100, "flounder")

pge_sprat_bomp = proj_growth(params_bomp, prj_ts, 100, "sprat")
pge_sprat_bom = proj_growth(params_bom, prj_ts, 100, "sprat")
pge_sprat_bop = proj_growth(params_bop, prj_ts, 100, "sprat")
pge_sprat_bmp = proj_growth(params_bmp, prj_ts, 100, "sprat")
pge_sprat_omp = proj_growth(params_omp, prj_ts, 100, "sprat")
pge_sprat_bo = proj_growth(params_bo, prj_ts, 100, "sprat")
pge_sprat_bm = proj_growth(params_bm, prj_ts, 100, "sprat")
pge_sprat_bp = proj_growth(params_bp, prj_ts, 100, "sprat")
pge_sprat_om = proj_growth(params_om, prj_ts, 100, "sprat")
pge_sprat_op = proj_growth(params_op, prj_ts, 100, "sprat")
pge_sprat_mp = proj_growth(params_mp, prj_ts, 100, "sprat")
pge_sprat_b = proj_growth(params_b, prj_ts, 100, "sprat")
pge_sprat_o = proj_growth(params_o, prj_ts, 100, "sprat")
pge_sprat_m = proj_growth(params_m, prj_ts, 100, "sprat")
pge_sprat_p = proj_growth(params_p, prj_ts, 100, "sprat")
pge_sprat_none = proj_growth(params_none, prj_ts, 100, "sprat")

pge_herring_bomp = proj_growth(params_bomp, prj_ts, 100, "herring")
pge_herring_bom = proj_growth(params_bom, prj_ts, 100, "herring")
pge_herring_bop = proj_growth(params_bop, prj_ts, 100, "herring")
pge_herring_bmp = proj_growth(params_bmp, prj_ts, 100, "herring")
pge_herring_omp = proj_growth(params_omp, prj_ts, 100, "herring")
pge_herring_bo = proj_growth(params_bo, prj_ts, 100, "herring")
pge_herring_bm = proj_growth(params_bm, prj_ts, 100, "herring")
pge_herring_bp = proj_growth(params_bp, prj_ts, 100, "herring")
pge_herring_om = proj_growth(params_om, prj_ts, 100, "herring")
pge_herring_op = proj_growth(params_op, prj_ts, 100, "herring")
pge_herring_mp = proj_growth(params_mp, prj_ts, 100, "herring")
pge_herring_b = proj_growth(params_b, prj_ts, 100, "herring")
pge_herring_o = proj_growth(params_o, prj_ts, 100, "herring")
pge_herring_m = proj_growth(params_m, prj_ts, 100, "herring")
pge_herring_p = proj_growth(params_p, prj_ts, 100, "herring")
pge_herring_none = proj_growth(params_none, prj_ts, 100, "herring")

# Sum across species
pge_cod = c(pge_cod_bomp, pge_cod_bom, pge_cod_bop, pge_cod_bmp, pge_cod_omp, pge_cod_bo, pge_cod_bm, pge_cod_bp, pge_cod_om, pge_cod_op, pge_cod_mp, pge_cod_b, pge_cod_o, pge_cod_m, pge_cod_p, pge_cod_none)
pge_flounder = c(pge_flounder_bomp, pge_flounder_bom, pge_flounder_bop, pge_flounder_bmp, pge_flounder_omp, pge_flounder_bo, pge_flounder_bm, pge_flounder_bp, pge_flounder_om, pge_flounder_op, pge_flounder_mp, pge_flounder_b, pge_flounder_o, pge_flounder_m, pge_flounder_p, pge_flounder_none)
pge_sprat = c(pge_sprat_bomp, pge_sprat_bom, pge_sprat_bop, pge_sprat_bmp, pge_sprat_omp, pge_sprat_bo, pge_sprat_bm, pge_sprat_bp, pge_sprat_om, pge_sprat_op, pge_sprat_mp, pge_sprat_b, pge_sprat_o, pge_sprat_m, pge_sprat_p, pge_sprat_none)
pge_herring = c(pge_herring_bomp, pge_herring_bom, pge_herring_bop, pge_herring_bmp, pge_herring_omp, pge_herring_bo, pge_herring_bm, pge_herring_bp, pge_herring_om, pge_herring_op, pge_herring_mp, pge_herring_b, pge_herring_o, pge_herring_m, pge_herring_p, pge_herring_none)

# Total and ranked error
pge_cod_weight = 0.25
pge_flounder_weight = 0.25 
pge_sprat_weight = 0.25
pge_herring_weight = 0.25
pge_error = pge_cod_weight*pge_cod + pge_flounder_weight*pge_flounder + pge_sprat_weight*pge_sprat + pge_herring_weight*pge_herring
pge_rank = (pge_error - min(pge_error)) / (max(pge_error) - min(pge_error))

# Calculate final weighted error
ssb_weight = 0.3
yield_weight = 0.1
growth_weight = 0.6
final_error_sim = ssb_weight*pbe_rank + yield_weight*pye_rank + growth_weight*pge_rank
model_labels[order(final_error_sim)]

# Create and save error table
proj_error = data.frame(model = model_labels)
proj_error$SSB = ssb_weight*pbe_rank
proj_error$yield = yield_weight*pye_rank
proj_error$growth = growth_weight*pge_rank
proj_error$rank = proj_error$SSB + proj_error$yield + proj_error$growth
proj_error$SSB_error = pbe_error
proj_error$yield_error = pye_error
proj_error$growth_error = pge_error
write.table(proj_error, "Calibration/proj_error.txt")