# Load and plot calibration results

# Created: August 29, 2022
# Last modified: April 1, 2023

# Set working directory
setwd("/Users/epdus/OneDrive/Breathless/Code/Packets/physoxiaMizer")

# Load packages
library(mizer)
library(RColorBrewer)
library(png)
library(xtable)
library(plotrix)
library(ggplot2)
library(ggpubr)

# Load common data
load("Data/oxy_model.rda")

# Contents (ctrl-f):
#	I. Source scripts
#	II. Load models and results with metabolic prioritization
#	III. Load models and results without metabolic prioritization
#	IV. Prepare and run simulations
#	V. Parameter estimates table
#	VI. Plot sensitivity analysis projection error
#	VII. Plot projection error comparison
# VIII. Plot energy flow
# IX. Plot scaling
# X. Plot comparative yield
# XI. Plot comparative biomass


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
#	rate_constantphys: rate scaling functions for constant physiological scaling
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
source("Code/rate_constantphys.R")
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
source("Code/scenarios.R")


########## II. Load models and results with metabolic prioritization ##########

# Load calibration and projection results
phys_cal_error = read.table("Calibration/cal_error.txt", header = T)
phys_cal_one = read.table("Sensitivity/Output/cal_one.txt", header = T)
phys_cal_two = read.table("Sensitivity/Output/cal_two.txt", header = T)
phys_cal_three = read.table("Sensitivity/Output/cal_three.txt", header = T)

# Load calibration and projection results
phys_proj_error = read.table("Calibration/proj_error.txt", header = T)
phys_proj_one = read.table("Sensitivity/Output/proj_one.txt", header = T)
phys_proj_two = read.table("Sensitivity/Output/proj_two.txt", header = T)
phys_proj_three = read.table("Sensitivity/Output/proj_three.txt", header = T)

# Benthos, occupancy, mortality, and physiological scaling
load("Calibration/params_full_bomp.rda")
params_bomp = params_full
load("Sensitivity/One/Calibration/params_full_bomp.rda")
params_bomp_one = params_full
load("Sensitivity/Two/Calibration/params_full_bomp.rda")
params_bomp_two = params_full
load("Sensitivity/Three/Calibration/params_full_bomp.rda")
params_bomp_three = params_full

# Benthos, occupancy, and mortality scaling
load("Calibration/params_full_bom.rda")
params_bom = params_full
load("Sensitivity/One/Calibration/params_full_bom.rda")
params_bom_one = params_full
load("Sensitivity/Two/Calibration/params_full_bom.rda")
params_bom_two = params_full
load("Sensitivity/Three/Calibration/params_full_bom.rda")
params_bom_three = params_full

# Benthos, occupancy, and physiological scaling
load("Calibration/params_full_bop.rda")
params_bop = params_full
load("Sensitivity/One/Calibration/params_full_bop.rda")
params_bop_one = params_full
load("Sensitivity/Two/Calibration/params_full_bop.rda")
params_bop_two = params_full
load("Sensitivity/Three/Calibration/params_full_bop.rda")
params_bop_three = params_full

# Benthos, mortality, and physiological scaling
load("Calibration/params_full_bmp.rda")
params_bmp = params_full
load("Sensitivity/One/Calibration/params_full_bmp.rda")
params_bmp_one = params_full
load("Sensitivity/Two/Calibration/params_full_bmp.rda")
params_bmp_two = params_full
load("Sensitivity/Three/Calibration/params_full_bmp.rda")
params_bmp_three = params_full

# Occupancy, mortality, and physiological scaling
load("Calibration/params_full_omp.rda")
params_omp = params_full
load("Sensitivity/One/Calibration/params_full_omp.rda")
params_omp_one = params_full
load("Sensitivity/Two/Calibration/params_full_omp.rda")
params_omp_two = params_full
load("Sensitivity/Three/Calibration/params_full_omp.rda")
params_omp_three = params_full

# Benthos and occupancy scaling
load("Calibration/params_full_bo.rda")
params_bo = params_full
load("Sensitivity/One/Calibration/params_full_bo.rda")
params_bo_one = params_full
load("Sensitivity/Two/Calibration/params_full_bo.rda")
params_bo_two = params_full
load("Sensitivity/Three/Calibration/params_full_bo.rda")
params_bo_three = params_full

# Benthos and mortality scaling
load("Calibration/params_full_bm.rda")
params_bm = params_full
load("Sensitivity/One/Calibration/params_full_bm.rda")
params_bm_one = params_full
load("Sensitivity/Two/Calibration/params_full_bm.rda")
params_bm_two = params_full
load("Sensitivity/Three/Calibration/params_full_bm.rda")
params_bm_three = params_full

# Benthos and physiological scaling
load("Calibration/params_full_bp.rda")
params_bp = params_full
load("Sensitivity/One/Calibration/params_full_bp.rda")
params_bp_one = params_full
load("Sensitivity/Two/Calibration/params_full_bp.rda")
params_bp_two = params_full
load("Sensitivity/Three/Calibration/params_full_bp.rda")
params_bp_three = params_full

# Occupancy and mortality scaling
load("Calibration/params_full_om.rda")
params_om = params_full
load("Sensitivity/One/Calibration/params_full_om.rda")
params_om_one = params_full
load("Sensitivity/Two/Calibration/params_full_om.rda")
params_om_two = params_full
load("Sensitivity/Three/Calibration/params_full_om.rda")
params_om_three = params_full

# Occupancy and physiological scaling
load("Calibration/params_full_op.rda")
params_op = params_full
load("Sensitivity/One/Calibration/params_full_op.rda")
params_op_one = params_full
load("Sensitivity/Two/Calibration/params_full_op.rda")
params_op_two = params_full
load("Sensitivity/Three/Calibration/params_full_op.rda")
params_op_three = params_full

# Mortality and physiological scaling
load("Calibration/params_full_mp.rda")
params_mp = params_full
load("Sensitivity/One/Calibration/params_full_mp.rda")
params_mp_one = params_full
load("Sensitivity/Two/Calibration/params_full_mp.rda")
params_mp_two = params_full
load("Sensitivity/Three/Calibration/params_full_mp.rda")
params_mp_three = params_full

# Benthos scaling
load("Calibration/params_full_b.rda")
params_b = params_full
load("Sensitivity/One/Calibration/params_full_b.rda")
params_b_one = params_full
load("Sensitivity/Two/Calibration/params_full_b.rda")
params_b_two = params_full
load("Sensitivity/Three/Calibration/params_full_b.rda")
params_b_three = params_full

# Occupancy scaling
load("Calibration/params_full_o.rda")
params_o = params_full
load("Sensitivity/One/Calibration/params_full_o.rda")
params_o_one = params_full
load("Sensitivity/Two/Calibration/params_full_o.rda")
params_o_two = params_full
load("Sensitivity/Three/Calibration/params_full_o.rda")
params_o_three = params_full

# Mortality scaling
load("Calibration/params_full_m.rda")
params_m = params_full
load("Sensitivity/One/Calibration/params_full_m.rda")
params_m_one = params_full
load("Sensitivity/Two/Calibration/params_full_m.rda")
params_m_two = params_full
load("Sensitivity/Three/Calibration/params_full_m.rda")
params_m_three = params_full

# Physiological scaling
load("Calibration/params_full_p.rda")
params_p = params_full
load("Sensitivity/One/Calibration/params_full_p.rda")
params_p_one = params_full
load("Sensitivity/Two/Calibration/params_full_p.rda")
params_p_two = params_full
load("Sensitivity/Three/Calibration/params_full_p.rda")
params_p_three = params_full

# None scaling (with left fish)
load("Calibration/params_full_none.rda")
params_none = params_full
load("Sensitivity/One/Calibration/params_full_none.rda")
params_none_one = params_full
load("Sensitivity/Two/Calibration/params_full_none.rda")
params_none_two = params_full
load("Sensitivity/Three/Calibration/params_full_none.rda")
params_none_three = params_full


########## III. Load models and results without metabolic prioritization ##########

# Load calibration and projection results
hyp_cal_error = read.table("Constant Phys/cal_error.txt", header = T)
hyp_cal_one = read.table("Constant Phys/Sensitivity/Output/cal_one.txt", header = T)
hyp_cal_two = read.table("Constant Phys/Sensitivity/Output/cal_two.txt", header = T)
hyp_cal_three = read.table("Constant Phys/Sensitivity/Output/cal_three.txt", header = T)

# Load calibration and projection results
hyp_proj_error = read.table("Constant Phys/proj_error.txt", header = T)
hyp_proj_one = read.table("Constant Phys/Sensitivity/Output/proj_one.txt", header = T)
hyp_proj_two = read.table("Constant Phys/Sensitivity/Output/proj_two.txt", header = T)
hyp_proj_three = read.table("Constant Phys/Sensitivity/Output/proj_three.txt", header = T)

# Benthos, occupancy, mortality, and physiological scaling without metabolic prioritization
load("Constant Phys/params_full_bomp.rda")
cphys_bomp = params_full
load("Constant Phys/Sensitivity/One/Calibration/params_full_bomp.rda")
cphys_bomp_one = params_full
load("Constant Phys/Sensitivity/Two/Calibration/params_full_bomp.rda")
cphys_bomp_two = params_full
load("Constant Phys/Sensitivity/Three/Calibration/params_full_bomp.rda")
cphys_bomp_three = params_full

# Benthos, occupancy, and mortality scaling without metabolic prioritization
load("Constant Phys/params_full_bom.rda")
cphys_bom = params_full
load("Constant Phys/Sensitivity/One/Calibration/params_full_bom.rda")
cphys_bom_one = params_full
load("Constant Phys/Sensitivity/Two/Calibration/params_full_bom.rda")
cphys_bom_two = params_full
load("Constant Phys/Sensitivity/Three/Calibration/params_full_bom.rda")
cphys_bom_three = params_full

# Benthos, occupancy, and physiological scaling without metabolic prioritization
load("Constant Phys/params_full_bop.rda")
cphys_bop = params_full
load("Constant Phys/Sensitivity/One/Calibration/params_full_bop.rda")
cphys_bop_one = params_full
load("Constant Phys/Sensitivity/Two/Calibration/params_full_bop.rda")
cphys_bop_two = params_full
load("Constant Phys/Sensitivity/Three/Calibration/params_full_bop.rda")
cphys_bop_three = params_full

# Benthos, mortality, and physiological scaling without metabolic prioritization
load("Constant Phys/params_full_bmp.rda")
cphys_bmp = params_full
load("Constant Phys/Sensitivity/One/Calibration/params_full_bmp.rda")
cphys_bmp_one = params_full
load("Constant Phys/Sensitivity/Two/Calibration/params_full_bmp.rda")
cphys_bmp_two = params_full
load("Constant Phys/Sensitivity/Three/Calibration/params_full_bmp.rda")
cphys_bmp_three = params_full

# Occupancy, mortality, and physiological scaling without metabolic prioritization
load("Constant Phys/params_full_omp.rda")
cphys_omp = params_full
load("Constant Phys/Sensitivity/One/Calibration/params_full_omp.rda")
cphys_omp_one = params_full
load("Constant Phys/Sensitivity/Two/Calibration/params_full_omp.rda")
cphys_omp_two = params_full
load("Constant Phys/Sensitivity/Three/Calibration/params_full_omp.rda")
cphys_omp_three = params_full

# Benthos and occupancy scaling without metabolic prioritization
load("Constant Phys/params_full_bo.rda")
cphys_bo = params_full
load("Constant Phys/Sensitivity/One/Calibration/params_full_bo.rda")
cphys_bo_one = params_full
load("Constant Phys/Sensitivity/Two/Calibration/params_full_bo.rda")
cphys_bo_two = params_full
load("Constant Phys/Sensitivity/Three/Calibration/params_full_bo.rda")
cphys_bo_three = params_full

# Benthos and mortality scaling without metabolic prioritization
load("Constant Phys/params_full_bm.rda")
cphys_bm = params_full
load("Constant Phys/Sensitivity/One/Calibration/params_full_bm.rda")
cphys_bm_one = params_full
load("Constant Phys/Sensitivity/Two/Calibration/params_full_bm.rda")
cphys_bm_two = params_full
load("Constant Phys/Sensitivity/Three/Calibration/params_full_bm.rda")
cphys_bm_three = params_full

# Benthos and physiological scaling without metabolic prioritization
load("Constant Phys/params_full_bp.rda")
cphys_bp = params_full
load("Constant Phys/Sensitivity/One/Calibration/params_full_bp.rda")
cphys_bp_one = params_full
load("Constant Phys/Sensitivity/Two/Calibration/params_full_bp.rda")
cphys_bp_two = params_full
load("Constant Phys/Sensitivity/Three/Calibration/params_full_bp.rda")
cphys_bp_three = params_full

# Occupancy and mortality scaling without metabolic prioritization
load("Constant Phys/params_full_om.rda")
cphys_om = params_full
load("Constant Phys/Sensitivity/One/Calibration/params_full_om.rda")
cphys_om_one = params_full
load("Constant Phys/Sensitivity/Two/Calibration/params_full_om.rda")
cphys_om_two = params_full
load("Constant Phys/Sensitivity/Three/Calibration/params_full_om.rda")
cphys_om_three = params_full

# Occupancy and physiological scaling without metabolic prioritization
load("Constant Phys/params_full_op.rda")
cphys_op = params_full
load("Constant Phys/Sensitivity/One/Calibration/params_full_op.rda")
cphys_op_one = params_full
load("Constant Phys/Sensitivity/Two/Calibration/params_full_op.rda")
cphys_op_two = params_full
load("Constant Phys/Sensitivity/Three/Calibration/params_full_op.rda")
cphys_op_three = params_full

# Mortality and physiological scaling without metabolic prioritization
load("Constant Phys/params_full_mp.rda")
cphys_mp = params_full
load("Constant Phys/Sensitivity/One/Calibration/params_full_mp.rda")
cphys_mp_one = params_full
load("Constant Phys/Sensitivity/Two/Calibration/params_full_mp.rda")
cphys_mp_two = params_full
load("Constant Phys/Sensitivity/Three/Calibration/params_full_mp.rda")
cphys_mp_three = params_full

# Benthos scaling without metabolic prioritization
load("Constant Phys/params_full_b.rda")
cphys_b = params_full
load("Constant Phys/Sensitivity/One/Calibration/params_full_b.rda")
cphys_b_one = params_full
load("Constant Phys/Sensitivity/Two/Calibration/params_full_b.rda")
cphys_b_two = params_full
load("Constant Phys/Sensitivity/Three/Calibration/params_full_b.rda")
cphys_b_three = params_full

# Occupancy scaling without metabolic prioritization
load("Constant Phys/params_full_o.rda")
cphys_o = params_full
load("Constant Phys/Sensitivity/One/Calibration/params_full_o.rda")
cphys_o_one = params_full
load("Constant Phys/Sensitivity/Two/Calibration/params_full_o.rda")
cphys_o_two = params_full
load("Constant Phys/Sensitivity/Three/Calibration/params_full_o.rda")
cphys_o_three = params_full

# Mortality scaling without metabolic prioritization
load("Constant Phys/params_full_m.rda")
cphys_m = params_full
load("Constant Phys/Sensitivity/One/Calibration/params_full_m.rda")
cphys_m_one = params_full
load("Constant Phys/Sensitivity/Two/Calibration/params_full_m.rda")
cphys_m_two = params_full
load("Constant Phys/Sensitivity/Three/Calibration/params_full_m.rda")
cphys_m_three = params_full

# Physiological scaling without metabolic prioritization
load("Constant Phys/params_full_p.rda")
cphys_p = params_full
load("Constant Phys/Sensitivity/One/Calibration/params_full_p.rda")
cphys_p_one = params_full
load("Constant Phys/Sensitivity/Two/Calibration/params_full_p.rda")
cphys_p_two = params_full
load("Constant Phys/Sensitivity/Three/Calibration/params_full_p.rda")
cphys_p_three = params_full

# None scaling (with left fish) without metabolic prioritization
load("Constant Phys/params_full_none.rda")
cphys_none = params_full
load("Constant Phys/Sensitivity/One/Calibration/params_full_none.rda")
cphys_none_one = params_full
load("Constant Phys/Sensitivity/Two/Calibration/params_full_none.rda")
cphys_none_two = params_full
load("Constant Phys/Sensitivity/Three/Calibration/params_full_none.rda")
cphys_none_three = params_full


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

cphys_bomp_sim = cphys_bomp
cphys_bomp_sim@other_params$benthic_oxygen = benthic_oxygen_sim
cphys_bomp_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_bom_sim = params_bom
params_bom_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_bom_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

cphys_bom_sim = cphys_bom
cphys_bom_sim@other_params$benthic_oxygen = benthic_oxygen_sim
cphys_bom_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_bop_sim = params_bop
params_bop_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_bop_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

cphys_bop_sim = cphys_bop
cphys_bop_sim@other_params$benthic_oxygen = benthic_oxygen_sim
cphys_bop_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_bmp_sim = params_bmp
params_bmp_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_bmp_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

cphys_bmp_sim = cphys_bmp
cphys_bmp_sim@other_params$benthic_oxygen = benthic_oxygen_sim
cphys_bmp_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_omp_sim = params_omp
params_omp_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_omp_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

cphys_omp_sim = cphys_omp
cphys_omp_sim@other_params$benthic_oxygen = benthic_oxygen_sim
cphys_omp_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_bo_sim = params_bo
params_bo_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_bo_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

cphys_bo_sim = cphys_bo
cphys_bo_sim@other_params$benthic_oxygen = benthic_oxygen_sim
cphys_bo_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_bm_sim = params_bm
params_bm_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_bm_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

cphys_bm_sim = cphys_bm
cphys_bm_sim@other_params$benthic_oxygen = benthic_oxygen_sim
cphys_bm_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_bp_sim = params_bp
params_bp_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_bp_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

cphys_bp_sim = cphys_bp
cphys_bp_sim@other_params$benthic_oxygen = benthic_oxygen_sim
cphys_bp_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_om_sim = params_om
params_om_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_om_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

cphys_om_sim = cphys_om
cphys_om_sim@other_params$benthic_oxygen = benthic_oxygen_sim
cphys_om_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_op_sim = params_op
params_op_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_op_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

cphys_op_sim = cphys_op
cphys_op_sim@other_params$benthic_oxygen = benthic_oxygen_sim
cphys_op_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_mp_sim = params_mp
params_mp_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_mp_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

cphys_mp_sim = cphys_mp
cphys_mp_sim@other_params$benthic_oxygen = benthic_oxygen_sim
cphys_mp_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_b_sim = params_b
params_b_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_b_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

cphys_b_sim = cphys_b
cphys_b_sim@other_params$benthic_oxygen = benthic_oxygen_sim
cphys_b_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_o_sim = params_o
params_o_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_o_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

cphys_o_sim = cphys_o
cphys_o_sim@other_params$benthic_oxygen = benthic_oxygen_sim
cphys_o_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_m_sim = params_m
params_m_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_m_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

cphys_m_sim = cphys_m
cphys_m_sim@other_params$benthic_oxygen = benthic_oxygen_sim
cphys_m_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_p_sim = params_p
params_p_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_p_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

cphys_p_sim = cphys_p
cphys_p_sim@other_params$benthic_oxygen = benthic_oxygen_sim
cphys_p_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

params_none_sim = params_none
params_none_sim@other_params$benthic_oxygen = benthic_oxygen_sim
params_none_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

cphys_none_sim = cphys_none
cphys_none_sim@other_params$benthic_oxygen = benthic_oxygen_sim
cphys_none_sim@other_params$pelagic_oxygen = pelagic_oxygen_sim

# Scale rates by oxygen and temperature: see otmscale.R
params_bomp_sim = occupancy(params_bomp_sim, t_tune)
params_bomp_sim = rate_scale(params_bomp_sim, t_tune)
params_bomp_sim = mscale(params_bomp_sim, t_tune)

cphys_bomp_sim = occupancy(cphys_bomp_sim, t_tune)
cphys_bomp_sim = rs_cphys(cphys_bomp_sim, t_tune)
cphys_bomp_sim = mscale(cphys_bomp_sim, t_tune)

params_bom_sim = occupancy(params_bom_sim, t_tune)

cphys_bom_sim = occupancy(cphys_bom_sim, t_tune)

params_bop_sim = occupancy(params_bop_sim, t_tune)
params_bop_sim = rate_scale(params_bop_sim, t_tune)
params_bop_sim = mscale(params_bop_sim, t_tune)

cphys_bop_sim = occupancy(cphys_bop_sim, t_tune)
cphys_bop_sim = rs_cphys(cphys_bop_sim, t_tune)
cphys_bop_sim = mscale(cphys_bop_sim, t_tune)

params_bmp_sim = rate_scale(params_bmp_sim, t_tune)
params_bmp_sim = mscale(params_bmp_sim, t_tune)

cphys_bmp_sim = rs_cphys(cphys_bmp_sim, t_tune)
cphys_bmp_sim = mscale(cphys_bmp_sim, t_tune)

params_omp_sim = occupancy(params_omp_sim, t_tune)
params_omp_sim = rate_scale(params_omp_sim, t_tune)
params_omp_sim = mscale(params_omp_sim, t_tune)

cphys_omp_sim = occupancy(cphys_omp_sim, t_tune)
cphys_omp_sim = rs_cphys(cphys_omp_sim, t_tune)
cphys_omp_sim = mscale(cphys_omp_sim, t_tune)

params_bo_sim = occupancy(params_bo_sim, t_tune)

cphys_bo_sim = occupancy(cphys_bo_sim, t_tune)

params_bp_sim = rate_scale(params_bp_sim, t_tune)
params_bp_sim = mscale(params_bp_sim, t_tune)

cphys_bp_sim = rs_cphys(cphys_bp_sim, t_tune)
cphys_bp_sim = mscale(cphys_bp_sim, t_tune)

params_om_sim = occupancy(params_om_sim, t_tune)

cphys_om_sim = occupancy(cphys_om_sim, t_tune)

params_op_sim = occupancy(params_op_sim, t_tune)
params_op_sim = rate_scale(params_op_sim, t_tune)
params_op_sim = mscale(params_op_sim, t_tune)

cphys_op_sim = occupancy(cphys_op_sim, t_tune)
cphys_op_sim = rs_cphys(cphys_op_sim, t_tune)
cphys_op_sim = mscale(cphys_op_sim, t_tune)

params_mp_sim = rate_scale(params_mp_sim, t_tune)
params_mp_sim = mscale(params_mp_sim, t_tune)

cphys_mp_sim = rs_cphys(cphys_mp_sim, t_tune)
cphys_mp_sim = mscale(cphys_mp_sim, t_tune)

params_o_sim = occupancy(params_o_sim, t_tune)

cphys_o_sim = occupancy(cphys_o_sim, t_tune)

params_p_sim = rate_scale(params_p_sim, t_tune)
params_p_sim = mscale(params_p_sim, t_tune)

cphys_p_sim = rs_cphys(cphys_p_sim, t_tune)
cphys_p_sim = mscale(cphys_p_sim, t_tune)

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

simcp_bomp = project(cphys_bomp_sim, t_max = 100, effort = effort_full)
simcp_bom = project(cphys_bom_sim, t_max = 100, effort = effort_full)
simcp_bop = project(cphys_bop_sim, t_max = 100, effort = effort_full)
simcp_bmp = project(cphys_bmp_sim, t_max = 100, effort = effort_full)
simcp_omp = project(cphys_omp_sim, t_max = 100, effort = effort_full)
simcp_bo = project(cphys_bo_sim, t_max = 100, effort = effort_full)
simcp_bm = project(cphys_bm_sim, t_max = 100, effort = effort_full)
simcp_bp = project(cphys_bp_sim, t_max = 100, effort = effort_full)
simcp_om = project(cphys_om_sim, t_max = 100, effort = effort_full)
simcp_op = project(cphys_op_sim, t_max = 100, effort = effort_full)
simcp_mp = project(cphys_mp_sim, t_max = 100, effort = effort_full)
simcp_b = project(cphys_b_sim, t_max = 100, effort = effort_full)
simcp_o = project(cphys_o_sim, t_max = 100, effort = effort_full)
simcp_m = project(cphys_m_sim, t_max = 100, effort = effort_full)
simcp_p = project(cphys_p_sim, t_max = 100, effort = effort_full)
simcp_none = project(cphys_none_sim, t_max = 100, effort = effort_full)


########## V. Parameter estimates table ##########

# Create matrix to house values
table_ests = matrix(nrow = 16, ncol = 12)

# Column names
colnames(table_ests) = c("Species", "Model", "Us", "as", "Ui", "ai", "Ua", "aa", "Ue", "ae", "Um", "am")

# Species
table_ests[,"Species"] = rep(c("Cod","Flounder"), each = 8)

# Models
table_ests[,"Model"] = c("BOMP", "BOP", "BMP", "OMP", "BP", "OP", "MP", "P")

# Create lists of mizer objects with metabolic prioritization
params = list(params_bomp, 
				params_bop, 
				params_bmp, 
				params_omp, 
				params_bp, 
				params_op, 
				params_mp, 
				params_p)
params_one = list(params_bomp_one, 
				params_bop_one, 
				params_bmp_one, 
				params_omp_one, 
				params_bp_one, 
				params_op_one, 
				params_mp_one, 
				params_p_one)
params_two = list(params_bomp_two, 
				params_bop_two, 
				params_bmp_two, 
				params_omp_two, 
				params_bp_two, 
				params_op_two, 
				params_mp_two, 
				params_p_two)
params_three = list(params_bomp_three, 
				params_bop_three, 
				params_bmp_three, 
				params_omp_three, 
				params_bp_three, 
				params_op_three, 
				params_mp_three, 
				params_p_three)
names(params) = c("BOMP", "BOP", "BMP", "OMP", "BP", "OP", "MP", "P")
names(params_one) = c("BOMP", "BOP", "BMP", "OMP", "BP", "OP", "MP", "P")
names(params_two) = c("BOMP", "BOP", "BMP", "OMP", "BP", "OP", "MP", "P")
names(params_three) = c("BOMP", "BOP", "BMP", "OMP", "BP", "OP", "MP", "P")

# Create lists of mizer objects with metabolic prioritization
cphys = list(cphys_bomp, 
				cphys_bop, 
				cphys_bmp, 
				cphys_omp, 
				cphys_bp, 
				cphys_op, 
				cphys_mp, 
				cphys_p)
cphys_one = list(cphys_bomp_one, 
				cphys_bop_one, 
				cphys_bmp_one, 
				cphys_omp_one, 
				cphys_bp_one, 
				cphys_op_one, 
				cphys_mp_one, 
				cphys_p_one)
cphys_two = list(cphys_bomp_two, 
				cphys_bop_two, 
				cphys_bmp_two, 
				cphys_omp_two, 
				cphys_bp_two, 
				cphys_op_two, 
				cphys_mp_two, 
				cphys_p_two)
cphys_three = list(cphys_bomp_three, 
				cphys_bop_three, 
				cphys_bmp_three, 
				cphys_omp_three, 
				cphys_bp_three, 
				cphys_op_three, 
				cphys_mp_three, 
				cphys_p_three)
names(cphys) = c("BOMP", "BOP", "BMP", "OMP", "BP", "OP", "MP", "P")
names(cphys_one) = c("BOMP", "BOP", "BMP", "OMP", "BP", "OP", "MP", "P")
names(cphys_two) = c("BOMP", "BOP", "BMP", "OMP", "BP", "OP", "MP", "P")
names(cphys_three) = c("BOMP", "BOP", "BMP", "OMP", "BP", "OP", "MP", "P")

# Extraction functions
get_rate = function(mod, rate) {
	
	# Get cod rates
	rng_cod = c(params[[mod]]@species_params["cod",rate],
				params_one[[mod]]@species_params["cod",rate],
				params_two[[mod]]@species_params["cod",rate],
				params_three[[mod]]@species_params["cod",rate])
				
	# Get flounder rates
	rng_fld = c(params[[mod]]@species_params["flounder",rate],
				params_one[[mod]]@species_params["flounder",rate],
				params_two[[mod]]@species_params["flounder",rate],
				params_three[[mod]]@species_params["flounder",rate])
	
	# Cod and flounder for model "mod" and rate "rate"			
	ret = paste(format(round(rng_cod[1],1),nsmall=1), " (", format(round(min(rng_cod),1),nsmall=1), ",", format(round(max(rng_cod),1),nsmall=1), ")", sep = "")
	ret[2] = paste(format(round(rng_fld[1],1),nsmall=1), " (", format(round(min(rng_fld),1),nsmall=1), ",", format(round(max(rng_fld),1),nsmall=1), ")", sep = "")
	
	return(ret)
}

# U and a values for search
table_ests[,"Us"] = c(t(sapply(names(params), get_rate, rate = "U_search")))
table_ests[,"as"] = c(t(sapply(names(params), get_rate, rate = "a_search")))

# U and a values for maximum consumption
table_ests[,"Ui"] = c(t(sapply(names(params), get_rate, rate = "U_maxin")))
table_ests[,"ai"] = c(t(sapply(names(params), get_rate, rate = "a_maxin")))

# U and a values for search
table_ests[,"Ua"] = c(t(sapply(names(params), get_rate, rate = "U_alpha")))
table_ests[,"aa"] = c(t(sapply(names(params), get_rate, rate = "a_alpha")))

# U and a values for search
table_ests[,"Ue"] = c(t(sapply(names(params), get_rate, rate = "U_erepro")))
table_ests[,"ae"] = c(t(sapply(names(params), get_rate, rate = "a_erepro")))

# U and a values for search
table_ests[,"Um"] = c(t(sapply(names(params), get_rate, rate = "U_met")))
table_ests[,"am"] = c(t(sapply(names(params), get_rate, rate = "a_met")))

# Convert to latex table
print(xtable(table_ests), include.rownames = F)


########## VI. Plot sensitivity analysis projection error ##########

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

# Plot colors
error_colors = hcl.colors(n = 3, palette = "Teal")

# Initialize plot
jpeg("Plots/sens_error_sim.jpg", width = 22, height = 9, units = 'cm', res = 600)
par(mfrow = c(1,3), oma = c(3,3,0,0), mar = c(5,2,2,2))

# Plot calibration for first repetition
rep1 = barplot(cbind(growth,yield,SSB) ~ model, data = phys_proj_one[order(phys_proj_one$rank),], ylim = c(0,1), xlab = "", names.arg = rep(NA,16), col = error_colors, cex.axis = 1.2)
axis(1, at = rep1, labels = model_labels[order (phys_proj_one$rank)], las = 2)
box()

# Y-axis label
mtext("Weighted error", 2, cex = 1.2, line = 3)

# Plot calibration for second repetition
rep2 = barplot(cbind(growth,yield,SSB) ~ model, data = phys_proj_two[order(phys_proj_two$rank),], ylim = c(0,1), xlab = "", names.arg = rep(NA,16), col = error_colors, cex.axis = 1.2)
axis(1, at = rep2, labels = model_labels[order (phys_proj_two$rank)], las = 2)
box()

# Plot calibration for third repetition
rep3 = barplot(cbind(growth,yield,SSB) ~ model, data = phys_proj_three[order(phys_proj_three$rank),], ylim = c(0,1), xlab = "", names.arg = rep(NA,16), col = error_colors, cex.axis = 1.2)
axis(1, at = rep3, labels = model_labels[order (phys_proj_three$rank)], las = 2)
box()

# Legend
par(fig = c(0,1,0,1), oma = c(0,1,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("bottom", bty = 'n', c("Growth", "Yield", "SSB"), col = error_colors, pch = 15, pt.cex = 2, xpd = TRUE, horiz = TRUE, cex = 1.2, seg.len=1)

# Finish plot
dev.off()


########## VII. Plot projection error comparison ##########

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

# Get error ranges for models with metabolic prioritization
phys_proj_rng = data.frame(model = model_labels, 
							rank1 = phys_proj_error$rank, 
							rank2 = phys_proj_one$rank, 
							rank3 = phys_proj_two$rank,
							rank4 = phys_proj_three$rank)
phys_proj_rng$min = apply(phys_proj_rng[,2:5],1,min)
phys_proj_rng$max = apply(phys_proj_rng[,2:5],1,max)

# Get mean errors for models without metabolic prioritization
hyp_proj_rng = data.frame(model = model_labels, 
							rank1 = hyp_proj_error$rank, 
							rank2 = hyp_proj_one$rank, 
							rank3 = hyp_proj_two$rank,
							rank4 = hyp_proj_three$rank)
hyp_proj_rng$min = apply(hyp_proj_rng[,2:5],1,min)
hyp_proj_rng$max = apply(hyp_proj_rng[,2:5],1,max)

# Plot colors
col.bp = c(hcl.colors(n = 3, palette = "Teal"), hcl.colors(n = 3, palette = "Peach"))

# Create data frame to house error of both constant and priority scaling
error = data.frame(Scaling = c(rep("Priority",16),rep("Constant",16)))
error$Model = rep(model_labels[order(phys_proj_error$rank)], 2)
error$Category = paste(error$Scaling, error$Model)
error$SSB = c(phys_proj_error$SSB[order(phys_proj_error$rank)], rep(0,16))
error$Yield = c(phys_proj_error$yield[order(phys_proj_error$rank)], rep(0,16))
error$Growth = c(phys_proj_error$growth[order(phys_proj_error$rank)], rep(0,16))
error$SSBcp = c(rep(0,16), hyp_proj_error$SSB[order(phys_proj_error$rank)])
error$Yieldcp = c(rep(0,16), hyp_proj_error$yield[order(phys_proj_error$rank)])
error$Growthcp = c(rep(0,16), hyp_proj_error$growth[order(phys_proj_error$rank)])
error$Min = c(phys_proj_rng$min[order(phys_proj_error$rank)], rep(0,16))
error$Max = c(phys_proj_rng$max[order(phys_proj_error$rank)], rep(0,16))
error$Mincp = c(rep(0,16), hyp_proj_rng$min[order(phys_proj_error$rank)])
error$Maxcp = c(rep(0,16), hyp_proj_rng$max[order(phys_proj_error$rank)])

# Re-order for priority scaling calibration rankings
error = error[rep(seq(16),each=2) + c(0,16), ]
rownames(error) = seq(32)

# Initialize plot
jpeg("Plots/proj_error_compare.jpg", width = 18, height = 14, units = 'cm', res = 600)
par(oma = c(0,0,0,7), mar = c(6,5,3,1))

# Barplot
ebar = barplot(cbind(SSB,Yield,Growth,SSBcp,Yieldcp,Growthcp) ~ Category, data = error, ylim = c(0,1.1), names.arg = rep(NA,32), xlab = "", space = c(0,rep(c(0,1),15),0), col = col.bp, cex.axis = 1.2)
axis(1, at = colMeans(matrix(ebar,nrow=2,ncol=16)), labels = model_labels[order (phys_proj_error$rank)], las = 2)
points(ebar[c(1,3,5,8,9,12,13,15,17,20,22,24,25,27,29,32)], rowSums(error[,4:9])[c(1,3,5,8,9,12,13,15,17,20,22,24,25,27,29,32)]+0.02, pch = 8, cex = 0.5)
box()

# Add SE error bars
arrows(x0 = ebar, y0 = error$Min, x1 = ebar, y1 = error$Max, code = 0)
arrows(x0 = ebar, y0 = error$Mincp, x1 = ebar, y1 = error$Maxcp, code = 0)

# Add outer margin label
mtext("Weighted error", 2, line = 3, cex = 1.2)

# Add legend
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend(0.65, 0.3, c("Growth","Yield","SSB"), bty = "n", pch = 15, pt.cex = 2, col = col.bp[3:1]); text(0.65, 0.33, "MP", adj = 0, font = 2)
legend(0.65, -0.1, c("Growth","Yield","SSB"), bty = "n", pch = 15, pt.cex = 2, col = col.bp[6:4]); text(0.65, -0.07, "No MP", adj = 0, font = 2)
par(mfrow = c(1,1))

# Finish plot
dev.off()


########## VIII. Plot energy flow ##########

# Oxygen values
oxy_benthic = as.list(seq(1, 3, length.out = 20))

# Time
t_max = 100

# All scenarios
phys.scenarios = lapply(oxy_benthic, oxy_scenario, params = params_bop, t = 100, p = T)
hyp.scenarios = lapply(oxy_benthic, oxy_scenario, params = cphys_bom, t = 100, p = F)

# Create data frame to house values
all.cod = data.frame(Oxygen = rep(unlist(oxy_benthic),each=6))
all.cod$Item = c("Benthos","Plankton","Cod","Flounder","Sprat","Herring")
all.cod$Flow_phys = NA
all.cod$Flow_hyp = NA
all.cod$Size_phys = NA
all.cod$Size_hyp = NA

# Change item factor level order
all.cod$Item = factor(all.cod$Item, levels = c("Benthos","Plankton","Cod","Flounder","Sprat","Herring"))

# Benthic oxygen and size
oxy.cod = vector(length = 29)
size.cod = vector(length = 29)
size.phys = vector(length = 29)
size.hyp = vector(length = 29)
for(i in 1:length(oxy.cod)) {
	oxy.cod[i] = params_bop@other_params$benthic_oxygen[i+72]
	size.cod[i] = max(subset(cod_smalk, cod_smalk$Year == i+1990)$IndWgt, na.rm = T) / 1000
	size.phys[i] = oxy_scenario(params_bop, 100, oxy.cod[i], p = T)$size["cod"]
	size.hyp[i] = oxy_scenario(cphys_bom, 100, oxy.cod[i], p = F)$size["cod"]
}
size.df = data.frame(Oxygen = oxy.cod, Size = size.cod, Size_phys = size.phys/1000, Size_hyp = size.hyp/1000)

# Extract size values for cod
for(i in 1:length(phys.scenarios)) {
	ind = seq(1,120,6)[i]
	all.cod$Flow_phys[ind:(ind+5)] = phys.scenarios[[i]]$cod[c(5,6,1,2,3,4)] / 1000
	all.cod$Size_phys[ind] = phys.scenarios[[i]]$size[1] / 1000
	all.cod$Flow_hyp[ind:(ind+5)] = hyp.scenarios[[i]]$cod[c(5,6,1,2,3,4)] / 1000
	all.cod$Size_hyp[ind] = hyp.scenarios[[i]]$size[1] / 1000
}

# Plot colors
col_benthos = rgb(38, 24, 95, maxColorValue=255)
col_plankton = rgb(0, 106, 168, maxColorValue = 255)
col_cod = rgb(0, 166, 174, maxColorValue = 255)
col_flounder = rgb(252, 255, 221, maxColorValue = 255)
col_sprat = rgb(205, 240, 203, maxColorValue = 255)
col_herring = rgb(119, 209, 181, maxColorValue = 255)
col_all = c(col_benthos, col_plankton, col_cod, col_flounder, col_sprat, col_herring)

# Stacked area plot
p1 = ggplot(all.cod, aes(x = Oxygen, y = Flow_phys, fill = Item)) +
		geom_area() + xlab("") + ylab("") + labs(title = "MP", fill = "") +
		scale_fill_manual(values=col_all) +
		geom_line(aes(x = Oxygen, y = 2.5*Size_phys),linewidth=1,color="dodgerblue",linetype="dashed") +
		geom_point(aes(x = Oxygen, y = 2.5*Size), data = size.df, shape = 17, color = "dodgerblue", show.legend = F, inherit.aes = F) +
		scale_x_continuous(limits=c(1,3), breaks = c(1,2,3)) +
		scale_y_continuous(name="",limits=c(0,50),sec.axis=sec_axis(trans=~./2.5, name="")) +
		theme_classic() +
		theme(panel.border=element_rect(colour="black",fill=NA,linewidth=1))

p2 = ggplot(all.cod, aes(x = Oxygen, y = Flow_hyp, fill = Item)) +
		geom_area() + xlab("") + ylab("") + labs(title = "No MP", fill = "") +
		geom_line(aes(x = Oxygen, y = 2.5*Size_hyp),linewidth=1,color="dodgerblue",linetype="dashed") +
		geom_point(aes(x = Oxygen, y = 2.5*Size), data = size.df, shape = 17, color = "dodgerblue", show.legend = F, inherit.aes = F) +
		scale_fill_manual(values=c(col_benthos,col_plankton,col_cod,col_flounder,col_sprat,col_herring)) +
		scale_x_continuous(limits=c(1,3), breaks = c(1,2,3)) +
		scale_y_continuous(name="",limits=c(0,50),sec.axis=sec_axis(trans=~./2.5, name="")) +
		theme_classic() +
		theme(panel.border=element_rect(colour="black",fill=NA,linewidth=1))

p = ggarrange(p1, p2, ncol = 1, common.legend = T, legend = "top", labels = c("a.","b."))

# Initialize plot
jpeg("Plots/cod_flowsize.jpg", width = 9, height = 16, units = 'cm', res = 600)

# Annotate the figure by adding a common labels
annotate_figure(p,
                bottom = text_grob("Oxygen (mL/L)", size = 12),
                left = text_grob("Energy flow (kg/year)", size = 12, rot = 90),
                right = text_grob("Body size (kg)", size = 12, rot = 270))

# Finish plot
dev.off()


########## IX. Plot scaling ##########

# Logistic and exponential scaling
logexp = function(x, mod, species, rate, met = FALSE) {
  U = species_params(mod)[species,paste("U_",rate,sep="")]
  k = species_params(mod)[species,paste("a_",rate,sep="")]
  p = ifelse(species == "cod", cod_pcrit, flounder_pcrit)
  occ = 1 / (1 + exp(-params_bop@species_params[species,]$U_hab*(x-p*params_bop@species_params[species,]$a_hab)))
  ox = 6 * (1-occ) + x*occ
  ret = 1 / (1 + exp(-U*(ox - p*k)))
  if(met) ret = 1 + exp(-U*(ox - p*k))
  return(ret)
}

# Oxygen prediction values
oxy = seq(1, 3, length.out = 100)

# Choose weight for cod and flounder (based on size at maturation)
cod_wt = params_bop@w[79]
flounder_wt = params_bop@w[66]

# Get cod and flounder pcrit
cod_pcrit = params_bop@other_params$P_crit["cod",79]
flounder_pcrit = params_bop@other_params$P_crit["flounder",66]

# Create scaling data frames
cod_scale = data.frame(oxy = oxy)
flounder_scale = data.frame(oxy = oxy)

# Predict cod activity scaling for all BOP U and k values
actmat = matrix(c(logexp(oxy, params_bop, "cod", "search"),
                  logexp(oxy, params_bop_one, "cod", "search"),
                  logexp(oxy, params_bop_two, "cod", "search"),
                  logexp(oxy, params_bop_three, "cod", "search")),
                ncol = 4)
cod_scale$Activity = logexp(oxy, params_bop, "cod", "search")
cod_scale$Activity_min = apply(actmat, 1, min)
cod_scale$Activity_max = apply(actmat, 1, max)

# Predict cod consumption scaling for all BOP U and k values
ingmat = matrix(c(logexp(oxy, params_bop, "cod", "maxin"),
                  logexp(oxy, params_bop_one, "cod", "maxin"),
                  logexp(oxy, params_bop_two, "cod", "maxin"),
                  logexp(oxy, params_bop_three, "cod", "maxin")),
                ncol = 4)
cod_scale$Ingestion = logexp(oxy, params_bop, "cod", "maxin")
cod_scale$Ingestion_min = apply(ingmat, 1, min)
cod_scale$Ingestion_max = apply(ingmat, 1, max)

# Maximum FCE
fce_bop = oxy_scenario(params_bop, t_max+1, 10, p = T)
fce_bop_one = oxy_scenario(params_bop_one, t_max+1, 10, p = T)
fce_bop_two = oxy_scenario(params_bop_two, t_max+1, 10, p = T)
fce_bop_three = oxy_scenario(params_bop_three, t_max+1, 10, p = T)
fcemax_cod = (sum(getHypoxiaDiet(fce_bop$sim, t_max+1, proportion = F)["cod",79,]) * (0.6*logexp(10, params_bop, "cod", "alpha")) - fce_bop$sim@params@metab["cod",79] * fce_bop$sim@params@other_params$metab_scale[t_max+1,"cod",79]) / sum(getHypoxiaDiet(fce_bop$sim, t_max+1, proportion = F)["cod",79,])
fcemax_cod_one = (sum(getHypoxiaDiet(fce_bop_one$sim, t_max+1, proportion = F)["cod",79,]) * (0.6*logexp(10, params_bop_one, "cod", "alpha")) - fce_bop_one$sim@params@metab["cod",79] * fce_bop_one$sim@params@other_params$metab_scale[t_max+1,"cod",79]) / sum(getHypoxiaDiet(fce_bop_one$sim, t_max+1, proportion = F)["cod",79,])
fcemax_cod_two = (sum(getHypoxiaDiet(fce_bop_two$sim, t_max+1, proportion = F)["cod",79,]) * (0.6*logexp(10, params_bop_two, "cod", "alpha")) - fce_bop_two$sim@params@metab["cod",79] * fce_bop_two$sim@params@other_params$metab_scale[t_max+1,"cod",79]) / sum(getHypoxiaDiet(fce_bop_two$sim, t_max+1, proportion = F)["cod",79,])
fcemax_cod_three = (sum(getHypoxiaDiet(fce_bop_three$sim, t_max+1, proportion = F)["cod",79,]) * (0.6*logexp(10, params_bop_three, "cod", "alpha")) - fce_bop_three$sim@params@metab["cod",79] * fce_bop_three$sim@params@other_params$metab_scale[t_max+1,"cod",79]) / sum(getHypoxiaDiet(fce_bop_three$sim, t_max+1, proportion = F)["cod",79,])
fcemax_flounder = (sum(getHypoxiaDiet(fce_bop$sim, t_max+1, proportion = F)["flounder",66,]) * (0.6*logexp(10, params_bop, "flounder", "alpha")) - fce_bop$sim@params@metab["flounder",66] * fce_bop$sim@params@other_params$metab_scale[t_max+1,"flounder",66]) / sum(getHypoxiaDiet(fce_bop$sim, t_max+1, proportion = F)["flounder",66,])
fcemax_flounder_one = (sum(getHypoxiaDiet(fce_bop_one$sim, t_max+1, proportion = F)["flounder",66,]) * (0.6*logexp(10, params_bop_one, "flounder", "alpha")) - fce_bop_one$sim@params@metab["flounder",66] * fce_bop_one$sim@params@other_params$metab_scale[t_max+1,"flounder",66]) / sum(getHypoxiaDiet(fce_bop_one$sim, t_max+1, proportion = F)["flounder",66,])
fcemax_flounder_two = (sum(getHypoxiaDiet(fce_bop_two$sim, t_max+1, proportion = F)["flounder",66,]) * (0.6*logexp(10, params_bop_two, "flounder", "alpha")) - fce_bop_two$sim@params@metab["flounder",66] * fce_bop_two$sim@params@other_params$metab_scale[t_max+1,"flounder",66]) / sum(getHypoxiaDiet(fce_bop_two$sim, t_max+1, proportion = F)["flounder",66,])
fcemax_flounder_three = (sum(getHypoxiaDiet(fce_bop_three$sim, t_max+1, proportion = F)["flounder",66,]) * (0.6*logexp(10, params_bop_three, "flounder", "alpha")) - fce_bop_three$sim@params@metab["flounder",66] * fce_bop_three$sim@params@other_params$metab_scale[t_max+1,"flounder",66]) / sum(getHypoxiaDiet(fce_bop_three$sim, t_max+1, proportion = F)["flounder",66,])

# Predict cod assimilation scaling for all BOP U and k values
foodmat_cod = matrix(nrow = length(oxy), ncol = 4)
metabmat_cod = matrix(nrow = length(oxy), ncol = 4)
foodmat_flounder = matrix(nrow = length(oxy), ncol = 4)
metabmat_flounder = matrix(nrow = length(oxy), ncol = 4)
for(i in 1:nrow(cod_scale)) {
  temp_bop = oxy_scenario(params_bop, t_max+1, cod_scale$oxy[i], p = T)
  temp_bop_one = oxy_scenario(params_bop_one, t_max+1, cod_scale$oxy[i], p = T)
  temp_bop_two = oxy_scenario(params_bop_two, t_max+1, cod_scale$oxy[i], p = T)
  temp_bop_three = oxy_scenario(params_bop_three, t_max+1, cod_scale$oxy[i], p = T)
  foodmat_cod[i,1] = sum(getHypoxiaDiet(temp_bop$sim, t_max+1, proportion = F)["cod",79,])
  foodmat_cod[i,2] = sum(getHypoxiaDiet(temp_bop_one$sim, t_max+1, proportion = F)["cod",79,])
  foodmat_cod[i,3] = sum(getHypoxiaDiet(temp_bop_two$sim, t_max+1, proportion = F)["cod",79,])
  foodmat_cod[i,4] = sum(getHypoxiaDiet(temp_bop_three$sim, t_max+1, proportion = F)["cod",79,])
  metabmat_cod[i,1] = temp_bop$sim@params@metab["cod",79] * temp_bop$sim@params@other_params$metab_scale[t_max+1,"cod",79]
  metabmat_cod[i,2] = temp_bop_one$sim@params@metab["cod",79] * temp_bop_one$sim@params@other_params$metab_scale[t_max+1,"cod",79]
  metabmat_cod[i,3] = temp_bop_two$sim@params@metab["cod",79] * temp_bop_two$sim@params@other_params$metab_scale[t_max+1,"cod",79]
  metabmat_cod[i,4] = temp_bop_three$sim@params@metab["cod",79] * temp_bop_three$sim@params@other_params$metab_scale[t_max+1,"cod",79]
  foodmat_flounder[i,1] = sum(getHypoxiaDiet(temp_bop$sim, t_max+1, proportion = F)["flounder",66,])
  foodmat_flounder[i,2] = sum(getHypoxiaDiet(temp_bop_one$sim, t_max+1, proportion = F)["flounder",66,])
  foodmat_flounder[i,3] = sum(getHypoxiaDiet(temp_bop_two$sim, t_max+1, proportion = F)["flounder",66,])
  foodmat_flounder[i,4] = sum(getHypoxiaDiet(temp_bop_three$sim, t_max+1, proportion = F)["flounder",66,])
  metabmat_flounder[i,1] = temp_bop$sim@params@metab["flounder",66] * temp_bop$sim@params@other_params$metab_scale[t_max+1,"flounder",66]
  metabmat_flounder[i,2] = temp_bop_one$sim@params@metab["flounder",66] * temp_bop_one$sim@params@other_params$metab_scale[t_max+1,"flounder",66]
  metabmat_flounder[i,3] = temp_bop_two$sim@params@metab["flounder",66] * temp_bop_two$sim@params@other_params$metab_scale[t_max+1,"flounder",66]
  metabmat_flounder[i,4] = temp_bop_three$sim@params@metab["flounder",66] * temp_bop_three$sim@params@other_params$metab_scale[t_max+1,"flounder",66]
}
alphamat = matrix(c(logexp(oxy, params_bop, "cod", "alpha"),
                    logexp(oxy, params_bop_one, "cod", "alpha"),
                    logexp(oxy, params_bop_two, "cod", "alpha"),
                    logexp(oxy, params_bop_three, "cod", "alpha")),
                  ncol = 4)
fcemat = (foodmat_cod * (alphamat*0.6) - metabmat_cod) / foodmat_cod
fcemat = sweep(fcemat, 2, c(fcemax_cod,fcemax_cod_one,fcemax_cod_two,fcemax_cod_three), "/")
cod_scale$FCE = fcemat[,1]
cod_scale$FCE_min = apply(fcemat, 1, min)
cod_scale$FCE_max = apply(fcemat, 1, max)

# Predict cod efficiency of reproduction scaling for all BOP U and k values
erpmat = matrix(c(logexp(oxy, params_bop, "cod", "erepro"),
                  logexp(oxy, params_bop_one, "cod", "erepro"),
                  logexp(oxy, params_bop_two, "cod", "erepro"),
                  logexp(oxy, params_bop_three, "cod", "erepro")),
                ncol = 4)
cod_scale$erp = logexp(oxy, params_bop, "cod", "erepro")
cod_scale$erp_min = apply(erpmat, 1, min)
cod_scale$erp_max = apply(erpmat, 1, max)

# Predict cod metabolic costs scaling for all BOP U and k values
metmat = matrix(c(logexp(oxy, params_bop, "cod", "met", met = T),
                  logexp(oxy, params_bop_one, "cod", "met", met = T),
                  logexp(oxy, params_bop_two, "cod", "met", met = T),
                  logexp(oxy, params_bop_three, "cod", "met", met = T)),
                ncol = 4)
cod_scale$met = logexp(oxy, params_bop, "cod", "met", met = T)
cod_scale$met_min = apply(metmat, 1, min)
cod_scale$met_max = apply(metmat, 1, max)

# Predict flounder activity scaling for all BOP U and k values
actmat = matrix(c(logexp(oxy, params_bop, "flounder", "search"),
                  logexp(oxy, params_bop_one, "flounder", "search"),
                  logexp(oxy, params_bop_two, "flounder", "search"),
                  logexp(oxy, params_bop_three, "flounder", "search")),
                ncol = 4)
flounder_scale$Activity = logexp(oxy, params_bop, "flounder", "search")
flounder_scale$Activity_min = apply(actmat, 1, min)
flounder_scale$Activity_max = apply(actmat, 1, max)

# Predict flounder consumption scaling for all BOP U and k values
ingmat = matrix(c(logexp(oxy, params_bop, "flounder", "maxin"),
                  logexp(oxy, params_bop_one, "flounder", "maxin"),
                  logexp(oxy, params_bop_two, "flounder", "maxin"),
                  logexp(oxy, params_bop_three, "flounder", "maxin")),
                ncol = 4)
flounder_scale$Ingestion = logexp(oxy, params_bop, "flounder", "maxin")
flounder_scale$Ingestion_min = apply(ingmat, 1, min)
flounder_scale$Ingestion_max = apply(ingmat, 1, max)

# Predict flounder FCE scaling for all BOP U and k values
alphamat = matrix(c(logexp(oxy, params_bop, "flounder", "alpha"),
                    logexp(oxy, params_bop_one, "flounder", "alpha"),
                    logexp(oxy, params_bop_two, "flounder", "alpha"),
                    logexp(oxy, params_bop_three, "flounder", "alpha")),
                  ncol = 4)
fcemat = (foodmat_flounder * alphamat - metabmat_flounder) / foodmat_flounder
fcemat = sweep(fcemat, 2, c(fcemax_flounder,fcemax_flounder_one,fcemax_flounder_two,fcemax_flounder_three), "/")
flounder_scale$FCE = fcemat[,1] / max(fcemat[,1])
flounder_scale$FCE_min = apply(fcemat, 1, min)
flounder_scale$FCE_max = apply(fcemat, 1, max)

# Predict flounder efficiency of reproduction scaling for all BOP U and k values
erpmat = matrix(c(logexp(oxy, params_bop, "flounder", "erepro"),
                  logexp(oxy, params_bop_one, "flounder", "erepro"),
                  logexp(oxy, params_bop_two, "flounder", "erepro"),
                  logexp(oxy, params_bop_three, "flounder", "erepro")),
                ncol = 4)
flounder_scale$erp = logexp(oxy, params_bop, "flounder", "erepro")
flounder_scale$erp_min = apply(erpmat, 1, min)
flounder_scale$erp_max = apply(erpmat, 1, max)

# Predict flounder metabolic costs scaling for all BOP U and k values
metmat = matrix(c(logexp(oxy, params_bop, "flounder", "met", met = T),
                  logexp(oxy, params_bop_one, "flounder", "met", met = T),
                  logexp(oxy, params_bop_two, "flounder", "met", met = T),
                  logexp(oxy, params_bop_three, "flounder", "met", met = T)),
                ncol = 4)
flounder_scale$met = logexp(oxy, params_bop, "flounder", "met", met = T)
flounder_scale$met_min = apply(metmat, 1, min)
flounder_scale$met_max = apply(metmat, 1, max)

# Rate colors
col.act = "#332288"
col.ing = "#117733"
col.fce = "#882255"
col.erp = "#DDCC77"
col.met = "#88CCEE"

# Initialize scaling plot
jpeg("Plots/all_scaling.jpg", width = 18, height = 15, units = 'cm', res = 600)
par(mfrow = c(2,2), mar = c(0.2,3,0.2,0), oma = c(5,2,3,10))

# Metabolic costs for cod
plot(met ~ oxy, data = cod_scale, axes = F,
     type = "l", lwd = 3, col = col.met,
     xlim = c(1,3), ylim = c(1,2), 
     xlab = "", ylab = "")
axis(2, at = seq(1,2,0.2), cex.axis = 1.2); box(); mtext("Cod", side = 3)
polygon(x = c(oxy,rev(oxy)),
        y = c(cod_scale$met_min,rev(cod_scale$met_max)),
        border = NA, col = adjustcolor(col.met, alpha.f = 0.5))

# Metabolic costs for flounder
plot(met ~ oxy, data = flounder_scale, axes = F,
     type = "l", lwd = 3, col = col.met,
     xlim = c(1,3), ylim = c(1,2), 
     xlab = "", ylab = "")
axis(2, at = seq(1,2,0.2), cex.axis = 1.2); box(); mtext("Flounder", side = 3)
polygon(x = c(oxy,rev(oxy)),
        y = c(flounder_scale$met_min,rev(flounder_scale$met_max)),
        border = NA, col = adjustcolor(col.met, alpha.f = 0.5))

# Down-scaling for cod
plot(erp ~ oxy, data = cod_scale, axes = F,
     type = "l", lwd = 3, col = col.erp,
     xlim = c(1,3), ylim = c(0,1), 
     xlab = "", ylab = "")
axis(1, at = seq(3), cex.axis = 1.2); axis(2, at = seq(0,1,0.2), labels = c(0,0.2,0.4,0.6,0.8,1), cex.axis = 1.2); box()
polygon(x = c(oxy,rev(oxy)),
        y = c(cod_scale$erp_min,rev(cod_scale$erp_max)),
        border = NA, col = adjustcolor(col.erp, alpha.f = 0.5))
lines(cod_scale$oxy, cod_scale$Activity,
      lwd = 3, col = col.act)
polygon(x = c(oxy,rev(oxy)),
        y = c(cod_scale$Activity_min,rev(cod_scale$Activity_max)),
        border = NA, col = adjustcolor(col.act, alpha.f = 0.5))
lines(cod_scale$oxy, cod_scale$Ingestion,
      lwd = 3, col = col.ing)
polygon(x = c(oxy,rev(oxy)),
        y = c(cod_scale$Ingestion_min,rev(cod_scale$Ingestion_max)),
        border = NA, col = adjustcolor(col.ing, alpha.f = 0.5))
lines(cod_scale$oxy, cod_scale$FCE,
      lwd = 3, col = col.fce)
polygon(x = c(oxy,rev(oxy)),
        y = c(cod_scale$FCE_min,rev(cod_scale$FCE_max)),
        border = NA, col = adjustcolor(col.fce, alpha.f = 0.5))

# Down-scaling for flounder
plot(erp ~ oxy, data = flounder_scale, axes = F,
     type = "l", lwd = 3, col = col.erp,
     xlim = c(1,3), ylim = c(0,1), 
     xlab = "", ylab = "")
axis(1, at = seq(3), cex.axis = 1.2); axis(2, at = seq(0,1,0.2), labels = c(0,0.2,0.4,0.6,0.8,1), cex.axis = 1.2); box()
polygon(x = c(oxy,rev(oxy)),
        y = c(flounder_scale$erp_min,rev(flounder_scale$erp_max)),
        border = NA, col = adjustcolor(col.erp, alpha.f = 0.5))
lines(flounder_scale$oxy, flounder_scale$Activity,
      lwd = 3, col = col.act)
polygon(x = c(oxy,rev(oxy)),
        y = c(flounder_scale$Activity_min,rev(flounder_scale$Activity_max)),
        border = NA, col = adjustcolor(col.act, alpha.f = 0.5))
lines(flounder_scale$oxy, flounder_scale$Ingestion,
      lwd = 3, col = col.ing)
polygon(x = c(oxy,rev(oxy)),
        y = c(flounder_scale$Ingestion_min,rev(flounder_scale$Ingestion_max)),
        border = NA, col = adjustcolor(col.ing, alpha.f = 0.5))
lines(flounder_scale$oxy, flounder_scale$FCE,
      lwd = 3, col = col.fce)
polygon(x = c(oxy,rev(oxy)),
        y = c(flounder_scale$FCE_min,rev(flounder_scale$FCE_max)),
        border = NA, col = adjustcolor(col.fce, alpha.f = 0.5))

# Add axis labels
mtext("Oxygen (mL/L)", side = 1, cex = 1.2, line = 3, outer = T)
mtext("Proportion of baseline", side = 2, cex = 1.2, outer = T)

# Add legend
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("right", bty = "n", c("Metabolic costs","Activity","Consumption","FCE","Egg survival"), cex = 1,
       lwd = 3, col = c(col.met,col.act,col.ing,col.fce,col.erp))
par(mfrow = c(1,1))

# Finish jpeg
dev.off()


########## X. Plot comparative yield ##########

# Analysis period
sim_ts = seq(1991,2019)

# Colors for plotting
color = c("dodgerblue","#E69F00","#009E73")

# Initialize plot
jpeg("Plots/yield_compare.jpg", width = 20, height = 18, units = 'cm', res = 600)
par(mfrow = c(2,2), oma = c(2,3,0,10), mar = c(4,4,3,1))

# Calibration and projection year indices
inds = (t_tune-length(sim_ts)+1):t_tune

# Observed catch for each species
yield_cod_obs = cod_catch[cod_catch$Year %in% sim_ts, ]$Catch/1000
yield_flounder_obs = flounder_catch[flounder_catch$Year %in% sim_ts, ]$Catch/1000
yield_sprat_obs = sprat_catch[sprat_catch$Year %in% sim_ts, ]$Catch/1000
yield_herring_obs = herring_catch[herring_catch$Year %in% sim_ts, ]$Catch/1000

# Calculate relative yield for metabolic prioritization models
yield_cod_bop = (getYield(sim_bop)[inds,"cod"]/1000/1000/1000) / yield_cod_obs
yield_flounder_bop = (getYield(sim_bop)[inds,"flounder"]/1000/1000/1000) / yield_flounder_obs
yield_sprat_bop = (getYield(sim_bop)[inds,"sprat"]/1000/1000/1000) / yield_sprat_obs
yield_herring_bop = (getYield(sim_bop)[inds,"herring"]/1000/1000/1000) / yield_herring_obs

# Calculate relative yield for models without metabolic prioritization
yield_cod_bom = (getYield(simcp_bom)[inds,"cod"]/1000/1000/1000) / yield_cod_obs
yield_flounder_bom = (getYield(simcp_bom)[inds,"flounder"]/1000/1000/1000) / yield_flounder_obs
yield_sprat_bom = (getYield(simcp_bom)[inds,"sprat"]/1000/1000/1000) / yield_sprat_obs
yield_herring_bom = (getYield(simcp_bom)[inds,"herring"]/1000/1000/1000) / yield_herring_obs

# Cod yield
plot(sim_ts, yield_cod_obs/yield_cod_obs, ylim = c(0,3), cex.main = 1.2, cex.axis = 1.2, col = color[1], type = 'l', lwd = 3, xlab = "", ylab = "", main = "Cod")
polygon(x = c(1991,2000,2000,1991), y = c(0,0,3,3), border = NA, col = rgb(0.5,0.5,0.5,0.2))
lines(sim_ts, yield_cod_bom, lwd = 3, col = color[2], lty = 5)
lines(sim_ts, yield_cod_bop, lwd = 3, col = color[3], lty = 3)

# Flounder yield
plot(sim_ts, yield_flounder_obs/yield_flounder_obs, ylim = c(0,3), cex.main = 1.2, cex.axis = 1.2, col = color[1], type = 'l', lwd = 3, xlab = "", ylab = "", main = "Flounder")
polygon(x = c(1991,2000,2000,1991), y = c(0,0,3,3), border = NA, col = rgb(0.5,0.5,0.5,0.2))
lines(sim_ts, yield_flounder_bom, lwd = 3, col = color[2], lty = 5)
lines(sim_ts, yield_flounder_bop, lwd = 3, col = color[3], lty = 3)

# Sprat yield
plot(sim_ts, yield_sprat_obs/yield_sprat_obs, ylim = c(0,3), cex.main = 1.2, cex.axis = 1.2, col = color[1], type = 'l', lwd = 3, xlab = "", ylab = "", main = "Sprat")
polygon(x = c(1991,2000,2000,1991), y = c(0,0,3,3), border = NA, col = rgb(0.5,0.5,0.5,0.2))
lines(sim_ts, yield_sprat_bom, lwd = 3, col = color[2], lty = 5)
lines(sim_ts, yield_sprat_bop, lwd = 3, col = color[3], lty = 3)

# Herring yield
plot(sim_ts, yield_herring_obs/yield_herring_obs, ylim = c(0,3), cex.main = 1.2, cex.axis = 1.2, col = color[1], type = 'l', lwd = 3, xlab = "", ylab = "", main = "Herring")
polygon(x = c(1991,2000,2000,1991), y = c(0,0,3,3), border = NA, col = rgb(0.5,0.5,0.5,0.2))
lines(sim_ts, yield_herring_bom, lwd = 3, col = color[2], lty = 5)
lines(sim_ts, yield_herring_bop, lwd = 3, col = color[3], lty = 3)

# Axis labels
mtext("Year", 1, outer = T, cex = 1.2, line = 0)
mtext("Yield at F (model/observed)", 2, outer = T, cex = 1.2, line = 1)

# Add legend
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("right", bty = 'n', c("Data","BOP w/ MP","BOM w/o MP"), col = color, lty = c(1,5,3), lwd = 3, cex = 1.2)
par(mfrow = c(1,1))

# Finish plot
dev.off()


########## XI. Plot comparative biomass ##########

# Analysis period
sim_ts = seq(1991,2019)

# Colors for plotting
color = c("dodgerblue","#E69F00","#009E73")

# Initialize plot
jpeg("Plots/ssb_compare.jpg", width = 20, height = 18, units = 'cm', res = 600)
par(mfrow = c(2,2), oma = c(2,3,0,10), mar = c(4,4,3,1))

# Calibration and projection year indicices
inds = (t_tune-length(sim_ts)+1):t_tune

# Observed SSB for each species
ssb_cod_obs = cod_ssb[cod_ssb$Year %in% sim_ts, ]$SSB/1000/1000/1000
ssb_flounder_obs = flounder_ssb[flounder_ssb$Year %in% sim_ts, ]$SSB/1000/1000/1000
ssb_sprat_obs = sprat_ssb[sprat_ssb$Year %in% sim_ts, ]$SSB/1000/1000/1000
ssb_herring_obs = herring_ssb[herring_ssb$Year %in% sim_ts, ]$SSB/1000/1000/1000

# Calculate relative ssb for metabolic prioritization models
ssb_cod_bop = rowSums((sweep(sim_bop@n[inds,1,], 2, params_bop@w*params_bop@dw, "*"))[,params_bop@w >= species_params(params_bop)$w_mat[1]])/1000/1000/1000 / ssb_cod_obs
ssb_flounder_bop = rowSums((sweep(sim_bop@n[inds,2,], 2, params_bop@w*params_bop@dw, "*"))[,params_bop@w >= species_params(params_bop)$w_mat[2]])/1000/1000/1000 / ssb_flounder_obs
ssb_sprat_bop = rowSums((sweep(sim_bop@n[inds,3,], 2, params_bop@w*params_bop@dw, "*"))[,params_bop@w >= species_params(params_bop)$w_mat[3]])/1000/1000/1000 / ssb_sprat_obs
ssb_herring_bop = rowSums((sweep(sim_bop@n[inds,4,], 2, params_bop@w*params_bop@dw, "*"))[,params_bop@w >= species_params(params_bop)$w_mat[4]])/1000/1000/1000 / ssb_herring_obs

# Calculate relative ssb for models without metabolic prioritization
ssb_cod_bom = rowSums((sweep(simcp_bom@n[inds,1,], 2, cphys_bom@w*cphys_bom@dw, "*"))[,cphys_bom@w >= species_params(cphys_bom)$w_mat[1]])/1000/1000/1000 / ssb_cod_obs
ssb_flounder_bom = rowSums((sweep(simcp_bom@n[inds,2,], 2, cphys_bom@w*cphys_bom@dw, "*"))[,cphys_bom@w >= species_params(cphys_bom)$w_mat[2]])/1000/1000/1000 / ssb_flounder_obs
ssb_sprat_bom = rowSums((sweep(simcp_bom@n[inds,3,], 2, cphys_bom@w*cphys_bom@dw, "*"))[,cphys_bom@w >= species_params(cphys_bom)$w_mat[3]])/1000/1000/1000 / ssb_sprat_obs
ssb_herring_bom = rowSums((sweep(simcp_bom@n[inds,4,], 2, cphys_bom@w*cphys_bom@dw, "*"))[,cphys_bom@w >= species_params(cphys_bom)$w_mat[4]])/1000/1000/1000 / ssb_herring_obs

# Cod SSB
plot(sim_ts, ssb_cod_obs/ssb_cod_obs, ylim = c(0,2), cex.main = 1.2, cex.axis = 1.2, col = color[1], type = 'l', lwd = 3, xlab = "", ylab = "", main = "Cod")
polygon(x = c(1991,2000,2000,1991), y = c(0,0,2,2), border = NA, col = rgb(0.5,0.5,0.5,0.2))
lines(sim_ts, ssb_cod_bom, lwd = 3, col = color[2], lty = 5)
lines(sim_ts, ssb_cod_bop, lwd = 3, col = color[3], lty = 3)

# Flounder SSB
plot(sim_ts, ssb_flounder_obs/ssb_flounder_obs, ylim = c(0,2), cex.main = 1.2, cex.axis = 1.2, col = color[1], type = 'l', lwd = 3, xlab = "", ylab = "", main = "Flounder")
polygon(x = c(1991,2000,2000,1991), y = c(0,0,2,2), border = NA, col = rgb(0.5,0.5,0.5,0.2))
lines(sim_ts, ssb_flounder_bom, lwd = 3, col = color[2], lty = 5)
lines(sim_ts, ssb_flounder_bop, lwd = 3, col = color[3], lty = 3)

# Sprat SSB
plot(sim_ts, ssb_sprat_obs/ssb_sprat_obs, ylim = c(0,2), cex.main = 1.2, cex.axis = 1.2, col = color[1], type = 'l', lwd = 3, xlab = "", ylab = "", main = "Sprat")
polygon(x = c(1991,2000,2000,1991), y = c(0,0,2,2), border = NA, col = rgb(0.5,0.5,0.5,0.2))
lines(sim_ts, ssb_sprat_bom, lwd = 3, col = color[2], lty = 5)
lines(sim_ts, ssb_sprat_bop, lwd = 3, col = color[3], lty = 3)

# Herring SSB
plot(sim_ts, ssb_herring_obs/ssb_herring_obs, ylim = c(0,2), cex.main = 1.2, cex.axis = 1.2, col = color[1], type = 'l', lwd = 3, xlab = "", ylab = "", main = "Herring")
polygon(x = c(1991,2000,2000,1991), y = c(0,0,2,2), border = NA, col = rgb(0.5,0.5,0.5,0.2))
lines(sim_ts, ssb_herring_bom, lwd = 3, col = color[2], lty = 5)
lines(sim_ts, ssb_herring_bop, lwd = 3, col = color[3], lty = 3)

# Axis labels
mtext("Year", 1, outer = T, cex = 1.2, line = 0)
mtext("SSB (model/observed)", 2, outer = T, cex = 1.2, line = 1)

# Add legend
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("right", bty = 'n', c("Data","BOP w/ MP","BOM w/o MP"), col = color, lty = c(1,5,3), lwd = 3, cex = 1.2)
par(mfrow = c(1,1))

# Finish plot
dev.off()