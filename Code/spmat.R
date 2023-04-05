# Create a matrix containing all species parameters

# Created: March 3, 2021
# Last modified: September 5, 2021 by EPD

# Contents (ctrl-f):
#	I. Load species values
#	II. Species and interaction matrix function


########## I. Load species values ##########

# Load life history values
cod_lh = read.table("Parameters/Cod/cod_lh.txt", header = T)
flounder_lh = read.table("Parameters/Flounder/flounder_lh.txt", header = T)
sprat_lh = read.table("Parameters/Sprat/sprat_lh.txt", header = T)
herring_lh = read.table("Parameters/Herring/herring_lh.txt", header = T)

# Load growth data
cod_ivb = read.table("Data/cod_ivb.txt", header = T)
flounder_ivb = read.table("Data/flounder_ivb.txt", header = T)
sprat_ivb = read.table("Data/sprat_ivb.txt", header = T)
herring_ivb = read.table("Data/herring_ivb.txt", header = T)

# Load growth data for projection years
cod_ivb_sim = read.table("Data/cod_ivb_sim.txt", header = T)
flounder_ivb_sim = read.table("Data/flounder_ivb_sim.txt", header = T)
sprat_ivb_sim = read.table("Data/sprat_ivb_sim.txt", header = T)
herring_ivb_sim = read.table("Data/herring_ivb_sim.txt", header = T)

# Load annual growth estimates
cod_ivb_cal = read.table("Parameters/Cod/cod_ivb_cal.txt", header = T)
sprat_ivb_cal = read.table("Parameters/Sprat/sprat_ivb_cal.txt", header = T)
herring_ivb_cal = read.table("Parameters/Herring/herring_ivb_cal.txt", header = T)

# Load annual growth estimates for projection years
cod_ivb_prj = read.table("Parameters/Cod/cod_ivb_prj.txt", header = T)
flounder_ivb_prj = read.table("Parameters/Flounder/flounder_ivb_prj.txt", header = T)
sprat_ivb_prj = read.table("Parameters/Sprat/sprat_ivb_prj.txt", header = T)
herring_ivb_prj = read.table("Parameters/Herring/herring_ivb_prj.txt", header = T)

# Create flounder annual data frame
flounder_ivb_cal = data.frame(Year = seq(1991,2000))
flounder_ivb_cal$k = rep(flounder_lh["k","stan"], nrow(flounder_ivb_cal))
flounder_ivb_cal$L_inf = rep(flounder_lh["L_inf","stan"], nrow(flounder_ivb_cal))
flounder_ivb_cal$W_inf = rep(flounder_lh["W_inf","stan"], nrow(flounder_ivb_cal))
flounder_ivb_cal$a = rep(flounder_lh["a","stan"], nrow(flounder_ivb_cal))
flounder_ivb_cal$b = rep(flounder_lh["b","stan"], nrow(flounder_ivb_cal))

# Load maturity values
cod_wmat = read.table("Parameters/Cod/cod_wmat.txt", header = T)
flounder_wmat = read.table("Parameters/Flounder/flounder_wmat.txt", header = T)
sprat_wmat = read.table("Parameters/Sprat/sprat_wmat.txt", header = T)
herring_wmat = read.table("Parameters/Herring/herring_wmat.txt", header = T)

# Load ppmr values
cod_diet = read.table("Parameters/Cod/cod_diet.txt", header = T)
flounder_diet = read.table("Parameters/Flounder/flounder_diet.txt", header = T)
sprat_diet = read.table("Parameters/Sprat/sprat_diet.txt", header = T)
herring_diet = read.table("Parameters/Herring/herring_diet.txt", header = T)

load("Parameters/Cod/cod_ppmr_lm.rda")
load("Parameters/Flounder/flounder_ppmr_lm.rda")
load("Parameters/Sprat/sprat_ppmr_lm.rda")
load("Parameters/Herring/herring_ppmr_lm.rda")

# Load food preference values
cod_preference = read.table("Parameters/Cod/cod_preference.txt", header = T)
flounder_preference = read.table("Parameters/Flounder/flounder_preference.txt", header = T)

# Load interaction matrix
interaction = as.matrix(read.table("Parameters/Interaction/interaction.txt", header = T))

# Load fishing mortality values
cod_selectivity = read.table("Parameters/Cod/cod_selectivity.txt", header = T)
flounder_selectivity = read.table("Parameters/Flounder/flounder_selectivity.txt", header = T)
sprat_selectivity = read.table("Parameters/Sprat/sprat_selectivity.txt", header = T)
herring_selectivity = read.table("Parameters/Herring/herring_selectivity.txt", header = T)

# Load calibration SSB
cod_ssb = read.table("Data/cod_ssb.txt", header = T)
flounder_ssb = read.table("Data/flounder_ssb.txt", header = T)
sprat_ssb = read.table("Data/sprat_ssb.txt", header = T)
herring_ssb = read.table("Data/herring_ssb.txt", header = T)
ssb_cal = read.table("Data/ssb_cal.txt", header = T)

# Load catch data
catch_cal = read.csv("Data/catch_cal.csv")


########## II. Species and interaction matrix function ##########

# Creates species_params data frame
#	rmax: maximum recruitment for cod, flounder, sprat, and herring, in that order
#	inter: the type of interaction desired; "preference" => based on preference data; "spatial" => based on overlap only
#	returns the species_params object for input to newMultispeciesParams
sp_mat = function(rmax) {
	
	# Enter parameters
	species_params = data.frame(species = c("cod", "flounder", "sprat", "herring"),
				w_inf = c(cod_lh["W_inf","stan"], flounder_lh["W_inf","stan"], sprat_lh["W_inf","stan"], herring_lh["W_inf","stan"]),
				w_mat = c(cod_wmat["L50",], flounder_wmat["L50",], sprat_wmat["L50",], herring_wmat["L50",]),
				k_vb = c(cod_lh["k","stan"], flounder_lh["k","stan"], sprat_lh["k","stan"], herring_lh["k","stan"]),
				w_min = rep(0.001, 4),
				beta = c(cod_diet$beta, flounder_diet$beta, sprat_diet$beta, herring_diet$beta),
				sigma = c(cod_diet$sigma, flounder_diet$sigma, sprat_diet$sigma, herring_diet$sigma),
				a = c(cod_lh["a","stan"], flounder_lh["a","stan"], sprat_lh["a","stan"], herring_lh["a","stan"]),
				b = c(cod_lh["b","stan"], flounder_lh["b","stan"], sprat_lh["b","stan"], herring_lh["b","stan"]),
				l25 = c(cod_selectivity$L25, NA, sprat_selectivity$L25, herring_selectivity$L25),
				l50 = c(cod_selectivity$L50, NA, sprat_selectivity$L50, herring_selectivity$L50),
				knife_edge_size = c(NA, flounder_selectivity$knife_edge_size, NA, NA),
				sel_func = c("sigmoid_length", "knife_edge", "sigmoid_length", "sigmoid_length"),
				catchability = c(cod_selectivity$Fmax, flounder_selectivity$Fbar, sprat_selectivity$Fmax, herring_selectivity$Fmax),
				R_max = rmax,
				biomass_observed = ssb_cal$ssb,
				biomass_cutoff = c(cod_wmat["L50",], flounder_wmat["L50",], sprat_wmat["L50",], herring_wmat["L50",]),
				yield_observed = catch_cal$catch,
				cA = c(1,1,1,1),
				cD = c(2,2,2,2),
				EA = c(1,1,1,1),
				ED = c(2,2,2,2),
				k = c(0.1,0.1,0.1,0.1),
				TD = c(5,5,5,5))
	
	# Available PP to species across both benthic and pelagic habitats
	species_params$interaction_resource = c(1, 1, 1, 1)

	return(list(species_params = species_params, interaction = interaction))
	
}