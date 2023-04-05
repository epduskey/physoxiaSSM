# Effects of hypoxia on gadiformes meta-analysis

# Created: 22 February 2023
# Last modified: 8 March 2023 by EPD

# Set working directory
setwd("/Users/epdus/OneDrive/Breathless/Code/Packets/physoxiaMizer/Meta-analysis")

# Load packages
library(rstan)
library(loo)
library(modelr)
library(HDInterval)

# Contents (ctrl-f):
# I. Load all data
# II. Organize data for models
# III. Run Bayesian models
# IV. Make predictions
# V. Plot model results


########## I. Load all data ##########

# Load in collated data
gadid = read.csv("Data/AllData.csv")
head(gadid)
str(gadid)

# Take a quick look at effects of oxygen on physiology and behavior
plot(response_activity ~ oxygen_mLL, data = gadid,
     xlab = "Oxygen (mL/L)", ylab = "Proportion max observed",
     xlim = c(0,10.5), ylim = c(-0.1,1.1),
     pch = 21, bg = "darkred", cex = 0.5/sqrt(gadid$SE_rel_activity))
points(gadid$oxygen_mLL, gadid$response_ingestion, pch = 22, bg = "darkblue", cex = 0.5/sqrt(gadid$SE_rel_ingestion))
points(gadid$oxygen_mLL, gadid$response_fce, pch = 23, bg = "darkgreen", cex = 0.5/sqrt(gadid$SE_rel_fce))
legend("bottomright", bty = "n", c("Activity","Ingestion","FCE"), pch = c(21,22,23), pt.bg = c("darkred","darkblue","darkgreen"), pt.cex = 1.2)


########## II. Organize data for models ##########

# Create new data frame with a single response variable
gadid_ro = data.frame(response = c(gadid$response_activity,gadid$response_ingestion,gadid$response_fce))

# Add response type as a factor
gadid_ro$type = factor(rep(c("activity","ingestion","fce"), each = nrow(gadid)))

# Add study as a factor
gadid_ro$study = factor(gadid$study)

# Add relative SE
gadid_ro$SE_rel = c(gadid$SE_rel_activity, gadid$SE_rel_ingestion, gadid$SE_rel_fce)

# Add oxygen levels
gadid_ro$oxygen = gadid$oxygen_mLL

# Remove missing values
gadid_ro = na.exclude(gadid_ro)

# Add study-response type as a factor
gadid_ro$studytype = factor(paste(gadid_ro$study,gadid_ro$type))

# Check out the results
str(gadid_ro)


########## III. Run Bayesian models ##########

# Source model scripts
source("Code/models.R")
source("Code/modeldata.R")
source("Code/modelinits.R")

# Compile models
gm.h0 = stan_model(model_code = stanH0, model_name = "stanH0")
gm.h1 = stan_model(model_code = stanH1, model_name = "stanH1")
gm.h2 = stan_model(model_code = stanH2, model_name = "stanH2")

# Run the null hypothesis model - no response type or age class effects
h0.stan = sampling(gm.h0,
  data = dat.h0,
  chains = 4,
  cores = 4,
  warmup = 5000,
  iter = 10000,
  refresh = 1000,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  init = function() {h0.inits})

# Run alternative hypothesis model - response type only
h1.stan = sampling(gm.h1,
                   data = dat.h1,
                   chains = 4,
                   cores = 4,
                   warmup = 5000,
                   iter  = 10000,
                   refresh = 1000,
                   control = list(adapt_delta = 0.99, max_treedepth = 15),
                   init = function() {h1.inits})

# Run alternative hypothesis model - response type with Michaelis-Menten type dynamics
h2.stan = sampling(gm.h2,
                   data = dat.h2,
                   chains = 4,
                   cores = 4,
                   warmup = 5000,
                   iter  = 10000,
                   refresh = 1000,
                   control = list(adapt_delta = 0.99, max_treedepth = 15),
                   init = function() {h2.inits})

# Check model diagnostics
check_hmc_diagnostics(h0.stan)
check_hmc_diagnostics(h1.stan)
check_hmc_diagnostics(h2.stan)

# Extract log likelihood for loo
h0_log_lik = extract_log_lik(h0.stan, merge_chains = F)
h1_log_lik = extract_log_lik(h1.stan, merge_chains = F)
h2_log_lik = extract_log_lik(h2.stan, merge_chains = F)

# LOO results for each model
( loo.h0 = loo(h0_log_lik, r_eff = relative_eff(exp(h0_log_lik)), pars = "log_lik") )
( loo.h1 = loo(h1_log_lik, r_eff = relative_eff(exp(h1_log_lik)), pars = "log_lik") )
( loo.h2 = loo(h2_log_lik, r_eff = relative_eff(exp(h2_log_lik)), pars = "log_lik") )

# Compare h0 and h1 using LOO
h0h1.compare = loo_compare(loo.h0, loo.h1)

# Compare h1 and h2 using LOO
h1h2.compare = loo_compare(loo.h1, loo.h2)

# Save the list of models and LOO results
results = list(h0.stan = h0.stan,
               h1.stan = h1.stan,
               h2.stan = h2.stan,
               loo.h0 = loo.h0,
               loo.h1 = loo.h1,
               loo.h2 = loo.h2,
               h0h1.compare = h0h1.compare,
               h1h2.compare = h1h2.compare)
save(object = results, file = "Output/meta_results.rda")

# # Load results
# load("Output/meta_results.rda")
# h0.stan = results$h0.stan
# h1.stan = results$h1.stan
# h2.stan = results$h2.stan
# loo.h0 = results$loo.h0
# loo.h1 = results$loo.h1
# loo.h2 = results$loo.h2
# h0h1.compare = results$h0h1.compare
# h1h2.compare = results$h1h2.compare


########## IV. Make predictions ##########

# Extract relevant parameters from all models
h0.pars = extract(h0.stan, pars = c("U","k","mu_sef"))
h1.pars = extract(h1.stan, pars = c("U","k","mu_sef"))
h2.pars = extract(h2.stan, pars = c("k","mu_sef"))

# Choose oxygen points for prediction
all.oxy = seq_range(gadid_ro$oxygen, n = 50)

# Make predictions with the null hypothesis model
h0.post = sapply(all.oxy, function(x) {(1/(1+exp(-h0.pars$U*(x-h0.pars$k))))*h0.pars$mu_sef})

# Make predictions with the alternative hypothesis model
h1.post.act = sapply(all.oxy, function(x) {(1/(1+exp(-h1.pars$U[,1]*(x-h1.pars$k[,1]))))*h1.pars$mu_sef[,1]})
h1.post.fce = sapply(all.oxy, function(x) {(1/(1+exp(-h1.pars$U[,2]*(x-h1.pars$k[,2]))))*h1.pars$mu_sef[,2]})
h1.post.ing = sapply(all.oxy, function(x) {(1/(1+exp(-h1.pars$U[,3]*(x-h1.pars$k[,3]))))*h1.pars$mu_sef[,3]})

# Make predictions with the Michaelis-Menten model
h2.post.act = sapply(all.oxy, function(x) {(x/(h2.pars$k[,1]+x))*h2.pars$mu_sef[,1]})
h2.post.fce = sapply(all.oxy, function(x) {(x/(h2.pars$k[,2]+x))*h2.pars$mu_sef[,2]})
h2.post.ing = sapply(all.oxy, function(x) {(x/(h2.pars$k[,3]+x))*h2.pars$mu_sef[,3]})

# Get mean for each rate
h0.mean = apply(h0.post, 2, mean)
h1.mean.act = apply(h1.post.act, 2, mean)
h1.mean.fce = apply(h1.post.fce, 2, mean)
h1.mean.ing = apply(h1.post.ing, 2, mean)
h2.mean.act = apply(h2.post.act, 2, mean)
h2.mean.fce = apply(h2.post.fce, 2, mean)
h2.mean.ing = apply(h2.post.ing, 2, mean)

# Get HDI for each rate
h0.hdi = apply(h0.post, 2, hdi)
h1.hdi.act = apply(h1.post.act, 2, hdi)
h1.hdi.fce = apply(h1.post.fce, 2, hdi)
h1.hdi.ing = apply(h1.post.ing, 2, hdi)
h2.hdi.act = apply(h2.post.act, 2, hdi)
h2.hdi.fce = apply(h2.post.fce, 2, hdi)
h2.hdi.ing = apply(h2.post.ing, 2, hdi)


########## V. Plot Bayesian model results ##########

# Rate colors
col.act = "#332288"
col.ing = "#117733"
col.fce = "#882255"
  
# Subset points for plotting among rate and age class
act.pts = subset(gadid_ro,type=="activity")
ing.pts = subset(gadid_ro,type=="ingestion")
fce.pts = subset(gadid_ro,type=="fce")

# Initialize jpeg
jpeg("meta_results.jpeg", width = 21, height = 9, units = "cm", res = 600)

# Initialize plot
par(mfrow = c(1,3), mar = c(5,3,3,0), oma = c(0,2,0,10))

# Initialize null hypothesis plot
plot(1, type = "n", 
     xlim = c(0.5,10.5), ylim = c(-0.15,1.15), 
     xlab = "", ylab = "", main = "H0",
     yaxt = "n", cex.axis = 1.3, cex.main = 1.3)
axis(side = 2, at = seq(0,1,0.2), cex.axis = 1.3)

# Add activity lines and points
lines(all.oxy, h0.mean, lwd = 3)
polygon(x = c(all.oxy,rev(all.oxy)),
        y = c(h0.hdi[1,],rev(h0.hdi[2,])),
        border = NA, col = adjustcolor("black", alpha.f = 0.5))
points(gadid_ro$oxygen, gadid_ro$response,
       col = "black", cex = 0.05/gadid_ro$SE_rel, lwd = 2)

# Initialize alternative hypothesis 1 plot
plot(1, type = "n", 
     xlim = c(0.5,10.5), ylim = c(-0.15,1.15), 
     xlab = "", ylab = "", main = "H1: Logistic",
     yaxt = "n", cex.axis = 1.3, cex.main = 1.3)
axis(side = 2, at = seq(0,1,0.2), cex.axis = 1.3)

# Add activity line and points
lines(all.oxy, h1.mean.act, lwd = 3, col = col.act)
polygon(x = c(all.oxy,rev(all.oxy)),
        y = c(h1.hdi.act[1,],rev(h1.hdi.act[2,])),
        border = NA, col = adjustcolor(col.act, alpha.f = 0.5))
points(act.pts$oxygen, act.pts$response,
       col = col.act, cex = 0.05/act.pts$SE_rel, lwd = 2)

# Add ingestion line and points
lines(all.oxy, h1.mean.ing, lwd = 3, col = col.ing)
polygon(x = c(all.oxy,rev(all.oxy)),
        y = c(h1.hdi.ing[1,],rev(h1.hdi.ing[2,])),
        border = NA, col = adjustcolor(col.ing, alpha.f = 0.5))
points(ing.pts$oxygen, ing.pts$response,
       col = col.ing, cex = 0.05/ing.pts$SE_rel, lwd = 2)

# Add FCE line and points
lines(all.oxy, h1.mean.fce, lwd = 3, col = col.fce)
polygon(x = c(all.oxy,rev(all.oxy)),
        y = c(h1.hdi.fce[1,],rev(h1.hdi.fce[2,])),
        border = NA, col = adjustcolor(col.fce, alpha.f = 0.5))
points(fce.pts$oxygen, fce.pts$response,
       col = col.fce, cex = 0.05/fce.pts$SE_rel, lwd = 2)

# Initialize alternative hypothesis 2 plot
plot(1, type = "n", 
     xlim = c(0.5,10.5), ylim = c(-0.15,1.15), 
     xlab = "", ylab = "", main = "H2: Michaelis-Menten",
     yaxt = "n", cex.axis = 1.3, cex.main = 1.3)
axis(side = 2, at = seq(0,1,0.2), cex.axis = 1.3)

# Add activity line and points
lines(all.oxy, h2.mean.act, lwd = 3, col = col.act)
polygon(x = c(all.oxy,rev(all.oxy)),
        y = c(h2.hdi.act[1,],rev(h2.hdi.act[2,])),
        border = NA, col = adjustcolor(col.act, alpha.f = 0.5))
points(act.pts$oxygen, act.pts$response,
       col = col.act, cex = 0.05/act.pts$SE_rel, lwd = 2)

# Add ingestion line and points
lines(all.oxy, h2.mean.ing, lwd = 3, col = col.ing)
polygon(x = c(all.oxy,rev(all.oxy)),
        y = c(h2.hdi.ing[1,],rev(h2.hdi.ing[2,])),
        border = NA, col = adjustcolor(col.ing, alpha.f = 0.5))
points(ing.pts$oxygen, ing.pts$response,
       col = col.ing, cex = 0.05/ing.pts$SE_rel, lwd = 2)

# Add FCE line and points
lines(all.oxy, h2.mean.fce, lwd = 3, col = col.fce)
polygon(x = c(all.oxy,rev(all.oxy)),
        y = c(h2.hdi.fce[1,],rev(h2.hdi.fce[2,])),
        border = NA, col = adjustcolor(col.fce, alpha.f = 0.5))
points(fce.pts$oxygen, fce.pts$response,
       col = col.fce, cex = 0.05/fce.pts$SE_rel, lwd = 2)

# Add axis labels
mtext("Oxygen (mL/L)", side = 1, cex = 1, line = -2, outer = T)
mtext("Proportion of baseline", side = 2, cex = 1, outer = T)

# Add legend
par(fig = c(0,1,0,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = T)
plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
legend("right", bty = "n", c("All rates","","Activity","Consumption","FCE"), cex = 1.2,
       lwd = 3, col = c("#000000","#FFFFFF",col.act,col.ing,col.fce))
par(mfrow = c(1,1))

# Finish jpeg
dev.off()