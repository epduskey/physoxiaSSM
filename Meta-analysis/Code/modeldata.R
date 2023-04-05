# Data for Bayesian gadid meta-analysis models

# Data for the null hypothesis model - no response type effects
dat.h0 = list(n_obs = nrow(gadid_ro),
              n_studs = length(unique(gadid_ro$study)),
              oxygen = gadid_ro$oxygen,
              response = gadid_ro$response,
              study = as.integer(gadid_ro$study),
              weights = 1/gadid_ro$SE_rel)

# Data for alternative hypothesis model - response type included
dat.h1 = list(n_obs = nrow(gadid_ro),
              n_studs = length(unique(gadid_ro$studytype)),
              n_types = length(unique(gadid_ro$type)),
              oxygen = gadid_ro$oxygen,
              response = gadid_ro$response,
              type = as.integer(gadid_ro$type),
              study = as.integer(gadid_ro$studytype),
              weights = 1/gadid_ro$SE_rel)

# Data for alternative hypothesis model - response type with Michaelis-Menten type dynamics
dat.h2 = list(n_obs = nrow(gadid_ro),
              n_studs = length(unique(gadid_ro$studytype)),
              n_types = length(unique(gadid_ro$type)),
              oxygen = gadid_ro$oxygen,
              response = gadid_ro$response,
              type = as.integer(gadid_ro$type),
              study = as.integer(gadid_ro$studytype),
              weights = 1/gadid_ro$SE_rel)