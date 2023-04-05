# Initial values for Bayesian gadid meta-analysis models

# Initial values for the null hypothesis model - no response type effects
h0.inits = list(U = 0.6, 
                k = 1.5, 
                mu_sef = 1,
                sigma_sef = 0.1,
                sigma = 0.1)

# Initial values for alternative hypothesis model - response type included
h1.inits = list(U = c(0.6,0.8,0.3), 
                k = c(1.5,3.0,-5.0), 
                mu_sef = c(1,1,1),
                sigma_sef = c(0.1,0.1,0.1),
                sigma = rep(0.1,3))

# Initial values for alternative hypothesis model - response type with Michaelis-Menten type dynamics
h2.inits = list(k = c(1.8,0.3,2.2), 
                mu_sef = c(1,1,1),
                sigma_sef = c(0.1,0.1,0.1),
                sigma = rep(0.1,3))