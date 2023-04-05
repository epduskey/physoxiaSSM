# All Bayesian models for gadid meta-analysis

# Write null hypothesis Stan model - no response type effects
stanH0 = "
data
{
    // Number of observations
    int n_obs;
    
    // Number of studies
    int n_studs;
    
    // Oxygen
    vector[n_obs] oxygen;
    
    // Response
    vector[n_obs] response;
    
    // Mapping variable for study
    array[n_obs] int study;
    
    // Regression weights
    vector[n_obs] weights;
}
parameters
{
    // U and k values
    real<lower=0> U;
    real k;
    
    // Study effect
    vector[n_studs] sef;
    
    // Study effect mean and standard deviation
    real<lower=0> mu_sef;
    real<lower=0> sigma_sef;
    
    // Population standard deviation
    real<lower=0> sigma;
}
transformed parameters
{
    // Mean vector
    vector[n_obs] mu;
    
    // Calculate mean
    for(i in 1:n_obs) {
      mu[i] = (1 / (1 + exp(-U*(oxygen[i] - k)))) * sef[study[i]];
    }
}
model
{
    // Priors
    U ~ normal(0,1000) T[0,];
    k ~ normal(0,1000);
    mu_sef ~ normal(0,1000);
    sigma_sef ~ gamma(0.001,0.001);
    sef ~ normal(mu_sef,sigma_sef);
    sigma ~ gamma(0.001,0.001);
    
    // Likelihood
    target += normal_lpdf(response | mu, sigma) * weights;
}
generated quantities 
{
    // Calculate log likelihood
    vector[n_obs] log_lik;
    for (i in 1:n_obs) {
      log_lik[i] = normal_lpdf(response[i] | mu[i], sigma);
    }
}
"

# Write alternative hypothesis Stan model - response type included
stanH1 = "
data
{
    // Number of observations
    int n_obs;
    
    // Number of studies for each rate combined
    int n_studs;
    
    // Number of response types
    int n_types;
    
    // Oxygen
    vector[n_obs] oxygen;
    
    // Response
    vector[n_obs] response;
    
    // Mapping variable for response type
    array[n_obs] int type;
    
    // Mapping variable for study
    array[n_obs] int study;
    
    // Regression weights
    vector[n_obs] weights;
}
parameters
{
    // U and k values
    vector<lower=0>[n_types] U;
    vector[n_types] k;
    
    // Study effect
    vector[n_studs] sef;
    
    // Study effect mean and standard deviation
    vector<lower=0>[n_types] mu_sef;
    vector<lower=0>[n_types] sigma_sef;
    
    // Population standard deviation by rate
    vector<lower=0>[n_types] sigma;
}
transformed parameters
{
    // Mean vector
    vector[n_obs] mu;
    
    // Calculate mean
    for(i in 1:n_obs) {
      mu[i] = (1 / (1 + exp(-U[type[i]]*(oxygen[i] - k[type[i]])))) * sef[study[i]];
    }
}
model
{
    // Priors
    for(i in 1:n_types) U[i] ~ normal(0,1000) T[0,];
    k ~ normal(0,1000);
    mu_sef ~ normal(0,1000);
    sigma_sef ~ gamma(0.001,0.001);
    sef[study] ~ normal(mu_sef[type],sigma_sef[type]);
    sigma ~ gamma(0.001,0.001);
    
    // Likelihood
    target += normal_lpdf(response | mu, sigma[type]) * weights;
}
generated quantities 
{
    // Calculate log likelihood
    vector[n_obs] log_lik;
    for (i in 1:n_obs) {
      log_lik[i] = normal_lpdf(response[i] | mu[i], sigma[type[i]]);
    }
}
"

# Write alternative hypothesis Stan model - response type with Michaelis-Menten type dynamics
stanH2 = "
data
{
    // Number of observations
    int n_obs;
    
    // Number of studies for each rate combined
    int n_studs;
    
    // Number of response types
    int n_types;
    
    // Oxygen
    vector[n_obs] oxygen;
    
    // Response
    vector[n_obs] response;
    
    // Mapping variable for response type
    array[n_obs] int type;
    
    // Mapping variable for study
    array[n_obs] int study;
    
    // Regression weights
    vector[n_obs] weights;
}
parameters
{
    // k values
    vector<lower=0>[n_types] k;
    
    // Study effect
    vector[n_studs] sef;
    
    // Study effect mean and standard deviation
    vector<lower=0>[n_types] mu_sef;
    vector<lower=0>[n_types] sigma_sef;
    
    // Population standard deviation by rate
    vector<lower=0>[n_types] sigma;
}
transformed parameters
{
    // Mean vector
    vector[n_obs] mu;
    
    // Calculate mean
    for(i in 1:n_obs) {
      mu[i] = (oxygen[i] / (k[type[i]] + oxygen[i])) * sef[study[i]];
    }
}
model
{
    // Priors
    for(i in 1:n_types) k[i] ~ normal(0,1000) T[0,];
    mu_sef ~ normal(0,1000);
    sigma_sef ~ gamma(0.001,0.001);
    sef[study] ~ normal(mu_sef[type],sigma_sef[type]);
    sigma ~ gamma(0.001,0.001);
    
    // Likelihood
    target += normal_lpdf(response | mu, sigma[type]) * weights;
}
generated quantities 
{
    // Calculate log likelihood
    vector[n_obs] log_lik;
    for (i in 1:n_obs) {
      log_lik[i] = normal_lpdf(response[i] | mu[i], sigma[type[i]]);
    }
}
"