// based on codes/comments by Guido Biele, Joseph Burling, Andrew Ellis, and potentially others @ Stan mailing lists
data {
  int<lower=0> Nu; // of upper boundary responses
  int<lower=0> Nl; // of lower boundary responses
  vector[Nu] RTu;    // upper boundary response times
  vector[Nl] RTl;    // lower boundary response times
  real minRT;      // minimum RT of the observed data
  real RTbound;    // lower bound or RT (e.g., 0.1 second)
  matrix[Nu,2] area_cold;
  matrix[Nl,2] area_warm;
}

parameters {
  vector<lower = 0>[2] b_cold_alpha;
  vector<lower = 0>[2] b_warm_alpha;
  
  vector<lower = 0>[2] b_cold_delta;
  vector<lower = 0>[2] b_warm_delta;
  
  vector<lower = 0>[2] b_cold_beta;
  vector<lower = 0>[2] b_warm_beta;
  
  
  real<lower=RTbound, upper=minRT> tau_cold;  // nondecision time
  real<lower=RTbound, upper=minRT> tau_warm;  // nondecision time


  real sigma_cold_alpha;
  real sigma_cold_delta;
  real sigma_cold_beta;
  
  real sigma_warm_alpha;
  real sigma_warm_delta;
  real sigma_warm_beta;
  
}



model {
  

  b_cold_delta ~ normal(0, 2);
  tau_cold ~ uniform(RTbound, minRT);
  b_cold_beta ~ normal(0,2);
  b_cold_alpha ~ normal(0,2);
  sigma_cold_alpha ~ normal(0,1);
  sigma_cold_delta ~ normal(0,1);
  sigma_cold_beta ~ normal(0,1);

  //alpha_warm ~ uniform(0, 5);
  
  
  
  b_warm_beta  ~ normal(0, 2);
  b_warm_delta ~ normal(0, 2);
  tau_warm ~ uniform(RTbound, minRT);
  b_warm_alpha ~ normal(0,2);
  sigma_warm_alpha ~ normal(0,1);
  sigma_warm_delta ~ normal(0,1);
  sigma_warm_beta ~ normal(0,1);


  
  RTu ~ wiener(exp(area_cold * b_cold_alpha+sigma_cold_alpha), tau_cold, inv_logit(area_cold * b_cold_beta+sigma_cold_beta), exp(area_cold * b_cold_delta+sigma_cold_delta));
  RTl ~ wiener(exp(area_warm * b_warm_alpha+sigma_warm_alpha), tau_warm, inv_logit(1- (area_warm * b_warm_beta+sigma_warm_beta)), -exp((area_warm * b_warm_delta+sigma_warm_delta)));
}

generated quantities {
  
  real prior_alpha; 
  real prior_beta;
  real prior_delta;
  real prior_tau;
  real log_lik;
  
  // For log likelihood calculation

  prior_alpha = uniform_rng(0, 5);
  prior_beta  = uniform_rng(0, 1);
  prior_delta = normal_rng(0, 2);
  prior_tau = uniform_rng(RTbound, minRT);

  // log_lik = wiener_lpdf(RTu | area_cold * b_cold, tau_cold, beta_cold, delta_cold);
  // log_lik += wiener_lpdf(RTl | area_warm * b_warm, tau_warm, 1-beta_warm, -delta_warm);

}

