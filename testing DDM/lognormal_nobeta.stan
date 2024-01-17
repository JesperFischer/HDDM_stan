data {
  
  int<lower=0> trials;
  real minRT;                        // minimum RT of the observed data
  
  vector[trials] RT;
  array[trials] int resp;
  
}

parameters {
  real theta_raw;
  real intercept;
  real sigma_raw;
  real ndt_raw;
}



transformed parameters{

  real ndt = inv_logit(ndt_raw)*minRT;
  real sigma = exp(sigma_raw);
  real theta = inv_logit(theta_raw);
  

  
}
  
  

model{
  
  target += normal_lpdf(ndt_raw | 0,1); //global mean for the tau parameter also inv_logit transformed
  
  
  target += normal_lpdf(sigma_raw | log(3), 0.6); //global mean for beta where its on the log scale as we exponentitate it.

  target += normal_lpdf(theta_raw | 0, 1); //global mean for the beta parameter also inv_logit transformed
  
  target += normal_lpdf(intercept | 0,5); //global mean for the tau parameter also inv_logit transformed
  
  
  real mu_rt;
  mu_rt = intercept;

    
  for(n in 1:trials){
    target += bernoulli_lpmf(resp[n] | theta);
    target += lognormal_lpdf(RT[n] - ndt | mu_rt, sigma);
    
  }
  
}


generated quantities{

  //prior posterior updates
  // 

  vector[trials] log_lik;
  
  for (n in 1:trials){
    
    log_lik[n] = bernoulli_lpmf(resp[n] | theta) + 
                 lognormal_lpdf(RT[n] - ndt | intercept, sigma);
  
  }
  
}
