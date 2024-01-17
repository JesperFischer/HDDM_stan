functions{
  real normal_lub_rng(real mu, real sigma, real lb, real ub) {
    real p_lb = normal_cdf(lb| mu, sigma);
    real p_ub = normal_cdf(ub| mu, sigma);
    real u = uniform_rng(p_lb, p_ub);
    real y = mu + sigma * inv_Phi(u);
    return y;
  }
}

// based on codes/comments by Guido Biele, Joseph Burling, Andrew Ellis, and potentially others @ Stan mailing lists
data {
  
  int<lower=0> trials;

  array[trials] int resp;
  
  vector[trials] u;
  
}

parameters {

  real<lower=0> zeta;  // boundary separation
  real <lower =0, upper = 1> lr;


  
}

transformed parameters{
  vector[trials+1] expect;
  vector[trials] p_resp;
    

  expect[1] = 0.5;

    for(i in 1:trials){
      p_resp[i] = (expect[i]^zeta)/((expect[i]^zeta)+(1-expect[i])^zeta);
      expect[i+1] = expect[i]+lr*(u[i]-expect[i]);
  

  }
     
     
}


model {
  
  
  target += beta_lpdf(lr | 1,1);
  
  target += normal_lpdf(zeta | 1, 10)-normal_lccdf(0 | 1, 10);
  
    
  for(i in 1:trials){
    target += bernoulli_lpmf(resp[i] | p_resp[i]); 
  }
}


