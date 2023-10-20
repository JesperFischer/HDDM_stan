//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> N;
  array[N] int resp;
  vector[N] u;
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  real<lower=0, upper = 1> alpha;
  real<lower=0> zeta;
}

transformed parameters{

  vector[N] p_resp;
  vector[N+1] m;

  m[1] = 0.5;
  
  for(i in 1:N){
    
    m[i+1] = m[i]+alpha*(u[i]-m[i]);
    
    p_resp[i] = (inv_logit(m[i])^zeta)/((inv_logit(m[i])^zeta)+(1-inv_logit(m[i]))^zeta);
  }
  
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  
  alpha ~ beta_proportion(0.5,5);
  zeta ~ lognormal(log(10),0.5);

  
  for(i in 1:N){
    resp[i] ~ bernoulli(p_resp[i]);
  }
}


generated quantities{
  vector[N] log_lik;
  
  for(i in 1:N){
    log_lik[i] = bernoulli_lpmf(resp[i] | p_resp[i]);
  }
  
}
