// functions{
//   real normal_lub_rng(real mu, real sigma, real lb, real ub) {
//     real p_lb = normal_cdf(lb| mu, sigma);
//     real p_ub = normal_cdf(ub| mu, sigma);
//     real u = uniform_rng(p_lb, p_ub);
//     real y = mu + sigma * inv_Phi(u);
//     return y;
//   }
// }

// based on codes/comments by Guido Biele, Joseph Burling, Andrew Ellis, and potentially others @ Stan mailing lists
data {
  
  int<lower=0> trials;
  real minRT;                        // minimum RT of the observed data
  
  vector[trials] RT;
  array[trials] int resp;
  
}

parameters {
  
  real raw_tau; // logit of non-decision time
  real raw_delta; // drift-rate
  real raw_alpha; // boundary separation
  real raw_beta; // starting point
}

transformed parameters{
  
  real tau = inv_logit(raw_tau)*minRT;
  real delta = raw_delta;
  real alpha = exp(raw_alpha);
  real beta = inv_logit(raw_beta);
  
  
}


model {
  
  int c;


  target += normal_lpdf(raw_tau | 0,1); //global mean for the tau parameter also inv_logit transformed
  
  target += normal_lpdf(raw_delta | 0,4); //global mean of delta evaluated on the real scale
  
  target += normal_lpdf(raw_alpha | log(3), 0.6); //global mean for beta where its on the log scale as we exponentitate it.

  target += normal_lpdf(raw_beta | 0, 1); //global mean for the beta parameter also inv_logit transformed
  

    for(n in 1:trials){
      c = resp[n];
      
      if(c == 1){
        target += wiener_lpdf(RT[n] | alpha, tau, beta, delta); 
      } else {
        target += wiener_lpdf(RT[n] | alpha, tau, 1-beta, -delta);
        }
    }
}


generated quantities{

  int c;
  real prior_alpha;
  real prior_delta;
  real prior_tau;
  real prior_beta;
  vector[trials] log_lik;

  prior_delta = normal_rng(0,4);
  prior_alpha = normal_rng(log(3),0.6);
  prior_tau = normal_rng(0,1);
  prior_beta = normal_rng(0,1);


  log_lik = rep_vector(0.0, trials);

  for(n in 1:trials){
    c = resp[n];
    if(c == 1){
        log_lik[n] += wiener_lpdf(RT[n] | alpha, tau, beta, delta); 
      } else {
        log_lik[n]  += wiener_lpdf(RT[n] | alpha, tau, 1-beta, -delta);
        }
      }
  }
