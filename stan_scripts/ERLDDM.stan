functions{
  real normal_lub_rng(real mu, real sigma, real lb, real ub) {
    real p_lb = normal_cdf(lb | mu, sigma);
    real p_ub = normal_cdf(ub | mu, sigma);
    real u = uniform_rng(p_lb, p_ub);
    real y = mu + sigma * inv_Phi(u);
    return y;
  }
}

// based on codes/comments by Guido Biele, Joseph Burling, Andrew Ellis, and potentially others @ Stan mailing lists
data {
  
  int<lower=0> trials;
  real minRT;                        // minimum RT of the observed data
  vector[trials] RT;
  
  vector[trials] u;
  
  array[trials+1] int resp;
  
  array[trials] real stim;
  array[trials] int cue;
  
  array[trials] real percept;


}

parameters {


  real<lower=0, upper=10> alpha;  // boundary separation
  real<lower=0, upper=1> beta;   // initial bias
  real delta;  // drift rate
  real tau_raw;  // nondecision time
  real <lower =0, upper = 1> lr;
  real  nu;
  real <lower =0> prec_per;
  real  sens_w;
  real  sens_c;
  
  
  
}

transformed parameters{
  array[trials+1] real<lower =0, upper =1> expect;
  array[trials] real uncert;
  array[trials] real<lower =0, upper =1> belief_to_cold;
  array[trials] real deltat;
  array[trials] real<lower =0, upper =1> mu_per;
  
  
  
  
  real tau = inv_logit(tau_raw) * minRT; // non-decision time at RT scale
  
  expect[1] = 0.5;
  for(i in 1:trials){

    if(cue[i] == 1)
      belief_to_cold[i] = expect[i];
    else
      belief_to_cold[i] = 1-expect[i];

    if(stim[i] == 1){

      mu_per[i] = inv_logit(sens_c+((1-nu)*stim[i]+nu*(belief_to_cold[i]-0.5)));
    
    }else{
      mu_per[i] = inv_logit(sens_w+(-((1-nu)*(-stim[i])+nu*(belief_to_cold[i]-0.5))));
    
    }
    
    
    uncert[i] = expect[i] - (1 - expect[i]);

    deltat[i] = delta * uncert[i];
    expect[i+1] = expect[i]+lr*(u[i]-expect[i]);

  }
  
  
}

model {
  int c;

  target += beta_proportion_lpdf(lr | 0.3,5);
  
  target += normal_lpdf(alpha | 1, 5)-normal_lccdf(0 | 1, 5);
  
  target += beta_proportion_lpdf(beta | 0.5, 10);
  
   target += normal_lpdf(delta | 0, 20);

  
  target += normal_lpdf(tau_raw | 0,1);
  
  target += lognormal_lpdf(prec_per | log(10),1);
  
  target += normal_lpdf(nu | 0,5);
  
  target += normal_lpdf(sens_c | 0,5);
  
  target += normal_lpdf(sens_w | 0,5);
  
  

    
    for(i in 1:trials){
      c = resp[i];
      target += beta_proportion_lpdf(percept[i] | mu_per[i],prec_per);
      if(c == 1){
        target += wiener_lpdf(RT[i] | alpha, tau, beta, deltat[i]); 
      } else {
        target += wiener_lpdf(RT[i] | alpha, tau, 1-beta, -deltat[i]);
        }
    }

}

generated quantities{
  
  vector[trials] log_lik;

  real prior_lr;
  real prior_alpha;
  real prior_beta;
  real prior_delta;
  real prior_tau;
  real prior_tau_raw;
  real prior_prec_per;
  real prior_nu;
  real prior_sens_w;
  real prior_sens_c;
  
  int c;
  
  vector[trials] pred_percept;

  
  
  prior_prec_per = lognormal_rng(log(10),1);
  prior_nu = normal_rng(0,5);
  prior_sens_w = normal_rng(0,5);
  prior_sens_c = normal_rng(0,5);
  
  prior_lr = beta_proportion_rng(0.3,5);
  
  prior_alpha = normal_lub_rng(1,5,0,10000);
  
  prior_beta = beta_proportion_rng(0.5,10);
  
  
  prior_delta = normal_rng(0, 20);

  prior_tau_raw = normal_rng(0,1);
  
  prior_tau = inv_logit(prior_tau_raw) * minRT;
  
  log_lik = rep_vector(0.0, trials);
  
  for(i in 1:trials){
      c = resp[i];
      pred_percept[i] = beta_proportion_rng(mu_per[i],prec_per);
      log_lik[i] += beta_proportion_lpdf(percept[i] | mu_per[i],prec_per);
    if(c == 1){
      log_lik[i] += wiener_lpdf(RT[i] | alpha, tau, beta, deltat[i]); 
    } else {
      log_lik[i] += wiener_lpdf(RT[i] | alpha, tau, 1-beta, -deltat[i]);
      }
  }
}
