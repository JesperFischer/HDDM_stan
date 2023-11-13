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
  
  int<lower = 0, upper = 1> run_estimation; // a switch to evaluate the likelihood
  
  vector[trials] u;
  
  array[trials+1] int resp;
  
  int<lower = 0, upper = 1> linear; // a switch to evaluate the likelihood
  
  
  array[trials] real stim;
  array[trials] int cue;
  
  array[trials] real percept;


}

parameters {
  // parameters of the DDM (parameter names in Ratcliffs DDM), from https://github.com/gbiele/stan_wiener_test/blob/master/stan_wiener_test.R
  // also see: https://groups.google.com/forum///!searchin/stan-users/wiener%7Csort:relevance/stan-users/-6wJfA-t2cQ/Q8HS-DXgBgAJ
  // alpha (a): Boundary separation or Speed-accuracy trade-off (high alpha means high accuracy). alpha > 0
  // beta (b): Initial bias Bias for either response (beta > 0.5 means bias towards "upper" response 'A'). 0 < beta < 1
  // delta (v): Drift rate Quality of the stimulus (delta close to 0 means ambiguous stimulus or weak ability). 0 < delta
  // tau (ter): Nondecision time + Motor response time + encoding time (high means slow encoding, execution). 0 < ter (in seconds)
  ///* upper boundary of tau must be smaller than minimum RT
  //to avoid zero likelihood for fast responses.
  //tau can for physiological reasone not be faster than 0.1 s.*/

  real<lower=0, upper=10> alpha;  // boundary separation
  real<lower=0, upper=1> beta;   // initial bias
  real delta;  // drift rate
  real tau_raw;  // nondecision time
  real <lower =0, upper = 1> lr;
  real <lower =0, upper = 1> nu;
  real <lower =0> prec_per;
  
  
}

transformed parameters{
  array[trials+1] real expect;
  array[trials] real uncert;
  array[trials] real belief_to_cold;
  array[trials] real deltat;
  array[trials] real mu_per;
  
  
  real tau = inv_logit(tau_raw) * minRT; // non-decision time at RT scale
  
  expect[1] = 0.5;
  for(i in 1:trials){

    if(cue[i] == 1)
      belief_to_cold[i] = expect[i];
    else
      belief_to_cold[i] = 1-expect[i];


    mu_per[i] = (1-nu)*stim[i]+nu*belief_to_cold[i];
    
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
  
  if(linear){
   target += normal_lpdf(delta | 0, 20);
  }else if(!linear){
   target += normal_lpdf(delta | 0, 20);
  }
  
  target += normal_lpdf(tau_raw | 0,1);
  
  target += lognormal_lpdf(prec_per | log(10),1);
  
  target += beta_proportion_lpdf(nu | 0.2,5);
  
  
  
  if(run_estimation==1){
    
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
  int c;
  
  vector[trials] pred_percept;

  
  
  prior_prec_per = lognormal_rng(log(10),1);
  prior_nu = beta_proportion_rng(0.2,5);
  
  prior_lr = beta_proportion_rng(0.3,5);
  
  prior_alpha = normal_lub_rng(1,5,0,10000);
  
  prior_beta = beta_proportion_rng(0.5,10);
  
  
  if(linear){
   prior_delta = normal_rng(0, 20);
  }else if(!linear){
   prior_delta = normal_rng(0, 20);
  }
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