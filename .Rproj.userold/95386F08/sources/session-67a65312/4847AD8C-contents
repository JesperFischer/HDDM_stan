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

  real minRT;                        // minimum RT of the observed data
  
  vector[trials] RT;
  
  array[trials] int resp;
  
  vector[trials] u;
  
  int<lower = 0, upper = 1> run_estimation;
  
  int<lower = 0, upper = 1> linear; // a switch to evaluate the likelihood
  
}

parameters {


  real<lower=0> alpha;  // boundary separation
  real<lower=0, upper=1> beta;   // initial bias
  real delta;  // drift rate
  real tau_raw;  // nondecision time
  real <lower =0, upper = 1> lr;


  
}

transformed parameters{
  vector[trials+1] expect;
  vector[trials] deltat;
    
    
      
  real tau = inv_logit(tau_raw) * minRT; // non-decision time at RT scale

  
  
  expect[1] = 0.5;
  if(!linear){
   for(i in 1:trials){
    if(expect[i] > 0.5){
      deltat[i] = (-(expect[i]*(1-expect[i]))*delta+0.25*delta);
    }else{
      deltat[i] = -(-(expect[i]*(1-expect[i]))*delta+0.25*delta);
    }
    expect[i+1] = expect[i]+lr*(u[i]-expect[i]);
   }
  }else if (linear){
    for(i in 1:trials){

      deltat[i] = delta *  (expect[i] - (1 - expect[i]));
      expect[i+1] = expect[i]+lr*(u[i]-expect[i]);

     }
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
  
  
  if(run_estimation==1){
    
    for(i in 1:trials){
      c = resp[i];
      

      if(c == 1){
        target += wiener_lpdf(RT[i] | alpha, tau, beta, deltat[i]); 
      } else {
        target += wiener_lpdf(RT[i] | alpha, tau, 1-beta, -deltat[i]);
        }
    }
  }
}


generated quantities{
  
  real prior_lr;
  real prior_alpha;
  real prior_beta;
  real prior_delta;
  real prior_tau;
  real prior_tau_raw;
  
  
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
  
  
  vector[trials] log_lik;
  int c;
  log_lik = rep_vector(0.0, trials);
  
  for(i in 1:trials){
    c = resp[i];
    if(c == 1){
      log_lik[i] += wiener_lpdf(RT[i] | alpha, tau, beta, deltat[i]); 
    } else {
      log_lik[i] += wiener_lpdf(RT[i] | alpha, tau, 1-beta, -deltat[i]);
      }
  }

  
}
