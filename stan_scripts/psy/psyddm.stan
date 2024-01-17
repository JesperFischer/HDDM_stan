functions{
  real normal_lub_rng(real mu, real sigma, real lb, real ub) {
    real p_lb = normal_cdf(lb, mu, sigma);
    real p_ub = normal_cdf(ub, mu, sigma);
    real u = uniform_rng(p_lb, p_ub);
    real y = mu + sigma * inv_Phi(u);
    return y;
  }
  
}

// based on codes/comments by Guido Biele, Joseph Burling, Andrew Ellis, and potentially others @ Stan mailing lists
data {
  
  int<lower=0> trials;
  real minRT;                        // minimum RT of the observed data
  real RT[trials];
  int<lower = 0, upper = 1> run_estimation; // a switch to evaluate the likelihood

  int<lower = 0, upper = 1> linear; // a switch to evaluate the likelihood
  array[trials] real x;
  array[trials] int resp;
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
  real<lower=0> alpha;  // boundary separation
  real<lower=0, upper=1> beta;   // initial bias
  real<lower=0> delta;  // drift rate
  real tau_raw;
  real psy_alpha;
  real<lower=0> psy_beta;


  
}

transformed parameters{
    
  real deltat[trials];
  real phi[trials];
  
  real tau = inv_logit(tau_raw) * minRT;
  

  if(!linear){
   for(i in 1:trials){
     phi[i] = 0.5+0.5*erf((x[i]-psy_alpha)/(psy_beta*sqrt(2)));
      
    if(phi[i] > 0.5){
      deltat[i] = (-(phi[i]*(1-phi[i]))*delta+0.25*delta);
    }else{
      deltat[i] = -(-(phi[i]*(1-phi[i]))*delta+0.25*delta);
    }
   }
  }else if (linear){
    for(i in 1:trials){
      phi[i] = 0.5+0.5*erf((x[i]-psy_alpha)/(psy_beta*sqrt(2)));
      
      deltat[i] = delta *  (phi[i] - (1 - phi[i]));
     
     }
  }
     
     
}


model {
  
  int c;

  
  psy_alpha ~ uniform(-40.5,40.5);
  
  target += normal_lpdf(psy_beta | 1, 5)-normal_lccdf(1 | 1, 5);;
  
  target += normal_lpdf(alpha | 1, 3)-normal_lccdf(1 | 1, 3);
  
  target += beta_proportion_lpdf(beta | 0.5, 10);
  
  if(linear){
   target += normal_lpdf(delta | 0, 5);
  }else if(!linear){
   target += normal_lpdf(delta | 0, 10);
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

  vector[trials] log_lik;
  int c;
  real prior_psy_alpha;
  real prior_psy_beta;
  real prior_alpha;
  real prior_beta;
  real prior_delta;
  real prior_tau_raw;
  
  prior_psy_alpha = uniform_rng(-40.5,40.5);
  prior_psy_beta = normal_lub_rng(1,5,0,100000);
  prior_alpha = normal_lub_rng(1,3,0,1000000);
  prior_beta = beta_proportion_rng(0.5,10);
  if(linear){
   prior_delta = normal_rng(0, 5);
  }else if(!linear){
   prior_delta = normal_rng(0, 10);
  }
  
  prior_tau_raw = normal_rng(0,1);
  
  
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