
// based on codes/comments by Guido Biele, Joseph Burling, Andrew Ellis, and potentially others @ Stan mailing lists
data {
  
  int<lower=0> trials;
  real minRT;                        // minimum RT of the observed data
  vector[trials] RT;
  
  int<lower = 0, upper = 1> run_estimation; // a switch to evaluate the likelihood
  
  vector[trials] u;
  
  array[trials+1] int resp;
  
  int<lower = 0, upper = 1> linear; // a switch to evaluate the likelihood
  
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
  real <lower =0, upper = 1> lr;


  
}

transformed parameters{
  vector[trials+1] expect;
  vector[trials] deltat;
    
  
  real tau = inv_logit(tau_raw) * minRT;
  
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
  
  target += normal_lpdf(alpha | 1, 5)-normal_lccdf(1 | 1, 5);
  
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
