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
  
  int<lower=0> S;
  
  array[trials] int S_id;
  
  vector[S] minRT;                        // minimum RT of the observed data
  
  vector[trials] RT;
  
  array[trials] int resp;
  
  array[S] int trial_per_par;
  
  vector[trials] u;
  
  int<lower = 0, upper = 1> run_estimation; // a switch to evaluate the likelihood
  
  int<lower = 0, upper = 1> linear; // a switch to evaluate the likelihood
  
}

transformed data{
  int idx = 1;
  int n_param = 5;
  
}

parameters {



  // hierarchical group level means 
  vector [n_param] gm;
  // hierarchical group level deviations
  vector<lower = 0>[n_param]  tau_u;
  // Subject-level estimate matrix 
  matrix[n_param, S] z_expo;
  // for the cholesky decomposition
  cholesky_factor_corr[n_param] L_u;
  


  
}

transformed parameters{
  vector[trials+1] expect;
  vector[trials] deltat;
    
    
  vector<lower=0>[S] alpha;
  vector<lower = 0, upper = 1> [S] beta;
  vector[S] delta;
  vector<lower=0>[S] tau;
  vector<lower=0, upper = 1> [S] lr;
  
  vector [S] inital_expect = rep_vector(0.5,S);
  
  
  
  delta = gm[1]+((diag_pre_multiply(tau_u, L_u) * z_expo)'[,1]);
  alpha = exp(gm[2]+((diag_pre_multiply(tau_u, L_u) * z_expo)'[,2]));
  beta = inv_logit(gm[3]+((diag_pre_multiply(tau_u, L_u) * z_expo)'[,3]));
  tau = inv_logit(gm[4]+((diag_pre_multiply(tau_u, L_u) * z_expo)'[,4])) .* minRT;
  lr = inv_logit(gm[5]+((diag_pre_multiply(tau_u, L_u) * z_expo)'[,5]));


  
  expect[trial_per_par] = inital_expect;
  
  for (n in 1:trials){
    if(!linear){
      if(expect[n] > 0.5){
      deltat[n] = (-(expect[n]*(1-expect[n]))*delta[S_id[n]]+0.25*delta[S_id[n]]);
    }else{
      deltat[n] = -(-(expect[n]*(1-expect[n]))*delta[S_id[n]]+0.25*delta[S_id[n]]);
      }
    }else if(linear){
      deltat[n] = delta[S_id[n]] *  (expect[n] - (1 - expect[n]));
    }
    expect[n+1] = expect[n]+lr[S_id[n]]*(u[n]-expect[n]);
  }
  
}


model {
  
 int c;


  target += normal_lpdf(gm[1] | 0,4); //global mean of delta
  target += normal_lpdf(gm[2] | log(3),0.6); //global mean of alpha evaluated on the log scale as we exponentitate it.
  
  target += normal_lpdf(gm[3] | 0, 1); //global mean for beta where its inv_logit transformed i.e. 0 here is 0.5 bias
  target += normal_lpdf(gm[4] | 0, 2); //global mean for the tau parameter also inv_logit transformed
  target += normal_lpdf(gm[5] | -1, 1); //global mean for the tau parameter also inv_logit transformed
  
  target += std_normal_lpdf(to_vector(z_expo));
  target += normal_lpdf(tau_u | 0, 2)-normal_lccdf(0 | 0, 2);
  target += lkj_corr_cholesky_lpdf(L_u | 2);



  if(run_estimation==1){
    
    for(n in 1:trials){
      c = resp[n];
      
      if(c == 1){
        target += wiener_lpdf(RT[n] | alpha[S_id[n]], tau[S_id[n]], beta[S_id[n]], deltat[n]); 
      } else {
        target += wiener_lpdf(RT[n] | alpha[S_id[n]], tau[S_id[n]], 1-beta[S_id[n]], -deltat[n]);
        }
    }
  }
}


generated quantities{
  
  vector[n_param] prior_gm;
  
  // hierarchical group level deviations
  vector<lower = 0>[n_param]  prior_tau_u;
  
  // Subject-level estimate matrix 
  matrix[n_param, S] prior_z_expo;
  // for the cholesky decomposition
  cholesky_factor_corr[n_param] prior_L_u;
  
  vector[trials] log_lik;
  
  int c;
  
  
  vector<lower=0>[S] prior_alpha;
  vector<lower = 0, upper = 1> [S] prior_beta;
  vector[S] prior_delta;
  vector<lower=0>[S] prior_tau;
  vector<lower=0>[S] prior_lr;
  
  
  prior_gm[1] = normal_rng(0,4);
  prior_gm[2] = normal_rng(log(3),0.6);
  prior_gm[3] = normal_rng(0,1);
  prior_gm[4] = normal_rng(0,2);
  prior_gm[5] = normal_rng(-1,1);
  
  
  for(p in 1:n_param){
    prior_tau_u[p] = normal_lub_rng(0, 2, 0, 100000000);
  }
  
  prior_L_u = lkj_corr_cholesky_rng(n_param,2);
  
  for(i in 1:n_param){
    for(s in 1:S){
    prior_z_expo[i,s] = std_normal_rng();
    prior_delta[s] = prior_gm[1]+((diag_pre_multiply(prior_tau_u, prior_L_u) * prior_z_expo)'[s,1]);
    prior_alpha[s] = exp(prior_gm[2]+((diag_pre_multiply(prior_tau_u, prior_L_u) * prior_z_expo)'[s,2]));
    prior_beta[s] = inv_logit(prior_gm[3]+((diag_pre_multiply(prior_tau_u, prior_L_u) * prior_z_expo)'[s,3]));
    prior_tau[s] = inv_logit(prior_gm[4]+((diag_pre_multiply(prior_tau_u, prior_L_u) * prior_z_expo)'[s,4]))* minRT[s];
    prior_lr[s] = inv_logit(prior_gm[5]+((diag_pre_multiply(prior_tau_u, prior_L_u) * prior_z_expo)'[s,5]));
    
    }
  }
  
  
  log_lik = rep_vector(0.0, trials);
  
  for(n in 1:trials){
    c = resp[n];
    if(c == 1){
      log_lik[n] += wiener_lpdf(RT[n] | alpha[S_id[n]], tau[S_id[n]], beta[S_id[n]], deltat[n]); 
    } else {
      log_lik[n] += wiener_lpdf(RT[n] | alpha[S_id[n]], tau[S_id[n]], 1-beta[S_id[n]], -deltat[n]);
      }
  }

  
}
