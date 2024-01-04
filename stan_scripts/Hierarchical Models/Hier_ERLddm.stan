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
  
  
  array[trials] real stim;
  array[trials] int cue;
  
  array[trials] real percept;


}

transformed data{
  int idx = 1;
  int n_param = 9;
  
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
  array[trials] real<lower =0, upper =1> belief_to_cold;
  
  array[trials] real<lower =0, upper =1> mu_per;
    
    
  vector<lower=0>[S] alpha;
  vector<lower = 0, upper = 1> [S] beta;
  vector[S] delta;
  vector<lower=0>[S] tau;
  vector<lower=0, upper = 1> [S] lr;
  vector[S] nu;
  vector[S] sens_w;
  vector[S] sens_c;
  vector<lower=0>[S] prec_per;
  
  vector [S] inital_expect = rep_vector(0.5,S);
  
  matrix[S, n_param] indi = (diag_pre_multiply(tau_u, L_u) * z_expo)';
  
  delta = gm[1]+(indi[,1]);
  alpha = exp(gm[2]+(indi[,2]));
  beta = inv_logit(gm[3]+(indi[,3]));
  tau = inv_logit(gm[4]+(indi[,4])) .* minRT;
  lr = inv_logit(gm[5]+(indi[,5]));
  
  nu = gm[6]+(indi[,6]);
  sens_w = gm[7]+(indi[,7]);
  sens_c = gm[8]+(indi[,8]);
  prec_per = exp(gm[9]+(indi[,9]));

  
  expect[trial_per_par] = inital_expect;
  
  for (i in 1:trials){
    
    if(cue[i] == 1)
      belief_to_cold[i] = expect[i];
    else
      belief_to_cold[i] = 1-expect[i];

    if(stim[i] == 1){

      mu_per[i] = inv_logit(sens_c[S_id[i]]+((1-nu[S_id[i]])*stim[i]+nu[S_id[i]]*(belief_to_cold[i]-0.5)));
    
    }else{
      mu_per[i] = inv_logit(sens_w[S_id[i]]+(-((1-nu[S_id[i]])*(-stim[i])+nu[S_id[i]]*(belief_to_cold[i]-0.5))));
    
    }
    
    deltat[i] = delta[S_id[i]] *  (expect[i] - (1 - expect[i]));
    
    expect[i+1] = expect[i]+lr[S_id[i]]*(u[i]-expect[i]);
  }
  
}


model {
  
 int c;


  target += normal_lpdf(gm[1] | 0,4); //global mean of delta
  target += normal_lpdf(gm[2] | log(3),0.6); //global mean of alpha evaluated on the log scale as we exponentitate it.
  
  target += normal_lpdf(gm[3] | 0, 1); //global mean for beta where its inv_logit transformed i.e. 0 here is 0.5 bias
  target += normal_lpdf(gm[4] | 0, 2); //global mean for the tau parameter also inv_logit transformed
  target += normal_lpdf(gm[5] | -1, 1); //global mean for the tau parameter also inv_logit transformed
  
  target += normal_lpdf(gm[6] | 0, 2); //global mean for the tau parameter also inv_logit transformed
  target += normal_lpdf(gm[7] | 0, 2); //global mean for the tau parameter also inv_logit transformed
  target += normal_lpdf(gm[8] | 0, 2); //global mean for the tau parameter also inv_logit transformed
  target += normal_lpdf(gm[9] | 0, 3); //global mean for the tau parameter also inv_logit transformed
  
  
  target += std_normal_lpdf(to_vector(z_expo));
  target += normal_lpdf(tau_u | 0, 5)-normal_lccdf(0 | 0, 5);
  target += lkj_corr_cholesky_lpdf(L_u | 2);


    
    for(n in 1:trials){
      c = resp[n];
      target += beta_proportion_lpdf(percept[n] | mu_per[n],prec_per[S_id[n]]);

      if(c == 1){
        target += wiener_lpdf(RT[n] | alpha[S_id[n]], tau[S_id[n]], beta[S_id[n]], deltat[n]); 
      } else {
        target += wiener_lpdf(RT[n] | alpha[S_id[n]], tau[S_id[n]], 1-beta[S_id[n]], -deltat[n]);
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
  vector[S] prior_nu;
  vector[S] prior_sens_w;
  vector[S] prior_sens_c;
  vector[S] prior_prec_per;


  prior_gm[1] = normal_rng(0,4);
  prior_gm[2] = normal_rng(log(3),0.6);
  prior_gm[3] = normal_rng(0,1);
  prior_gm[4] = normal_rng(0,2);
  prior_gm[5] = normal_rng(-1,1);
  
  prior_gm[6] = normal_rng(0,2);
  prior_gm[7] = normal_rng(0,2);
  prior_gm[8] = normal_rng(0,2);
  prior_gm[9] = normal_rng(0,3);


  for(p in 1:n_param){
    prior_tau_u[p] = normal_lub_rng(0, 5, 0, 100000000);
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
    prior_nu[s] = (prior_gm[6]+((diag_pre_multiply(prior_tau_u, prior_L_u) * prior_z_expo)'[s,6]));
    prior_sens_w[s] = (prior_gm[7]+((diag_pre_multiply(prior_tau_u, prior_L_u) * prior_z_expo)'[s,7]));
    prior_sens_c[s] = (prior_gm[8]+((diag_pre_multiply(prior_tau_u, prior_L_u) * prior_z_expo)'[s,8]));
    prior_prec_per[s] = exp(prior_gm[9]+((diag_pre_multiply(prior_tau_u, prior_L_u) * prior_z_expo)'[s,9]));
    

    }
  }


  log_lik = rep_vector(0.0, trials);

  for(n in 1:trials){
    c = resp[n];
    
    log_lik[n] += beta_proportion_lpdf(percept[n] | mu_per[n],prec_per[S_id[n]]);

    if(c == 1){
      log_lik[n] += wiener_lpdf(RT[n] | alpha[S_id[n]], tau[S_id[n]], beta[S_id[n]], deltat[n]);
    } else {
      log_lik[n] += wiener_lpdf(RT[n] | alpha[S_id[n]], tau[S_id[n]], 1-beta[S_id[n]], -deltat[n]);
      }
  }

  
}
