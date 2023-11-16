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
  
  real<lower=0>raw_tau; // logit of non-decision time
  real raw_delta; // drift-rate
  real<lower=0> raw_alpha; // boundary separation
  real<lower=0, upper=1> raw_beta; // starting point
}

transformed parameters{
  
  real tau = inv_logit(raw_tau)*minRT;
  real delta = raw_delta;
  real alpha = exp(raw_alpha);
  real beta = inv_logit(raw_beta);
  
  
}


model {
  
  int c;


  target += normal_lpdf(raw_tau | 0,1); //global mean of delta
  target += normal_lpdf(raw_delta | log(3),0.6); //global mean of delta evaluated on the log scale as we exponentitate it.
  
  target += normal_lpdf(raw_alpha | 0, 1); //global mean for beta where its inv_logit transformed i.e. 0 here is 0.5 bias
  target += normal_lpdf(raw_beta | 0, 2); //global mean for the tau parameter also inv_logit transformed
  

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

  real prior_alpha;
  real prior_delta;
  real prior_tau;
  real prior_beta;
  

  // hierarchical group level deviations
  vector<lower = 0>[4]  prior_tau_u;

  // Subject-level estimate matrix
  matrix[4, S] prior_z_expo;
  // for the cholesky decomposition
  cholesky_factor_corr[4] prior_L_u;

  vector[trials] log_lik;
  int c;


  vector<lower=0>[S] prior_alpha;
  vector<lower = 0, upper = 1> [S] prior_beta;
  vector[S] prior_delta;
  vector<lower=0>[S] prior_tau;




  prior_gm[1] = normal_rng(0,4);
  prior_gm[2] = normal_rng(log(3),0.6);
  prior_gm[3] = normal_rng(0,1);
  prior_gm[4] = normal_rng(0,2);



  for(p in 1:4){
    prior_tau_u[p] = normal_lub_rng(0, 3, 0, 100000000);
  }
  prior_L_u = lkj_corr_cholesky_rng(4,2);

  for(i in 1:4){
    for(s in 1:S){
    prior_z_expo[i,s] = std_normal_rng();
    prior_delta[s] = prior_gm[1]+((diag_pre_multiply(prior_tau_u, prior_L_u) * prior_z_expo)'[s,1]);
    prior_alpha[s] = exp(prior_gm[2]+((diag_pre_multiply(prior_tau_u, prior_L_u) * prior_z_expo)'[s,2]));
    prior_beta[s] = inv_logit(prior_gm[3]+((diag_pre_multiply(prior_tau_u, prior_L_u) * prior_z_expo)'[s,3]));
    prior_tau[s] = inv_logit(prior_gm[4]+((diag_pre_multiply(prior_tau_u, prior_L_u) * prior_z_expo)'[s,4]))* minRT[s];
    }
  }

  log_lik = rep_vector(0.0, trials);

  for(n in 1:trials){
    c = resp[n];
    if(c == 1){
      log_lik[n] += wiener_lpdf(RT[n] | alpha[S_id[n]], tau[S_id[n]], beta[S_id[n]], delta[S_id[n]]);
    } else {
      log_lik[n] += wiener_lpdf(RT[n] | alpha[S_id[n]], tau[S_id[n]], 1-beta[S_id[n]], -delta[S_id[n]]);
      }
  }


}
