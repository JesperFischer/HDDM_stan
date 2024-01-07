data {
  int<lower=1> T;
  int<lower=1> S; //n sessions
  int<lower=1> P; //n participants
  vector[S*P] minRT;
   //Matrix (responses)
  matrix<lower=0>[S*P,T]  RT;
  array[S*P,T] int resp;
  
  vector[S*P] real_alpha;
  vector[S*P] real_beta;
  vector[S*P] real_delta;
  vector[S*P] real_tau;
}
transformed data{
  int N = 4;
  
}

parameters {

    // Group means (of parameters)
  vector[N] mu_g;
  // Between participant scales
  vector<lower = 0>[N]  tau_b;
  // Between participant cholesky decomposition
  cholesky_factor_corr[N] L_b;
  
  // Participant deviation (from the group level) 
  matrix[N, P] z_p;  
  
  // Within participant scales
  vector<lower = 0>[N]  tau_w;
  // Within participant cholesky decomposition
  cholesky_factor_corr[N] L_w;
  // Session deviation (from own group level)
  matrix[N, S*P] z_s; 
  
}

transformed parameters{
  
  
  matrix[N, P] delta_mu_p;  
  matrix[N, P] mu_p;  
  matrix[N, S*P] mu_p_rep;  
  matrix[N, S*P] fp_s;

  vector[S*P] tau;
  vector[S*P] delta;
  vector[S*P] alpha;
  vector[S*P] beta;  

  


  ///Recomposition (deviations from the group levels)
  
  delta_mu_p = diag_pre_multiply(tau_b, L_b) * z_p;

  
  for(idx in 1:P){
    mu_p[,idx] = mu_g + delta_mu_p[,idx];
  }
  
  
  // adding another as there are 2 sesssions
  
  mu_p_rep=append_col(append_col(append_col(mu_p,mu_p),mu_p),mu_p);
  
  // deviations from individual level


  fp_s = mu_p_rep + diag_pre_multiply(tau_w, L_w) * z_s;
  
  

  
  alpha = to_vector(exp(fp_s[1,]));
  tau = to_vector(inv_logit(fp_s[2,])) .* minRT;
  delta = to_vector(fp_s[3,]);
  beta = to_vector(inv_logit(fp_s[4,]));

  
  
  
}


model {
  
  int c;
  mu_g[1] ~ normal(0,1);
  mu_g[2] ~ normal(0,1);
  mu_g[3] ~ normal(0,1);
  mu_g[4] ~ normal(0,1);
  
  //Between participant scales
  tau_b[1] ~ normal(0,1);
  tau_b[2] ~ normal(0,1);
  tau_b[3] ~ normal(0,1);
  tau_b[4] ~ normal(0,1);

  //Between participant cholesky decomposition
  L_b ~ lkj_corr_cholesky(2);
  
  // Participant deviation 
  to_vector(z_p) ~ std_normal();  
  
  // Within participant scales
  tau_w[1] ~ normal(0,1);
  tau_w[2] ~ normal(0,1);
  tau_w[3] ~ normal(0,1);
  tau_w[4] ~ normal(0,1);

  // Within participant cholesky decomposition
  L_w ~ lkj_corr_cholesky(2);
  
  // Session deviation 
  to_vector(z_s) ~ std_normal();  
  
  for(p in 1:P){
    for(s in 1:S){
      for(n in 1:T){
          c = resp[s*p,n];
  
        if(c == 1){
          target += wiener_lpdf(RT[s*p,n] | alpha[s*p], tau[s*p], beta[s*p], delta[s*p]);
        } else {
          target += wiener_lpdf(RT[s*p,n] | alpha[s*p], tau[s*p], 1-beta[s*p], -delta[s*p]);
          }
        }
      }
    }
}

generated quantities{
  
  matrix[N,N] Corr_B;
  matrix[N,N] Corr_W;

  real resid_alpha_var = variance(real_alpha-alpha);
  real resid_tau_var = variance(real_tau-tau);
  real resid_delta_var = variance(real_delta-delta);
  real resid_beta_var = variance(real_beta-beta);

  vector[N] ICC = square(tau_b) ./ (square(tau_b) + square(tau_w));
  
  vector[N] ICC_alpha = square(tau_b) ./ (square(tau_b) + square(tau_w) + resid_alpha_var);
  
  vector[N] ICC_tau = square(tau_b) ./ (square(tau_b) + square(tau_w) + resid_tau_var);

  vector[N] ICC_delta = square(tau_b) ./ (square(tau_b) + square(tau_w) + resid_delta_var);
  
  vector[N] ICC_beta = square(tau_b) ./ (square(tau_b) + square(tau_w) + resid_beta_var);  

  Corr_B = L_b * L_b';
  Corr_W = L_w * L_w';

  
}
