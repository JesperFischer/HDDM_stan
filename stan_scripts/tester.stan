// Stan code for RL-DDM model
// written by Jihyun Hur
// model = hba_two_v (model 2)

data {
  int<lower=1> N;                    // subject
  int<lower=1> T;                    // trial
  int<lower=1> Steps;                // total number of steps
  int<lower=0, upper=1> C[T, N];     // choice
  int<lower=0, upper=1> E[T, N];     // escape
  int<lower=1, upper=4> cond[T, N];  // condition
  int<lower=0, upper=1> fd[T, N];    // feedback
  int<lower=1, upper=2> grp[N];      // group
  real rt[T, N];                     // response time
  real minRT[N];                     // minimum RT for each subject of the observed data
  real RTbound;                      // minimum threshold for response time of a go action
  real RTmax;                        // max RT
}

transformed data {
  real<lower=0> n_sui;
  real<lower=0> n_nosui;
  vector[N] is_sui;
  vector[N] is_nosui;
  
  n_sui = 0;
  n_nosui = 0;
  
  is_sui = rep_vector(0, N);
  is_nosui = rep_vector(0, N);
  
  // group size
  for (i in 1:N) {
    if (grp[i] == 1) {
      n_sui += 1;
      is_sui[i] = 1;
    } else {
      n_nosui += 1;
      is_nosui[i] = 1;
    }
  }
}

parameters {
  // hyperparameters
  vector[7] mu_pr;
  vector<lower=0>[7] sigma;
  
  // subject-level raw parameters (for Matt trick)
  vector[N] tau_pr;                   
  vector[N] b1_pr;                    
  vector[N] b2_pr;
  vector[N] b3_pr;
  vector[N] w_pr;                    
  vector[N] alpha_pr;                 
  vector[N] omega_pr;               
}

transformed parameters {
  // transform subject-level raw parameters
  vector<lower=RTbound, upper=max(minRT)>[N] tau;             // non-decision time
  vector<lower=-20, upper=20>[N] b1;                          // constant go bias in escape cond
  vector<lower=-20, upper=20>[N] b2;                          // constant go bias in avoid cond
  vector<lower=-20, upper=20>[N] b3;                          // drift rate differential action value weight
  vector<lower=0.001, upper=0.999>[N] w;                      // starting point for both conditions
  vector<lower=0, upper=1>[N] alpha;                          // learning rate for state-action values
  vector<lower=0>[N] omega;                                   // decision threshold
  
  for (i in 1:N) {
    tau[i] = Phi_approx(mu_pr[1] + sigma[1] * tau_pr[i]) * (minRT[i] - RTbound) + RTbound; 
    w[i] = Phi_approx(mu_pr[5] + sigma[5] * w_pr[i]);
    alpha[i] = Phi_approx(mu_pr[6] + sigma[6] * alpha_pr[i]);
  }
  b1 = mu_pr[2] + sigma[2] * b1_pr;
  b2 = mu_pr[3] + sigma[3] * b2_pr;
  b3 = mu_pr[4] + sigma[4] * b3_pr;
  omega = exp(mu_pr[7] + sigma[7] * omega_pr);
}

model {
  real step_size;
  real temp_rt;
  real prob_density;
  real nogo_lik;
  
  // hyperparameters
  mu_pr ~ normal(0, 1);
  sigma ~ normal(0, 0.5);                           
  
  // individual parameters for non-centered parameterization
  tau_pr ~ normal(0, 1);                   
  b1_pr ~ normal(0, 1);                    
  b2_pr ~ normal(0, 1);                     
  b3_pr ~ normal(0, 1);                    
  w_pr ~ normal(0, 1);                    
  alpha_pr ~ normal(0, 1);                 
  omega_pr ~ normal(0, 1);
  
  // subject loop
  for (i in 1:N) {
    // assign individual variables
    int c;
    int s;
    
    // initial values for state-action values
    matrix[4, 2] Q;          
    Q = rep_matrix(0.0, 4, 2);                  
    
    for (t in 1:T) {
      real v;
      c = C[t, i];                               
      s = cond[t, i];                            
      
      // different intercepts for each condition
      v = (b1[i] * E[t, i]) + (b2[i] * (1 - E[t, i])) + (b3[i] * (Q[s, 2] - Q[s, 1]));
      
      // update posterior
      if (c == 1) {
        rt[t, i] ~ wiener(omega[i], tau[i], w[i], v); 
      } else {
        nogo_lik = 0.0;
        
        for (n in 1:Steps) {
          step_size = (n * 1.0) / Steps;
          temp_rt = step_size * (RTmax - tau[i]) + tau[i];
          prob_density = exp(wiener_lpdf(temp_rt | omega[i], tau[i], 1-w[i], -v));
          nogo_lik += ((RTmax - tau[i]) / Steps) * prob_density;
        }
        
      if (is_nan(nogo_lik) == 1) {
       nogo_lik = 2.220446e-16;
      }
      
      target += log(nogo_lik);
      }
      
      // update Q-value
      Q[s, c+1] += alpha[i]*(fd[t, i] - Q[s, c+1]);
    }
  }
}

generated quantities {
  // for group level parameters
  real<lower=RTbound, upper=max(minRT)> mu_tau;         // nondecision time
  real<lower=-20, upper=20>             mu_b1;          // constant go bias in escape cond
  real<lower=-20, upper=20>             mu_b2;          // constant go bias in avoid cond
  real<lower=-20, upper=20>             mu_b3;          // drift rate differential action value weight
  real<lower=0, upper=1>                mu_w;           // starting point for both conditions
  real<lower=0, upper=1>                mu_alpha;       // learning rate for state-action values
  real<lower=0>                         mu_omega;       // decision threshold
  /*
  real<lower=RTbound, upper=max(minRT)> mu_tau_sui;           
  real<lower=RTbound, upper=max(minRT)> mu_tau_nosui;
  */
  real<lower=-20, upper=20>             mu_b1_sui;            
  real<lower=-20, upper=20>             mu_b1_nosui;          
  real<lower=-20, upper=20>             mu_b2_sui;            
  real<lower=-20, upper=20>             mu_b2_nosui;
  real<lower=-20, upper=20>             mu_b3_sui;            
  real<lower=-20, upper=20>             mu_b3_nosui;
  real<lower=0, upper=1>                mu_w_sui;            
  real<lower=0, upper=1>                mu_w_nosui;
  /*
  real<lower=0, upper=1>                mu_alpha_sui;         
  real<lower=0, upper=1>                mu_alpha_nosui;       
  real<lower=0>                         mu_omega_sui;         
  real<lower=0>                         mu_omega_nosui;         
  */
  real log_lik[N];
  real nogo_lik;
  real step_size;
  real temp_rt;
  real prob_density;
  
  // calculate group means of each parameter
  /*
  mu_tau_sui = sum(is_sui .* tau) / n_sui;
  mu_tau_nosui = sum(is_nosui .* tau) / n_nosui;
  */
  mu_b1_sui = sum(is_sui .* b1) / n_sui;
  mu_b1_nosui = sum(is_nosui .* b1) / n_nosui;
  mu_b2_sui = sum(is_sui .* b2) / n_sui;
  mu_b2_nosui = sum(is_nosui .* b2) / n_nosui;
  mu_b3_sui = sum(is_sui .* b3) / n_sui;
  mu_b3_nosui = sum(is_nosui .* b3) / n_nosui;
  mu_w_sui = sum(is_sui .* w) / n_sui;
  mu_w_nosui = sum(is_nosui .* w) / n_nosui;
  /*
  mu_alpha_sui = sum(is_sui .* alpha) / n_sui;
  mu_alpha_nosui = sum(is_nosui .* alpha) / n_nosui;
  mu_omega_sui = sum(is_sui .* omega) / n_sui;
  mu_omega_nosui = sum(is_nosui .* omega) / n_nosui;
  */
  // assign group level parameter values
  mu_tau = Phi_approx(mu_pr[1]) * (mean(minRT) - RTbound) + RTbound; 
  mu_b1 = mu_pr[2];
  mu_b2 = mu_pr[3];
  mu_b3 = mu_pr[4];
  mu_w = Phi_approx(mu_pr[5]);
  mu_alpha = Phi_approx(mu_pr[6]);
  mu_omega = exp(mu_pr[7]);

 {// for log likelihood calculation (local)
  for (i in 1:N) {
    // assign individual variables
    int c;
    int s;
    matrix[4, 2] Q;
    
    // initial values for state-action values
    Q = rep_matrix(0.0, 4, 2);      
    log_lik[i] = 0.0;
    
    for (t in 1:T) {
      real v;
      c = C[t, i];                                                                     
      s = cond[t, i];                                                                  
      
      v = (b1[i] * E[t, i]) + (b2[i] * (1 - E[t, i])) + (b3[i] * (Q[s, 2] - Q[s, 1]));
      
      // calculate likelihood
      if (c == 1) {
        log_lik[i] += wiener_lpdf(rt[t, i] | omega[i], tau[i], w[i], v);
        } else {
          nogo_lik = 0.0;
          
          for (n in 1:Steps) {
            step_size = (n * 1.0) / Steps;
            temp_rt = step_size * (RTmax - tau[i]) + tau[i];
            prob_density = exp(wiener_lpdf(temp_rt | omega[i], tau[i], 1-w[i], -v));
            nogo_lik += ((RTmax - tau[i]) / Steps) * prob_density;
            }
            
            if (is_nan(nogo_lik) == 1) {
              nogo_lik = 2.220446e-16;
            }
              
          log_lik[i] += log(nogo_lik);
        }
      
      // update Q-value
      Q[s, c+1] += alpha[i]*(fd[t, i] - Q[s, c+1]);
      }
    }
  }
} 