data {
  int<lower=0> N;
  int<lower=0> n_points;
  array[N] int id;
  vector[N] sima;
  vector[N] simb;
  
}

transformed data{
  int n_param = 4;
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {

  // hierarchical group level means 
  vector [n_param] gm;
  // hierarchical group level deviations
  vector<lower = 0>[n_param]  tau_u;
  
  matrix[n_param, n_points] z_expo;
  // for the cholesky decomposition
  cholesky_factor_corr[n_param] L_u;
  
}

transformed parameters{

  vector[n_points] mua;
  vector[n_points] mub;

  vector<lower=0>[n_points] sigmaa;
  vector<lower=0>[n_points] sigmab;

  matrix[n_points, n_param] indi = (diag_pre_multiply(tau_u, L_u) * z_expo)';
  
  mua = gm[1]+(indi[,1]);
  mub = gm[2]+(indi[,2]);
  sigmaa = exp(gm[3]+(indi[,3]));
  sigmab = exp(gm[4]+(indi[,4]));
  
}


model {
  
  target += normal_lpdf(gm | 0, 5);
  target += lognormal_lpdf(tau_u | 0, 2);
  
  target += std_normal_lpdf(to_vector(z_expo));

  target += lkj_corr_cholesky_lpdf(L_u | 2);
  
  
  for(n in 1:N){
    sima[n] ~ normal(mua[id[n]], sigmaa[id[n]]);
    simb[n] ~ normal(mub[id[n]], sigmab[id[n]]);
    
  }
}

generated quantities{
  
  real r = (sum(mua-mean(mua)*mub-mean(mub)))/(sqrt((sum((mua-mean(mua))^2)*sum((mub-mean(mub))^2))));
  
  matrix[n_param,n_param] corr = L_u * L_u';
  
}


