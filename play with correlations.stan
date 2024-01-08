
data {
  int<lower=0> N;
  int<lower=0> n_points;
  array[N] int id;
  vector[N] sima;
  vector[N] simb;
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[n_points] mua;
  vector[n_points] mub;
  
  vector<lower=0>[n_points] sigmaa;
  vector<lower=0>[n_points] sigmab;
  
}


model {
  
  target += normal_lpdf(mua | 0, 5);
  target += normal_lpdf(mub | 0, 5);
  target += lognormal_lpdf(sigmab | 0, 2);
  target += lognormal_lpdf(sigmaa | 0, 2);
  
  for(n in 1:N){
    sima[n] ~ normal(mua[id[n]], sigmaa[id[n]]);
    simb[n] ~ normal(mub[id[n]], sigmab[id[n]]);
    
  }
}

generated quantities{
  

  real correlation;
  real sum_squared_diff_alpha;
  real sum_product_diff_alpha;

  array[n_points] real diff_s1_alpha;
  array[n_points] real diff_s2_alpha;
  array[n_points] real sq_diff_s1_alpha;
  array[n_points] real sq_diff_s2_alpha;
  array[n_points] real product_diff_alpha;
 
 
   for (s in 1:n_points){
    diff_s1_alpha[s] = mean(mua)-mua[s];
    diff_s2_alpha[s] = mean(mub)-mub[s];
    
    product_diff_alpha[s] = diff_s1_alpha[s] * diff_s2_alpha[s];
    
    sq_diff_s1_alpha[s] = (mean(mua)-mua[s])^2;
    sq_diff_s2_alpha[s] = (mean(mub)-mub[s])^2;
  }
  
  sum_squared_diff_alpha = sum(sq_diff_s1_alpha) * sum(sq_diff_s2_alpha);
  
  sum_product_diff_alpha = sum(product_diff_alpha);

  correlation = sum_product_diff_alpha / sqrt(sum_squared_diff_alpha);

  
}

