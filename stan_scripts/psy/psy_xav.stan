data {
  
  int<lower=0> trials;
  
  
  array[trials] real x;
  array[trials] int r;
}

parameters {


  real psy_alpha;
  real<lower=0> psy_beta;
  real<lower=0, upper = 0.5> psy_lambda;


  
}

transformed parameters{
    
  array[trials] real phi;
  
   for(i in 1:trials){
     phi[i] = psy_lambda + (1 - 2 * psy_lambda) *  (0.5+0.5*erf((x[i]-psy_alpha)/(psy_beta*sqrt(2))));
     
  }
     
     
}


model {
  
  target += normal_lpdf(psy_alpha | 0,20);
  target += normal_lpdf(psy_beta | 1, 20)-normal_lccdf(1 | 1, 20);
  
  target += beta_proportion_lpdf(psy_lambda | 0.05, 5) / 2;
  
  
  
    for(i in 1:trials){
      
      target += bernoulli_lpmf(r[i] | phi[i]); 

  }
}

