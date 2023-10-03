functions {
  real fs_cdf(real t, real a) {
    if (a < 1) {
      reject("a must be >= 1, found a = ", a);
    }
    
    return erfc(inv_sqrt(2 * a * t));
  }
  
  vector make_vars(real mu) {
    real mu2 = pow(mu, 2);
    real t_tilde = 0.12 + 0.5 * exp(-mu2 / 3);
    real a = (3 + sqrt(9 + 4 * mu2)) / 6;
    real sqrtamu = sqrt((a - 1) * mu2 / a);
    real fourmu2pi = (4 * mu2 + pi() ^ 2) / 8;
    real Cf1s = sqrt(a) * exp(-sqrtamu);
    real Cf1l = pi() / (4 * fourmu2pi);
    real CF1st = Cf1s * fs_cdf(t_tilde | a);
    real F1lt = -expm1(-t_tilde * fourmu2pi);
    real F1inf = CF1st + Cf1l * (1 - F1lt);
    
    return [mu2, //.......1
            t_tilde, //.. 2
            a, //.........3
            sqrtamu, //...4
            fourmu2pi, //.5
            Cf1s, //......6
            Cf1l, //......7
            CF1st, //.....8
            F1lt, //......9
            F1inf]'; //...10
  }
  
  int acceptt_rng(real t_star, real ft, real c) {
    if (c <= 0.06385320297074884) {
      reject("c is ", c);
    }
    if (is_nan(c)) {
      reject("c is nan!");
    }
    real z = ft * uniform_rng(0, 1);
    real b = exp(-c);
    int k = 3;
    
    while (1) {
      if (z > b) {
        return 0;
      }
      b -= k * exp(-c * k ^ 2);
      if (z < b) {
        return 1;
      }
      k += 2;
      b += k * exp(-c * k ^ 2);
      k += 2;
    }
    
    return 0;
  }
  
  real sample_small_mu_rng(vector vars) {
    real t_star;
    real pi_sq = pi() ^ 2;
    
    real mu2 = vars[1];
    real a = vars[3];
    real sqrtamu = vars[4];
    real fourmu2pi = vars[5];
    real Cf1s = vars[6];
    real Cf1l = vars[7];
    real CF1st = vars[8];
    real F1lt = vars[9];
    real F1inf = vars[10];
    
    int counter_outer = 0;
    while (1) {
      real p = F1inf * uniform_rng(0, 1);
      
      if (p <= CF1st) {
        t_star = 1. / (2 * a * pow(inv_erfc(p / Cf1s), 2));
        while (0.5 * t_star <= 0.06385320297074884) {
          p = uniform_rng(0.06385320297074884, CF1st);
          t_star = 1. / (2 * a * pow(inv_erfc(p / Cf1s), 2));
        }
        real ft = exp(-1. / (2 * a * t_star) - sqrtamu + mu2 * t_star);
        if (acceptt_rng(t_star, ft, 0.5 * t_star) == 1) {
          return t_star;
        }
      } else {
        t_star = -log1p(-(p - CF1st) / Cf1l - F1lt) / fourmu2pi;
        real pisqt = pi_sq * t_star / 8;
        while (pisqt <= 0.06385320297074884) {
          p = uniform_rng(CF1st, F1inf);
          t_star = -log1p(-(p - CF1st) / Cf1l - F1lt) / fourmu2pi;
          pisqt = pi_sq * t_star / 8;
        }
        if (acceptt_rng(t_star, exp(-pisqt), pisqt) == 1) {
          return t_star;
        }
      }
    }
    return 0;
  }
  
  real inverse_gaussian_rng(real mu, real mu_sq) {
    real v = pow(std_normal_rng(), 2);
    real z = uniform_rng(0, 1);
    real x = mu + 0.5 * mu_sq * v - 0.5 * mu * sqrt(4 * mu * v + mu_sq * v ^ 2);
    if (z <= (mu / (mu + x))) {
      return x;
    } else {
      return mu_sq / x;
    }
  }
  
  real sample_large_mu_rng(vector vars) {
    real mu2 = vars[1];
    real t_tilde = vars[2];
    real a = vars[3];
    real sqrtamu = vars[4];
    real fourmu2pi = vars[5];
    real Cf1s = vars[6];
    real Cf1l = vars[7];
    real CF1st = vars[8];
    real F1lt = vars[9];
    real F1inf = vars[10];
    
    real invabsmu = inv_sqrt(mu2);
 
    if (t_tilde >= 0.63662) {
      Cf1l = -log(pi() * 0.25) - 0.5 * log(2 * pi());
      Cf1s = 0;
    } else {
      Cf1l = -pi() ^ 2 * t_tilde / 8 + (3. / 2.) * log(t_tilde) + 0.5 * inv(t_tilde);
      Cf1s = Cf1l + 0.5 * log(2 * pi()) + log(pi() * 0.25);
    }
    
    while (1) {
      real t_star = inverse_gaussian_rng(invabsmu, inv(mu2));
      if (is_nan(t_star)) {
        reject("t_star is nan! ", mu2);
      }
      real one2t = 0.5 * inv(t_star);
      if (t_star <= 2.5) {
        real expone2t = exp(Cf1s - one2t);
        if (expone2t == 0) {
          expone2t = 1e-8;
        }
        if (acceptt_rng(t_star, expone2t, one2t) == 0 || invabsmu < 0.000666) {
          return t_star;
        }
      } else {
        real expone2t = exp(-log(pi() / 4) - 0.5 * log(2 * pi()) - one2t - (3. / 2.) * log(t_star));
        if (acceptt_rng(t_star, expone2t, pi() ^ 2 * t_star / 8) == 0) {
          return t_star;
        }
      }
    }
    return 0;
  }
  
  real fast_pt_rng(real alpha, real tau, real beta, real delta) {
    real absmu = abs(delta) ;
    vector[10] vars = make_vars(absmu);
    real pt;
    
    if (absmu < 1) {
      pt = sample_small_mu_rng(vars);
    } else {
      pt = sample_large_mu_rng(vars);
    }
    
    return pt;
  }
  
  vector wiener_rng(real alpha, real tau, real beta, real delta) {
    real t = 0 ;
    real sign_delta = delta > 0 ? 1 : -1;
    real x =  beta * alpha ;
    real mu = abs(delta);
    real hit_bound;
    vector[2] out;
    int counter = 0;
    
    if (beta == 0 || beta == 1) {
      return [tau, beta]';
    }

    while (1) {
      real mutheta;
      real xlo =  x ;
      real xhi = alpha - x ;
      // lower bound is 0 
      // upper bound is alpha in stan parmeterization
      
      // symmetric case, [x - xup, x + xup]
      if (abs(xlo - xhi) < 1e-6) {
        mutheta = xhi * mu;
        real pt = fast_pt_rng(alpha, tau, beta, xhi * abs(delta));
        hit_bound = sign_delta == 1 ? inv_logit( 2 * mutheta ) : 1 - inv_logit( 2 * mutheta );
        real bound = uniform_rng(0, 1) < hit_bound ? 1 : 0;
        return [ tau + t  + ( square(xhi)  * pt), bound]';
      // x is closer to upper bound, [x - xup, x + xup]
      } else if (xlo > xhi) {
        mutheta = xhi * mu;
        t += ( square(xhi ) * fast_pt_rng(alpha, tau, beta, xhi* abs(delta)))  ;
        hit_bound = sign_delta == 1 ? inv_logit( 2 * mutheta ) : 1 - inv_logit( 2 * mutheta );
        if (uniform_rng(0, 1) < hit_bound ) {
          return [tau + t, 1]';
        }
        x -= xhi ;
      } else {
       // x is closer to lower bound, [x - xlo, x + xlo]
        mutheta = xlo * mu ;
        t +=  ( square(xlo ) * fast_pt_rng(alpha, tau,  beta, xlo* abs(delta) )) ;
        hit_bound = sign_delta == 1 ? inv_logit( 2 * mutheta ) : 1 - inv_logit( 2 * mutheta );
        if (uniform_rng(0, 1) > hit_bound) {
          out[1] = tau + t;
          out[2] = 0 ;
          break;
        }
        x += xlo ;
      }
    }
    return out;
  }
}
data {
  int N;
  array[N] int u;
  real<lower=0> alpha_in;
  real<lower=0> tau_in;
  real<lower=0, upper=1> beta_in;
  real lr_in;
  real delta_in;
} 
transformed data {
  int N_lower = 0;
  int N_upper = 0;

  array[N+1] real expect;
  array[N] real uncertain;
  array[N] real deltat;

  array[N] real rt;
  array[N] real idx;
  array[N] int idx2;

  
  expect[1] = 0.5;
  int counter = 0;
  for (n in 1 : N) {
    
    //uncertain[n] = expect[n] * (1 - expect[n]);
    uncertain[n] = expect[n] - (1 - expect[n]);
    
    //deltat[n] = delta_in * uncertain[n];

    deltat[n] = uncertain[n];
    
    expect[n+1] = expect[n]+lr_in*(u[n]-expect[n]);
    
    vector[2] tmp = wiener_rng(alpha_in, tau_in, beta_in, deltat[n]);
    
    rt[n] = tmp[1];
    idx[n] = tmp[2];
    
    if (tmp[2] == 0) {
      N_lower += 1;
    } else {
      N_upper += 1;
    }
    if (idx[n] >= 0.5) {
      idx2[n] = 1;
      }else {
      idx2[n] = 0;
    }
  } 
  array[N_lower] int id_lower;
  array[N_upper] int id_upper;
  array[2] int cnt;
  cnt[1] = 0;
  cnt[2] = 0;
  for (n in 1 : N) {
    if (idx[n] == 0) {
      cnt[1] += 1;
      id_lower[cnt[1]] = n;
    } else {
      cnt[2] += 1;
      id_upper[cnt[2]] = n;
    }
  }
  
}
parameters {
  real tau_raw; // logit of non-decision time
  real delta; // drift-rate
  real<lower=0> alpha; // boundary separation
  real<lower=0, upper=1> beta; // starting point
  real<lower=0, upper=1> lr; // learning rate
  
}
transformed parameters {
  real tau = inv_logit(tau_raw) * min(rt); // non-decision time at RT scale
  
  array[N+1] real  <lower = 0, upper = 1> expect2;
  array[N] real uncertain2;
  array[N] real deltat2;
  array[N] real rt2;
  rt2 = rt;
  array[N] real idx3;
  idx3 = idx2;
  

  expect2[1] = 0.5;

  for (n in 1 : N) {
    
    //uncertain2[n] = expect2[n] * (1 - expect2[n]);
    uncertain2[n] = expect2[n] - (1 - expect2[n]);
    
    //deltat2[n] = delta * uncertain2[n];
    deltat2[n] = uncertain2[n];
    
    
    expect2[n+1] = expect2[n]+lr*(u[n]-expect2[n]);
  
  }
}

model {
  delta ~ normal(0, 3);
  alpha ~ normal(1, 2);
  beta ~ normal(0, 2);
  lr ~ beta_proportion(0.2,3);
  
  tau_raw ~ normal(-.2, .48);


  
  for(i in 1:N){
    target += bernoulli_lpmf(idx2[i] | expect2[i]);
    
    if(idx2[i] == 1){
    target += wiener_lpdf(rt[i] | alpha, tau, beta, deltat2[i]);
    }else{
    target += wiener_lpdf(rt[i] | alpha, tau, 1 - beta, deltat2[i]);
    }

  }
}


