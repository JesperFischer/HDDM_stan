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


// based on codes/comments by Guido Biele, Joseph Burling, Andrew Ellis, and potentially others @ Stan mailing lists
data {
  int<lower=0> Nu; // of upper boundary responses
  int<lower=0> Nl; // of lower boundary responses
  array[Nu] real RTu;    // upper boundary response times
  array[Nl] real RTl;    // lower boundary response times
  real minRT;      // minimum RT of the observed data
  int<lower = 0, upper = 1> run_estimation; // a switch to evaluate the likelihood


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

  real<lower=0, upper=10> alpha;  // boundary separation
  real<lower=0, upper=1> beta;   // initial bias
  real delta;  // drift rate
  real tau_raw;  // nondecision time
}

transformed parameters{
  
   real tau = inv_logit(tau_raw) * minRT; // non-decision time at RT scale
}

model {
  
  target += normal_lpdf(alpha | 0, 3)-normal_lccdf(0 | 0, 3);
  
  target += beta_proportion_lpdf(beta | 0.5, 5);
  
  delta ~ normal(0, 3);
  
  tau_raw ~ normal(0,1);

  
  if(run_estimation==1){
    RTu ~ wiener(alpha, tau, beta, delta);
    RTl ~ wiener(alpha, tau, 1-beta, -delta);
  }
}

generated quantities {
  
  real prior_alpha; 
  real prior_beta;
  real prior_delta;
  real prior_tau;
  //real log_lik;
  
  
  array[Nu + Nl] vector[2] out;
  
  for (n in 1 : (Nu + Nl)) {
    out[n] = wiener_rng(alpha, tau, beta, delta);
  }

  
  // For log likelihood calculation

  // For posterior predictive check (Not implementeed yet)
  //vector[Nu] y_pred_upper;
  //vector[Nl] y_pred_lower;

  prior_alpha = uniform_rng(0, 5);
  prior_beta  = uniform_rng(0, 1);
  prior_delta = normal_rng(0, 2);
  prior_tau = inv_logit(normal_rng(0,1)) * minRT;
  

  //log_lik = wiener_lpdf(RTu | alpha, tau, beta, delta)+wiener_lpdf(RTl | alpha, tau, 1-beta, -delta);


}

