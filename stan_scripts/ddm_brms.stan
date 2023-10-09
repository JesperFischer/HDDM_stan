// generated with brms 2.19.0
functions {
  /* Wiener diffusion log-PDF for a single response
   * Args:
   *   y: reaction time data
   *   dec: decision data (0 or 1)
   *   alpha: boundary separation parameter > 0
   *   tau: non-decision time parameter > 0
   *   beta: initial bias parameter in [0, 1]
   *   delta: drift rate parameter
   * Returns:
   *   a scalar to be added to the log posterior
   */
   real wiener_diffusion_lpdf(real y, int dec, real alpha,
                              real tau, real beta, real delta) {
     if (dec == 1) {
       return wiener_lpdf(y | alpha, tau, beta, delta);
     } else {
       return wiener_lpdf(y | alpha, tau, 1 - beta, - delta);
     }
   }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  int<lower=0,upper=1> dec[N];  // decisions
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  real min_Y = min(Y);
}
parameters {
  real Intercept;  // temporary intercept for centered predictors
  real Intercept_bs;  // temporary intercept for centered predictors
  real Intercept_ndt;  // temporary intercept for centered predictors
  real Intercept_bias;  // temporary intercept for centered predictors
}
transformed parameters {
  real lprior = 0;  // prior contributions to the log posterior
  lprior += student_t_lpdf(Intercept | 3, 0.8, 2.5);
  lprior += gamma_lpdf(Intercept_bs | 1, 1);
  lprior += uniform_lpdf(Intercept_ndt | 0, min_Y);
  lprior += beta_lpdf(Intercept_bias | 1, 1);
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] bs = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] ndt = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] bias = rep_vector(0.0, N);
    mu += Intercept;
    bs += Intercept_bs;
    ndt += Intercept_ndt;
    bias += Intercept_bias;
    for (n in 1:N) {
      target += wiener_diffusion_lpdf(Y[n] | dec[n], bs[n], ndt[n], bias[n], mu[n]);
    }
  }
  // priors including constants
  target += lprior;
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept;
  // actual population-level intercept
  real b_bs_Intercept = Intercept_bs;
  // actual population-level intercept
  real b_ndt_Intercept = Intercept_ndt;
  // actual population-level intercept
  real b_bias_Intercept = Intercept_bias;
}