data{
  int<lower=1> N;                  // num obs
  int<lower=1> N_test;             // num obs in test sample
  int<lower=1> J;                  // num firms
  int<lower=1> K;                  // num coefficients
  int<lower=1, upper=J> GroupID[N]; // GroupID for obs, e.g. FirmID or Industry-YearID
  vector[N] TA;                    // Outcome total accruals
  matrix[N, K] x;                  // Predictors InAt, ChRev, PPE, etc.

  matrix[N_test, K] x_test;                        // Predictors InAt, ChRev, PPE, etc. for test sample
  int<lower=1, upper=J> GroupID_test[N_test]; // GroupID for  for test sample obs
}
parameters{
  matrix[K, J] z;                  // standard normal sampler
  cholesky_factor_corr[K] L_Omega; // hypprior coefficient correlation
  vector<lower=0>[K] tau;          // hypprior coefficient scales
  vector[K] mu_b;                  // hypprior mean coefficients
  real<lower=0> sigma;             // error-term scale
}
transformed parameters{
  matrix[J, K] b;                  // coefficient vector
  // The multivariate non-centered version:
  b = (rep_matrix(mu_b, J) + diag_pre_multiply(tau,L_Omega) * z)';
}
model{
  to_vector(z) ~ normal(0, 1);
  L_Omega ~ lkj_corr_cholesky(2);
  mu_b  ~ normal(0, 2.5);
  sigma ~ exponential(1);
  tau ~ exponential(1);
  TA ~ normal(rows_dot_product(b[GroupID] , x), sigma);
}
generated quantities {
  vector[N] y_fit_wo;
  vector[N_test] y_fit_wo_test;  // out of sample predictive distributions
  vector[N] log_lik;

  for ( i in 1:N ) {
    y_fit_wo[i] = 0;
    for (k in 1:K) {
      y_fit_wo[i] = y_fit_wo[i] + b[GroupID[i], k] * x[i, k];
    }
    log_lik[i] = normal_lpdf(TA[i] | y_fit_wo[i], sigma);
  }

  // out of sample predictive distributions
  for ( i in 1:N_test ) {
    y_fit_wo_test[i] = 0;
    for (k in 1:K) {
      y_fit_wo_test[i] = y_fit_wo_test[i] + b[GroupID_test[i], k] * x_test[i, k];
    }
  }
}

