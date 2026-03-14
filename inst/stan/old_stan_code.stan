data {
  int<lower=1> n; // number of observations
  int<lower=2> G; // number of categories
  int<lower=1> p; // number of predictors (including intercept)
  array[n] int<lower=1,upper=G> z;  // categorical outcomes
  matrix[n, p] X; // design matrix
  vector[p] sd_beta; // prior standard deviation for beta
}

parameters {
  matrix[p, G-1] beta;
}

transformed parameters {
  matrix[n, G] eta;
  eta[, 1:(G-1)] = X * beta;
  eta[, G] = rep_vector(0, n);
}

model {

  for(g in 1:(G-1)) {
    beta[, g] ~ normal(0, sd_beta);
  }

  for (i in 1:n) {
    z[i] ~ categorical_logit(eta[i]');
  }
}

generated quantities {
  real log_posterior = 0;
  for(g in 1:(G-1)) {
    log_posterior += normal_lpdf(beta[, g] | 0, sd_beta);
  }

  for (i in 1:n) {
    log_posterior += categorical_logit_lpmf(z[i] | eta[i]');
  }
}
