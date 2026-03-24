data {
  int<lower=1> n;
  int<lower=2> G;
  int<lower=1> p;
  array[n] int<lower=1, upper=G> z;
  matrix[n, p] X;
  vector<lower=0>[p] sd_beta;
}

transformed data {
  vector[G] intercept = rep_vector(0, G);
}

parameters {
  matrix[p, G - 1] beta;
}

transformed parameters {
  matrix[p, G] beta_full;
  beta_full[, 1:(G - 1)] = beta;
  beta_full[, G] = rep_vector(0, p);
}

model {
  for (g in 1:(G - 1)) {
    beta[, g] ~ normal(0, sd_beta);
  }

  z ~ categorical_logit_glm(X, intercept, beta_full);
}

generated quantities {
  real log_like = categorical_logit_glm_lpmf(z | X, intercept, beta_full);
  real log_prior = 0;

  for(g in 1:(G - 1)) {
    log_prior += normal_lpdf(beta[, g] | 0, sd_beta);
  }

  real log_like_penal = log_like + log_prior;

}
