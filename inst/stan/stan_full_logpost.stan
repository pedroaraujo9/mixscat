functions {

  real dirichlet_multinomial_w_lpmf(array[] int w, int M, real alpha) {

    int n = size(w);
    array[M] int counts = rep_array(0, M);

    for (i in 1:n)
    counts[w[i]] += 1;

    real logp = lgamma(M * alpha) - lgamma(n + M * alpha);

    for (m in 1:M)
    logp += lgamma(counts[m] + alpha) - lgamma(alpha);

    return logp;
  }

}

data {
  int<lower=1> n;
  int<lower=1> n_id;
  int<lower=2> G;
  int<lower=1> p;
  int<lower=1> M;
  real<lower=0> dirichlet_param;
  array[n] int<lower=1, upper=G> z;
  array[n_id] int<lower=1, upper=M> w;
  matrix[n, p] X;
  vector<lower=0>[p] sd_beta;
}

transformed data {
  vector[G] intercept = rep_vector(0, G);
}

parameters {
  matrix[p, G - 1] beta;
  simplex[M] pw;
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

  pw ~ dirichlet(rep_vector(dirichlet_param, M));
  w ~ categorical(pw);
  z ~ categorical_logit_glm(X, intercept, beta_full);

}

generated quantities {
  real log_like = 0;
  log_like += categorical_logit_glm_lpmf(z | X, intercept, beta_full); // Z
  log_like += categorical_lpmf(w | pw); // W
  log_like += dirichlet_lpdf(pw | rep_vector(dirichlet_param, M)); // pw
  //log_like += dirichlet_multinomial_w_lpmf(w | M, dirichlet_param);

  // beta
  for(g in 1:(G - 1)) {
    log_like += normal_lpdf(beta[, g] | 0, sd_beta);
  }

}
