log_prior_w = function(w, M, alpha) {

  n <- length(w)

  # counts per category
  counts <- tabulate(w, nbins = M)

  # log Dirichlet–multinomial
  logp <- lgamma(M * alpha) - lgamma(n + M * alpha) +
    sum(lgamma(counts + alpha) - lgamma(alpha))

  return(logp)
}

compute_logpost = function(beta, w, pw, model_data, beta_sd, dirichlet_param) {

  M = length(pw)
  Z = model_data$data$Z
  B = model_data$spline$B
  n_basis = model_data$spline$n_basis
  G = model_data$dims$G

  p = predict_prob_cpp(M, w, B, beta)
  log_like = sum(log(rowSums(Z * p)))
  log_beta_prior = sum(dnorm(beta[, -G], mean = 0, sd = beta_sd, log = T))

  if(M > 1) {

    log_pw_prior = extraDistr::ddirichlet(pw, alpha = rep(dirichlet_param, M), log = T)

  }else{

    log_pw_prior = 0

  }

  log_w = sum(log(pw[w]))
  logpost = log_like + log_beta_prior +  log_w + log_pw_prior

  out = c(
    "logpost" = logpost,
    "loglike" = log_like,
    "logpenal" = log_like + log_beta_prior,
    "beta_logprior" = log_beta_prior
  )

  return(out)

}
