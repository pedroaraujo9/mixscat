log_prior_w = function(w, M, alpha) {

  n <- length(w)

  # counts per category
  counts <- tabulate(w, nbins = M)

  # log Dirichlet–multinomial
  logp <- lgamma(M * alpha) - lgamma(n + M * alpha) +
    sum(lgamma(counts + alpha) - lgamma(alpha))

  return(logp)
}

compute_logpost = function(beta, w, pw, model_data, sd_beta) {

  M = length(pw)
  Z = model_data$data$Z
  B = model_data$spline$B
  dirichlet_param = rep(model_data$spline$dirichlet_param, M)
  n_basis = model_data$spline$n_basis
  G = model_data$dims$G

  p = predict_prob_cpp(M, w, B, beta)
  log_like = sum(log(rowSums(Z * p)))
  log_beta_prior = sum(dnorm(beta[, -G], mean = 0, sd = sd_beta, log = T))

  empty_clusters = which(table(factor(w, levels = 1:M)) == 0)

  if(length(empty_clusters) == 0) {

    beta_active = beta

  }else{

    for(m in empty_clusters) {
      idx = ((m-1)*(n_basis) + 1):(m*(n_basis))
      beta_active = beta
      beta_active[idx, ] = NA
    }

  }

  log_beta_active = sum(dnorm(beta_active[, G], mean = 0, sd = sd_beta, log = T), na.rm = T)

  if(M > 1) {

    log_pw_prior = extraDistr::ddirichlet(pw, alpha = rep(dirichlet_param, M), log = T)

  }else{

    log_pw_prior = 0

  }

  log_w = sum(log(pw[w]))
  log_w_int = log_prior_w(w, M, model_data$spline$dirichlet_param[1])

  logpost = log_like + log_beta_prior +  log_w + log_pw_prior
  #logpost2 = log_like + log_beta_active + log_w_int

  out = c(
    "logpost" = logpost,
    #"logpost2" = logpost2,
    "logpenal" = log_like + log_beta_prior,
    "loglike" = log_like,
    "logpenal_active" = log_like + log_beta_active
  )

  return(out)

}
