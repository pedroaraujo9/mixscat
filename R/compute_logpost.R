compute_logpost = function(beta, w, pw, model_data, sd_beta) {

  M = length(pw)
  Z = model_data$data$Z
  B = model_data$spline$B
  dirichlet_param = rep(model_data$spline$dirichlet_param, M)
  n_basis = model_data$spline$n_basis
  G = model_data$dims$G

  W = fast_dummy_dense(w, M)
  W[, 1] = 1
  p = sofmax_cpp(kronecker(W, B) %*% beta)
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

  logpost = log_like + log_beta_prior +  log_w + log_pw_prior

  out = c(
    "logpost" = logpost,
    "logpenal" = log_like + log_beta_prior,
    "loglike" = log_like,
    "logpenal_active" = log_like + log_beta_active
  )

  return(out)

}
