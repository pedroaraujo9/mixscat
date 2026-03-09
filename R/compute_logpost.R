compute_logpost = function(beta, w, pw, model_data) {

  M = length(pw)
  Z = model_data$data$Z
  B = model_data$spline$B
  sd_beta = model_data$spline$sd_beta
  w_dirichlet = rep(model_data$spline$w_dirichlet, M)

  p = predict_prob_cpp(M, w, B, beta)
  log_like = sum(log(rowSums(Z * p)))
  log_beta_prior = sum(dnorm(beta, mean = 0, sd = sd_beta, log = T))

  if(M > 1) {

    log_pw_prior = extraDistr::ddirichlet(pw, alpha = w_dirichlet, log = T)

  }else{

    log_pw_prior = 0

  }

  logpost = log_like + log_beta_prior + log_pw_prior

  out = c(
    "logpost" = logpost,
    "logpenal" = log_like + log_beta_prior,
    "loglike" = log_like
  )

  return(out)

}
