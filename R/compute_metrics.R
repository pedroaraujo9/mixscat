compute_metrics = function(sample_list, lambda, model_data) {

  n = model_data$dims$n_id
  M = ncol(sample_list$pw)
  G = model_data$dims$G

  model_data$dims$M = M

  S = model_data$spline$S
  S = S * lambda
  S[1, 1] = model_data$spline$intercept_penalty
  model_data$spline$S_expand = S_expand = kronecker(diag(M), S)

  sd_beta = compute_beta_sd_matrix(model_data, lambda, M)

  beta = sample_list$beta |> compute_post_stat()
  pw = sample_list$pw |> colMeans()
  w = sample_list$w |> comp_class()

  pm_logpost = compute_logpost(
    beta = beta, w = w, pw = pw, model_data = model_data, sd_beta = sd_beta
  )

  logpenal = sample_list$logpost[, "logpenal"]
  loglike = sample_list$logpost[, "loglike"]
  avg_loglike = mean(loglike)
  avg_logpenal = mean(logpenal)

  loglike_complex = 2*(pm_logpost["loglike"] - avg_loglike)
  logpenal_complex = 2*(pm_logpost["logpenal"] - avg_logpenal)

  penal_BICM = -2 * avg_logpenal + (log(n)-1) * logpenal_complex
  BICM = -2 * avg_loglike + (log(n)-1) * loglike_complex
  ASW = ifelse(
    M == 1,
    NA,
    fpc::cluster.stats(d = model_data$data$z_dist, clustering = w)$avg.silwidth
  )

  out = data.frame(
    M = M,
    ASW = ASW,
    BICM = BICM,
    penal_BICM = penal_BICM,
    avg_loglike = avg_loglike,
    avg_logpenal = avg_logpenal,
    loglike_complex = loglike_complex,
    logpenal_complex = logpenal_complex,
    number_param = M * model_data$spline$n_basis * (G - 1) + M - 1
  )

  return(out)

}
