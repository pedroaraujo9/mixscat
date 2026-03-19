update_w = function(beta, z, pw, model_data) {

  G = model_data$dims$G
  M = model_data$dims$M
  n_id = model_data$dims$n_id
  n_basis = model_data$spline$n_basis

  id_unique = model_data$data$id_unique
  id = model_data$data$id
  id = factor(id, levels = id_unique)
  time_seq = model_data$data$time_seq
  B = model_data$spline$B
  Z = model_data$data$Z

  ll = matrix(nrow = n_id, ncol = M)


  for(m in 1:M) {

    idx = ((m-1)*(n_basis) + 1):(m*(n_basis))
    beta_group = beta[idx, , drop = FALSE]
    prob_group = compute_prob_group(B, beta_group, time_seq-1)

    log_pz = log(rowSums(Z * prob_group))

    ll[, m] = fast_aggregate_sum(log_pz, id)[, 1] + log(pw[m])

  }

  ll = ll - matrix(
    mclust::logsumexp(ll), nrow = n_id, ncol = M, byrow = F
  )

  ll = ll # / temperature

  w_post_prob = exp(ll)

  w = as.integer(extraDistr::rcatlp(n = n_id, log_prob = ll) + 1)
  names(w) = id_unique

  out = list(w = w, w_post_prob = w_post_prob)

  return(out)

}

update_w2 = function(beta, z, pw, model_data, temperature = 1) {

  G = model_data$dims$G
  M = model_data$dims$M
  n_id = model_data$dims$n_id
  n_basis = model_data$spline$n_basis
  alpha = model_data$spline$dirichlet_param[1]

  id_unique = model_data$data$id_unique
  id = model_data$data$id
  id = factor(id, levels = id_unique)

  time_seq = model_data$data$time_seq
  B = model_data$spline$B
  Z = model_data$data$Z

  counts = tabulate(w, nbins = M)

  w_post_prob = matrix(nrow = n_id, ncol = M)

  for(i in 1:n_id) {

    ll = numeric(M)

    counts[w[i]] = counts[w[i]] - 1

    for(m in 1:M) {

      idx = ((m-1)*(n_basis) + 1):(m*(n_basis))
      beta_group = beta[idx, , drop = FALSE]

      prob_group = compute_prob_group(B, beta_group, time_seq-1)

      log_pz = log(rowSums(Z * prob_group))

      ll[m] = fast_aggregate_sum(log_pz, id)[i,1] + log(counts[m] + alpha)

    }

    ll = ll - mclust::logsumexp(ll)
    ll = ll / temperature
    w_post_prob[i, ] = exp(ll)

    w[i] = extraDistr::rcatlp(1, ll) + 1

    counts[w[i]] = counts[w[i]] + 1
  }

  names(w) = id_unique

  out = list(w = w, w_post_prob = w_post_prob)

  return(out)

}
