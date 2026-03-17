single_run = function(M,
                      w = NULL,
                      model_data,
                      lambda,
                      init_list,
                      iters,
                      burn_in,
                      thin,
                      verbose,
                      temperature_vec = rep(1, iters),
                      seed = NULL) {

  if(!is.null(seed)) set.seed(seed)

  init_time = Sys.time()

  model_data$dims$M = M
  S = model_data$spline$S
  S = lambda * S
  S[1, 1] = model_data$spline$intercept_penalty
  model_data$spline$S_expand = kronecker(diag(M), S)
  sd_beta = compute_beta_sd_matrix(model_data, lambda, M)

  sample_list = create_sample_list(
    M = M,
    iters = iters,
    burn_in = burn_in,
    thin = thin,
    model_data = model_data,
    init_list = init_list,
    seed = seed
  )

  if(is.null(w)) {

    w = sample_list$w[1, ]
    update_w_iter = TRUE

  }else{

    update_w_iter = FALSE

  }

  beta = sample_list$beta[1,,]
  pw = sample_list$pw[1,]
  w_post_prob = sample_list$w_post_prob[1,,]

  i = 1

  for(iter in 1:iters) {

    if(verbose) cat(iter, "\r")

    if(iter %in% sample_list$iters_vec) {

      sample_list$beta[i,,] = beta
      sample_list$w[i, ] = w
      sample_list$pw[i, ] = pw
      sample_list$w_post_prob[i,,] = w_post_prob

      sample_list$logpost[i, ] = compute_logpost(
        beta = beta,
        w = w,
        pw = pw,
        model_data = model_data,
        sd_beta = sd_beta
      )

      i = i + 1

    }

    update = update_chain(
      beta = beta,
      w = w,
      pw = pw,
      model_data = model_data,
      update_w_iter = update_w_iter,
      temperature = temperature_vec[iter]
    )

    beta = update$beta
    w = update$w
    w_post_prob = update$w_post_prob
    pw = update$pw

  }

  out = list(
    sample_list = sample_list,
    init = init_list
  )

  end_time = Sys.time()
  out$run_time = end_time - init_time

  return(out)


}
