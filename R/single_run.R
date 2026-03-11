single_run = function(M,
                      w = NULL,
                      model_data,
                      lambda,
                      init_list,
                      iters,
                      burn_in,
                      thin,
                      verbose,
                      seed = NULL) {

  if(!is.null(seed)) set.seed(seed)

  init_time = Sys.time()

  model_data$dims$M = M

  S = model_data$spline$S
  S = S * lambda
  S[1, 1] = model_data$spline$intercept_penalty
  model_data$spline$S_expand = S_expand = kronecker(diag(M), S)

  model_data$spline$sd_beta = matrix(
    1/sqrt(diag(S_expand)),
    ncol = 4,
    nrow = length(diag(S_expand)),
    byrow = F
  )

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
    update_w = TRUE

  }else{

    update_w = FALSE

  }

  beta = sample_list$beta[1,,]
  pw = sample_list$pw[1,]

  i = 1

  for(iter in 1:iters) {

    if(verbose) cat(iter, "\r")

    update = update_chain(
      beta = beta,
      w = w,
      pw = pw,
      model_data = model_data,
      update_w = update_w
    )

    beta = update$beta
    w = update$w
    w_post_prob = update$w_post_prob
    pw = update$pw

    if(iter %in% sample_list$iters_vec) {

      sample_list$beta[i,,] = beta
      sample_list$w[i, ] = w
      sample_list$pw[i, ] = pw
      sample_list$w_post_prob[i,,] = w_post_prob

      sample_list$logpost[i, ] = compute_logpost(
        beta = beta,
        w = w,
        pw = pw,
        model_data = model_data
      )

      i = i + 1

    }

  }

  out = list(
    sample_list = sample_list,
    init = init_list
  )

  end_time = Sys.time()
  out$run_time = end_time - init_time

  return(out)


}
