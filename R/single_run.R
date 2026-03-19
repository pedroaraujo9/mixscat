#' Single run
#'
#' @description
#' Run MCMC for one single chain
#'
#'
#' @export
single_run = function(M,
                      w = NULL,
                      model_data,
                      lambda,
                      intercept_penalty,
                      dirichlet_param,
                      init_list,
                      iters,
                      burn_in,
                      thin,
                      verbose,
                      add_logpost = TRUE,
                      seed = NULL) {

  if(!is.null(seed)) set.seed(seed)

  init_time = Sys.time()

  model_data$dims$M = M

  # define prior precision for beta
  beta_precision_matrix = create_beta_precision_matrix(
   model_data = model_data,
   lambda = lambda, intercept_penalty = intercept_penalty, M = M
  )

  beta_sd = create_beta_sd_matrix(
    model_data = model_data,
    lambda = lambda, intercept_penalty = intercept_penalty, M = M
  )

  # create list to store posterior samples
  sample_list = create_sample_list(
    M = M,
    iters = iters,
    burn_in = burn_in,
    thin = thin,
    model_data = model_data,
    init_list = init_list,
    seed = seed
  )

  # itinialise parameters
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

  # run mcmc
  for(iter in 1:iters) {

    if(verbose) cat(iter, "\r")

    if(iter %in% sample_list$iters_vec) {

      sample_list$beta[i,,] = beta
      sample_list$w[i, ] = w
      sample_list$pw[i, ] = pw
      sample_list$w_post_prob[i,,] = w_post_prob

      if(add_logpost == TRUE) {
        sample_list$logpost[i, ] = compute_logpost(
          beta = beta,
          w = w,
          pw = pw,
          dirichlet_param = dirichlet_param,
          model_data = model_data,
          beta_sd = beta_sd
        )
      }

      i = i + 1

    }

    update = update_chain(
      beta = beta,
      w = w,
      pw = pw,
      model_data = model_data,
      dirichlet_param = dirichlet_param,
      update_w_iter = update_w_iter,
      beta_precision_matrix = beta_precision_matrix
    )

    beta = update$beta
    w = update$w
    w_post_prob = update$w_post_prob
    pw = update$pw

  }

  # save results
  out = list(
    sample_list = sample_list,
    init = init_list
  )

  end_time = Sys.time()
  out$run_time = end_time - init_time

  return(out)


}
