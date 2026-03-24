#' Single run
#'
#' @description
#' Run MCMC for one single chain
#'
#'
#' @export
single_run = function(M,
                      w = NULL,
                      pw = NULL,
                      model_data,
                      lambda = NULL,
                      a_lambda = 1,
                      b_lambda = 1,
                      intercept_penalty,
                      dirichlet_param,
                      init_list,
                      iters,
                      burn_in,
                      thin,
                      verbose,
                      seed = NULL) {

  if(!is.null(seed)) set.seed(seed)

  init_time = Sys.time()

  model_data$dims$M = M
  G = model_data$dims$G
  n_basis = model_data$spline$n_basis
  s = diag(model_data$spline$S)

  model_data$temp = list()
  model_data$temp$s = s
  model_data$temp$s_matrix = matrix(s, nrow = n_basis, ncol = G, byrow = F)
  model_data$temp$basis_idx_list = lapply(1:M, function(m){
    idx = ((m-1)*(n_basis) + 1):(m*(n_basis))
    idx[-1]
  })

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
    fixed_w = FALSE

  }else{

    fixed_w = TRUE
    sample_list$w = NULL
    sample_list$w_post_prob = NULL

  }

  if(is.null(pw)) {

    pw = sample_list$pw[1, ]
    fixed_pw = FALSE

  }else{

    fixed_pw = TRUE
    sample_list$pw = NULL

  }

  if(is.null(lambda)) {

    lambda = sample_list$lambda[1, ]
    fixed_lambda = FALSE
    beta_precision_matrix = NULL
    beta_sd = NULL

  }else{

    fixed_lambda = TRUE
    sample_list$lambda = NULL
    lambda = rep(lambda, M)

    beta_precision_matrix = create_beta_precision_matrix(
      model_data = model_data,
      lambda = lambda, intercept_penalty = intercept_penalty, 
      M = M
    )

    beta_sd = create_beta_sd_matrix(
      model_data = model_data,
      lambda = lambda, intercept_penalty = intercept_penalty, 
      M = M
    )

  }

  beta = sample_list$beta[1,,]
  w_post_prob = sample_list$w_post_prob[1,,]
  i = 1

  # run mcmc
  for(iter in 1:iters) {

    if(verbose) cat(iter, "\r")

    if(iter %in% sample_list$iters_vec) {

      sample_list$beta[i,,] = beta

      if(fixed_w == FALSE) {
        sample_list$w[i, ] = w
        sample_list$w_post_prob[i,,] = w_post_prob
      }

      if(fixed_pw == FALSE) sample_list$pw[i, ] = pw
      if(fixed_lambda == FALSE) sample_list$lambda[i, ] = lambda
      i = i + 1

    }

    update = update_chain(
      beta = beta,
      w = w,
      pw = pw,
      model_data = model_data,
      dirichlet_param = dirichlet_param,
      lambda = lambda,
      a_lambda = a_lambda,
      b_lambda = b_lambda,
      fixed_lambda = fixed_lambda,
      fixed_w = fixed_w,
      fixed_pw = fixed_pw,
      beta_precision_matrix = beta_precision_matrix
    )

    beta = update$beta
    w = update$w
    w_post_prob = update$w_post_prob
    pw = update$pw
    if(fixed_lambda == FALSE) lambda = update$lambda

  }

  # save results
  out = list(
    sample_list = sample_list,
    init = init_list,
    beta_precision_matrix = beta_precision_matrix,
    beta_sd = beta_sd
  )

  end_time = Sys.time()
  out$run_time = end_time - init_time

  return(out)


}
