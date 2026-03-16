find_init_w = function(M,
                       model_data,
                       n_init,
                       init_mcmc_iters,
                       lambda,
                       seed = NULL,
                       verbose = TRUE) {

  if(!is.null(seed)) set.seed(seed)

  # load stan model to compute log-posterior mode for a fixed w
  model_path = system.file(
    "stan", "stan_full_logpost.rds", package = "mixscat"
  )

  stan_model = readRDS(model_path)

  # w initialization using Ward's method
  w_init_ward = get_w_ward(M = M, z_dist = model_data$data$z_dist)

  # posterior mode for w_init_ward
  w_init_opt = find_post_map(
    M = M,
    w = w_init_ward,
    stan_model = stan_model,
    init_list = NULL,
    lambda = lambda,
    model_data = model_data
  )

  beta_init_ward = w_init_opt$par$beta
  pw_init_ward = w_init_opt$par$pw

  # runs different inits from ward's method and pick
  # the one with the highest log-posterior mode
  init_runs = lapply(1:n_init, function(i){

    if(verbose) cat("Init run:", i, "\r")

    run = single_run(
      M = M,
      w = NULL,
      model_data = model_data,
      lambda = lambda,
      init_list = list(
        beta = cbind(beta_init_ward, 0),
        pw = pw_init_ward,
        w = w_init_ward
      ),
      iters = init_mcmc_iters,
      burn_in = 0,
      thin = 1,
      verbose = FALSE,
      seed = NULL
    )

    w = run$sample_list$w[init_mcmc_iters, ]

    opt_run = find_post_map(
      M = M,
      w = w,
      stan_model = stan_model,
      init_list = NULL,
      lambda = lambda,
      model_data = model_data
    )

    logpost = opt_run$par$log_like

    list(
      logpost = logpost,
      beta = opt_run$par$beta,
      pw = opt_run$par$pw,
      w = w
    )

  })

  idx = init_runs %>% purrr::map_dbl(~.x$logpost) %>% which.max()
  w_init = init_runs[[idx]]$w
  beta_init = init_runs[[idx]]$beta
  pw_init = init_runs[[idx]]$pw

  out = list(
    w = w_init,
    beta = cbind(beta_init, 0),
    pw = pw_init,
    ari_init = mclust::adjustedRandIndex(w_init_ward, w_init)
  )

  return(out)

}




