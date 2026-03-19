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

find_best_hierarchical_w = function(M, model_data, stan_model, lambda, intercept_penalty, dirichlet_param) {


  d = as.dist(model_data$data$z_dist)
  link = "ward.D"

  w_init = hclust(d, method = link) |> cutree(k = M)

  post_map = find_post_map(
    M = M,
    w = w_init,
    stan_model = stan_model,
    init_list = NULL,
    lambda = lambda,
    dirichlet_param =  dirichlet_param,
    intercept_penalty = intercept_penalty,
    model_data = model_data
  )

  init_list = list(
    w = w_init,
    beta = cbind(post_map$par$beta, 0),
    pw = post_map$par$pw
  )

  out = list(
    init_list = init_list,
    logpost = post_map$par$log_like
  )

  return(out)

}

find_init_overfitted = function(M,
                                submodels,
                                model_data,
                                lambda,
                                intercept_penalty,
                                dirichlet_param,
                                init_iters,
                                init_thin,
                                init_burn_in,
                                n_cores,
                                seed,
                                verbose = TRUE) {

  if(!is.null(seed)) set.seed(seed)

  init_time = Sys.time()

  # load stan model to find posterior mode for fixed w
  model_path = system.file(
    "stan", "stan_full_logpost.rds", package = "mixscat"
  )

  stan_model = readRDS(model_path)

  # set number of cores for parallel chains
  if(n_cores > 1) {
    future::plan(future::multisession, workers = n_cores)
  }else{
    future::plan(future::sequential, split = TRUE)
  }

  if(verbose) cat("Number of workers:", future::nbrOfWorkers(), "\n")

  # run inits
  init_runs = future.apply::future_lapply(submodels, function(m){

    if(verbose) cat("Running init for M =", m, "\n")

    # start with ward initialization
    w_init = find_best_hierarchical_w(
      M = m,
      model_data = model_data,
      intercept_penalty = intercept_penalty,
      dirichlet_param = 1,
      stan_model = stan_model,
      lambda = lambda
    )

    init_list = w_init$init_list

    # update it with a short chain
    init_run = single_run(
      M = m,
      w = NULL,
      model_data = model_data,
      lambda = lambda,
      intercept_penalty = intercept_penalty,
      dirichlet_param = 1,
      init_list = init_list,
      iters = init_iters,
      burn_in = init_burn_in,
      thin = init_thin,
      verbose = TRUE,
      add_logpost = FALSE,
      seed = NULL
    )

    # get posterior mode for the w(0) for overfitted model
    w = init_run$sample_list$w |> comp_class()

    post_map = find_post_map(
      M = M,
      w = w,
      stan_model = stan_model,
      init_list = NULL,
      lambda = lambda,
      intercept_penalty = intercept_penalty,
      dirichlet_param = dirichlet_param,
      model_data = model_data
    )

    init_list = list(
      w = w,
      beta = cbind(post_map$par$beta, 0),
      pw = post_map$par$pw |> as.numeric()
    )

    logpost_map = compute_logpost(
      beta = init_list$beta,
      w = w,
      pw = init_list$pw,
      model_data = model_data,
      dirichlet_param = dirichlet_param,
      beta_sd = create_beta_sd_matrix(
        model_data = model_data,
        lambda = lambda,
        intercept_penalty = intercept_penalty
      )
    )

    out = list(
      logpost_map = logpost_map,
      init_list = init_list
    )

    return(out)

  },

  future.seed = TRUE,
  future.packages = c("Rcpp", "dplyr", "tidyr", "purrr", "posterior"),
  future.stdout = TRUE)

  logpost_map = init_runs |> purrr::map_dbl(~.x$logpost_map["logpost"])
  idx = which.max(logpost_map)
  init_list = init_runs[[idx]]$init_list

  out = list(
    init_runs = init_runs,
    logpost_map = logpost_map,
    init_list = init_list,
    runtime = Sys.time() - init_time
  )

  return(out)

}

