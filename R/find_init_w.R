find_init_w = function(M,
                       model_data,
                       seed,
                       init_control = list(
                         lambda_init = 1,
                         n_init = 30,
                         init_iters = 30,
                         init_burn_in = 15,
                         init_thin = 2,
                         init_final_run = 100,
                         verbose = FALSE
                       )) {

  set.seed(seed)

  # run different initializations
  init_runs = lapply(1:init_control$n_init, function(i){

    if(verbose) cat(paste0("init iter: ", i, "\r"))

    # short chain
    run = single_run(
      M = M,
      w = NULL,
      model_data = model_data,
      lambda = init_control$lambda_init,
      init_list = NULL,
      iters = init_control$init_iters,
      burn_in = init_control$init_burn_in,
      thin = init_control$init_thin,
      seed = NULL,
      verbose = FALSE
    )

    w = run$sample_list$w |> comp_class()

    # compute log-posterior at the mode
    opt_w = find_map(
      M = M,
      X = kronecker(create_dummy(w, M), model_data$spline$B),
      lambda = init_control$lambda_init,
      stan_model = readRDS(system.file(
        "stan", "stan_for_opt.rds", package = "mixscat"
      )),
      model_data = model_data
    )

    run$logpost = opt_w$opt_fit$par$log_posterior
    run$beta = cbind(opt_w$opt_fit$par$beta, 0)
    run$w = w

    return(run)

  })

  # select w with the highest mode
  logpost = purrr::map_dbl(init_runs, ~{.x$logpost})
  max_logpost = which.max(logpost)

  w0 = init_runs[[max_logpost]]$w
  beta0 = init_runs[[max_logpost]]$beta
  pw0 = w0 |> factor(levels = 1:M) |> table() |> prop.table()
  pw0[pw0 == 0] = 1e-300

  out = list(
    w0 = w0,
    beta0 = beta0,
    pw0 = pw0
  )

  return(out)

}











