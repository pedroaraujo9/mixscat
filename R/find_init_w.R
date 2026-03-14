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
                         verbose = TRUE
                       )) {

  set.seed(seed)

  model_path = system.file(
    "stan", "stan_full_logpost.rds", package = "mixscat"
  )

  stan_model = readRDS(model_path)

  G = model_data$dims$G
  n_basis = model_data$spline$n_basis
  n_id = model_data$dims$n_id
  id_unique = model_data$data$id_unique
  S = model_data$spline$S
  S = init_control$lambda_init * S
  S[1, 1] = model_data$spline$intercept_penalty * S[1, 1]
  model_data$spline$S_expand = kronecker(diag(M), S)

  # run different initializations
  init_runs = lapply(1:init_control$n_init, function(i){

    if(init_control$verbose) cat(paste0("init iter: ", i, "\r"))

    beta = matrix(
      rnorm(M * n_basis * (G - 1), sd = 0.01),
      nrow = M * n_basis,
      ncol = G - 1
    )

    beta = cbind(beta, 0)

    pw = rep(1/M, M)
    w = sample(1:M, size = n_id, replace = TRUE)
    names(w) = id_unique
    model_data$dims$M = M

    for(i in 1:init_iters) {

      update = update_chain(
        beta = beta,
        w = w,
        pw = pw,
        model_data = model_data,
        update_w_iter = TRUE,
        temperature = 1
      )

      beta = update$beta
      w = update$w
      pw = update$pw

    }

    return(w)

  })

  ws = init_runs %>% do.call(rbind, .)

  ari_matrix = lapply(1:nrow(ws), function(j){
    lapply(1:nrow(ws), function(i){
      ari(ws[i, ], ws[j, ])
    }) %>% do.call(c, .)
  }) %>% do.call(rbind, .)

  opt = find_post_map(
    M = M,
    w = w,
    stan_model = stan_model,
    init_list = NULL,
    lambda = init_control$lambda_init,
    model_data = model_data
  )



  # select w with the highest mode
  logpost = purrr::map_dbl(init_runs, ~{.x$logpost})
  max_logpost = which.max(logpost)

  out = init_runs[[max_logpost]]

  return(out)

}

find_init_w2 = function(M,
                        model_data,
                        seed,
                        init_control = list(
                          lambda_init = 1,
                          n_init = 30,
                          init_iters = 30,
                          init_burn_in = 15,
                          init_thin = 2,
                          init_final_run = 100,
                          verbose = TRUE
                        )) {

  set.seed(seed)

  # run different initializations
  init_runs = lapply(1:init_control$n_init, function(i){

    if(init_control$verbose) cat(paste0("init iter: ", i, "\r"))

    run1 = single_run(
      M = M,
      w = NULL,
      model_data = model_data,
      lambda = init_control$lambda_init,
      init_list = NULL,
      iters = init_control$init_iters,
      burn_in = init_control$init_burn_in,
      thin = init_control$init_thin,
      verbose = FALSE,
      seed = NULL
    )

    w = run1$sample_list$w[nrow(run1$sample_list$w), ]

    run2 = single_run(
      M = M,
      w = w,
      model_data = model_data,
      lambda = init_control$lambda_init,
      init_list = NULL,
      iters = init_control$init_final_run,
      burn_in = init_control$init_final_run/2,
      thin = 1,
      verbose = FALSE,
      seed = NULL
    )

    logpost = run2$sample_list$logpost[, 1] |> mean()

    out = list(
      logpost = logpost,
      w = w
    )

    return(out)

  })

  # select w with the highest mode
  logpost = purrr::map_dbl(init_runs, ~{.x$logpost})
  max_logpost = which.max(logpost)

  out = init_runs[[max_logpost]]

  return(out)

}


find_init_overfitted = function(M,
                                model_data,
                                seed,
                                init_control = list(
                                  lambda_init = 1,
                                  n_init = 30,
                                  init_iters = 30,
                                  init_burn_in = 15,
                                  init_thin = 2,
                                  init_final_run = 100,
                                  verbose = TRUE
                                )) {

  set.seed(seed)

  model_path = system.file(
    "stan", "stan_full_logpost.rds", package = "mixscat"
  )

  stan_model = readRDS(model_path)

  z_dist = model_data$data$z_dist
  cl = hclust(as.dist(z_dist), method = "ward.D2")

  # run different initializations
  init_runs = lapply(1:M, function(m){

    cat(paste0("init iter: ", m, "\r"))

    w_init = cutree(cl, k = m)
    init_list = NULL

    opt = find_post_map(
      M = M,
      w = w_init,
      stan_model = stan_model,
      init_list = init_list,
      lambda = init_control$lambda_init,
      model_data = model_data
    )

    out = list(
      logpost = opt$par$log_like,
      beta = opt$par$beta,
      pw = opt$par$pw,
      w = w_init,
      M = m
    )

    return(out)

  })

  # select w with the highest mode
  logpost = purrr::map_dbl(init_runs, ~{.x$logpost})
  max_logpost = which.max(logpost)
  out = init_runs[[max_logpost]]

  return(out)

}





