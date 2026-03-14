find_map = function(M,
                    X,
                    lambda,
                    model_data) {

  n_basis = model_data$spline$n_basis
  intercept_sd = model_data$spline$intercept_penalty
  z = model_data$data$z
  G = model_data$dims$G

  sd_beta = compute_beta_sd_matrix(model_data, lambda, M)
  sd_beta = sd_beta[, 1]

  data_list = list(
    n = model_data$dims$n,
    p = ncol(X),
    G = G,
    z = z,
    X = X,
    sd_beta = sd_beta
  )

  model_path = system.file(
    "stan", "stan_for_opt.rds", package = "mixscat"
  )

  stan_model = readRDS(model_path)

  opt_fit = rstan::optimizing(
    object = stan_model,
    init = list(beta = matrix(0, nrow = M * n_basis, ncol = G - 1)),
    data = data_list,
    algorithm = "LBFGS",
    as_vector = FALSE,
    iter = 1000,
    seed = 1
  )

  out = list(
    opt_fit = opt_fit,
    data_list = data_list
  )

  return(out)

}

calibrate_lambda = function(w,
                            M,
                            model_data,
                            lambda_grid) {



  X = kronecker(create_dummy(w, M), model_data$spline$B)
  n_basis = model_data$spline$n_basis

  obj_fun = function(lambda) {


    fit = find_map(
      M = M,
      X = X,
      lambda = lambda,
      model_data = model_data
    )

    score = fit$opt_fit$par$log_like
    beta = fit$opt_fit$par$beta

    return(list(Score = score, beta = beta))

  }

  log_posterior = purrr::map_dbl(lambda_grid, ~{
    obj_fun(.x)$Score
  })

  best_lambda = lambda_grid[which.max(log_posterior)]

  out = list(
    fit = data.frame(lambda = lambda_grid, log_posterior = log_posterior),
    best_lambda = best_lambda
  )

  out$lambda_plot = (
    ggplot2::ggplot(out$fit, ggplot2::aes(x=lambda, y=log_posterior)) +
      ggplot2::geom_point() +
      ggplot2::geom_line()
  )

  return(out)

}


find_post_map = function(M,
                         w,
                         stan_model,
                         init_list,
                         lambda,
                         model_data) {

  n_basis = model_data$spline$n_basis
  intercept_sd = model_data$spline$intercept_penalty
  z = model_data$data$z
  G = model_data$dims$G

  sd_beta = compute_beta_sd_matrix(model_data, lambda, M)
  sd_beta = sd_beta[, 1]

  X = kronecker(create_dummy(w, M), model_data$spline$B)
  n_basis = model_data$spline$n_basis

  data_list = list(
    n = model_data$dims$n,
    n_id = model_data$dims$n_id,
    w = w,
    M = M,
    p = ncol(X),
    G = G,
    dirichlet_param = model_data$spline$dirichlet_param,
    z = z,
    X = X,
    sd_beta = sd_beta
  )

  if(is.null(init_list)) {

    init_list = list(
      beta = matrix(0, nrow = M * n_basis, ncol = G - 1),
      pw = rep(1/M, M)
    )

  }

  opt_fit = rstan::optimizing(
    object = stan_model,
    init = init_list,
    data = data_list,
    algorithm = "LBFGS",
    as_vector = FALSE,
    iter = 20000,
    seed = 1
  )

  return(opt_fit)

}


calibrate_lambda2 = function(w,
                             M,
                             model_data,
                             lambda_grid) {


  X = kronecker(create_dummy(w, M), model_data$spline$B)
  n_basis = model_data$spline$n_basis

  model_path = system.file(
    "stan", "stan_for_opt.rds", package = "mixscat"
  )

  stan_model = readRDS(model_path)

  obj_fun = function(lambda) {

    find_map(
      M = M,
      w = w,
      stan_model = stan_model,
      init_list = NULL,
      lambda = lambda,
      model_data = model_data
    )

    score = fit$opt_fit$par$log_like
    beta = fit$opt_fit$par$beta

    beta_sd = compute_beta_sd_matrix(model_data, lambda, Me)
    score = score + sum(dnorm(beta, sd = beta_sd, log = TRUE), na.rm = T)

    return(list(Score = score, beta = fit$opt_fit$par$beta))

  }

  log_posterior = purrr::map_dbl(lambda_grid, ~{
    obj_fun(.x)$Score
  })

  best_lambda = lambda_grid[which.max(log_posterior)]

  out = list(
    fit = data.frame(lambda = lambda_grid, log_posterior = log_posterior),
    best_lambda = best_lambda
  )

  out$lambda_plot = (
    ggplot2::ggplot(out$fit, ggplot2::aes(x=lambda, y=log_posterior)) +
      ggplot2::geom_point() +
      ggplot2::geom_line()
  )

  return(out)

}














