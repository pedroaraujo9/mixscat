find_map = function(M,
                    w,
                    lambda,
                    stan_model,
                    model_data) {

  n_basis = model_data$spline$n_basis
  intercept_sd = model_data$spline$intercept_penalty
  z = model_data$data$z
  G = model_data$dims$G
  X = kronecker(create_dummy(w, M), model_data$spline$B)

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

  opt_fit = rstan::optimizing(
    object = stan_model,
    init = list(beta = matrix(0, nrow = M * n_basis, ncol = G - 1)),
    data = data_list,
    algorithm = "LBFGS",
    as_vector = FALSE,
    iter = 3000
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


  model_path = system.file(
    "stan", "stan_for_opt.rds", package = "mixscat"
  )

  stan_model = readRDS(model_path)

  obj_fun = function(lambda) {

    fit = find_map(
      M = M,
      w = w,
      lambda = lambda,
      stan_model = stan_model,
      model_data = model_data
    )

    score = fit$opt_fit$par$log_like

    return(list(Score = score))

  }

  log_posterior = purrr::map_dbl(lambda_grid, ~{
    obj_fun(.x)$Score
  })

  best_lambda = lambda_grid[which.max(log_posterior)]

  out = list(
    fit = data.frame(lambda = lambda_grid, log_posterior = log_posterior),
    best_lambda = best_lambda
  )

  # out$lambda_plot = (
  #   ggplot2::ggplot(out$fit, ggplot2::aes(x=lambda, y=log_posterior)) +
  #     ggplot2::geom_point() +
  #     ggplot2::geom_line()
  # )

  return(out)

}


find_post_map = function(M,
                         w,
                         stan_model,
                         init_list,
                         lambda,
                         intercept_penalty,
                         dirichlet_param,
                         model_data) {


  n_basis = model_data$spline$n_basis
  z = model_data$data$z
  G = model_data$dims$G

  sd_beta = create_beta_sd_matrix(
    model_data = model_data,
    lambda = lambda,
    intercept_penalty = intercept_penalty,
    M = M
  )

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
    dirichlet_param = dirichlet_param,
    z = z,
    X = X,
    sd_beta = sd_beta
  )

  if(is.null(init_list)) {

    init_list = list(
      beta = matrix(0, nrow = M * n_basis, ncol = G - 1),
      pw = array(rep(1/M, M))
    )

  }

  opt_fit = rstan::optimizing(
    object = stan_model,
    init = init_list,
    data = data_list,
    algorithm = "LBFGS",
    as_vector = FALSE,
    iter = 2000,
    seed = 100
  )

  return(opt_fit)

}

calibrate_lambda_post_map = function(w,
                                     M,
                                     model_data,
                                     lambda_grid) {


  X = kronecker(create_dummy(w, M), model_data$spline$B)
  n_basis = model_data$spline$n_basis

  model_path = system.file(
    "stan", "stan_full_logpost.rds", package = "mixscat"
  )

  stan_model = readRDS(model_path)

  obj_fun = function(lambda) {

    fit = find_post_map(
      M = M,
      w = w,
      stan_model = stan_model,
      init_list = NULL,
      lambda = lambda,
      model_data = model_data
    )

    score = fit$par$log_like

    return(list(Score = score))

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


# calibrate_lambda = function(w,
#                             M,
#                             model_data,
#                             lambda_grid) {
#
#
#   X = kronecker(create_dummy(w, M), model_data$spline$B)
#   n_basis = model_data$spline$n_basis
#
#   model_path = system.file(
#     "stan", "stan_for_opt.rds", package = "mixscat"
#   )
#
#   stan_model = readRDS(model_path)
#
#   obj_fun = function(lambda) {
#
#     fit = find_map(
#       M = M,
#       w = w,
#       stan_model = stan_model,
#       init_list = NULL,
#       lambda = lambda,
#       model_data = model_data,
#       seed = 1
#     )
#
#     score = fit$par$log_like
#
#     return(list(Score = score))
#
#   }
#
#   log_posterior = purrr::map_dbl(lambda_grid, ~{
#     obj_fun(.x)$Score
#   })
#
#   best_lambda = lambda_grid[which.max(log_posterior)]
#
#   out = list(
#     fit = data.frame(lambda = lambda_grid, log_posterior = log_posterior),
#     best_lambda = best_lambda
#   )
#
#   out$lambda_plot = (
#     ggplot2::ggplot(out$fit, ggplot2::aes(x=lambda, y=log_posterior)) +
#       ggplot2::geom_point() +
#       ggplot2::geom_line()
#   )
#
#   return(out)
#
# }
#
#
#
#
#
#
#
#
#





