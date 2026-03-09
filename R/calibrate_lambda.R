find_map = function(M,
                    X,
                    lambda,
                    stan_model,
                    model_data) {


  set.seed(1)
  n_basis = model_data$spline$n_basis
  intercept_sd = model_data$spline$intercept_penalty
  z = model_data$data$z
  G = model_data$dims$G

  S = model_data$spline$S
  S = lambda * S
  S[1, 1] = intercept_sd
  S = kronecker(diag(M), S)

  data_list = list(
    n = model_data$dims$n,
    p = ncol(X),
    G = G,
    z = z,
    X = X,
    sd_beta = sqrt(1/diag(S))
  )

  opt_fit = rstan::optimizing(
    object = stan_model,
    init = list(beta = matrix(0, nrow = M * n_basis, ncol = G - 1)),
    data = data_list,
    algorithm = "LBFGS",
    as_vector = FALSE,
    iter = 1000
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
    "stan", "stan_for_opt.rds", package = "bmixscat"
  )

  stan_model = readRDS(model_path)

  X = kronecker(create_dummy(w, M), model_data$spline$B)

  obj_fun = function(lambda) {
    fit = find_map(
      M = M,
      X = X,
      lambda = lambda,
      stan_model = stan_model,
      model_data = model_data
    )

    return(list(Score = fit$opt_fit$par$log_posterior, beta = fit$opt_fit$par$beta))

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

