#' Fit mixture of categorical distribution with a spline term for
#' each component of the mixture.
#'
#' @description
#' Fits a Bayesian mixture of categorical distributions with a spline term for
#' each component of the mixture.
#'
#' @param M Integer vector specifying the number of groups.
#' @param z Numeric vector or with response variables.
#' @param id Vector of subject identifiers.
#' @param time Vector of time points.
#' @param n_basis Integer specifying the number of basis functions to use.
#' @param iters Integer specifying the number of MCMC iterations.
#' @param thin Integer specifying the thinning interval for MCMC samples.
#' @param burn_in Integer specifying the number of burn-in iterations to discard.
#' @param chains Integer specifying the number of MCMC chains to run.
#' @param seed Integer seed for random number generation to ensure reproducibility.
#' @param lambda Numeric value for the regularization parameter. If NULL, will be
#'   calibrated automatically. Default is NULL.
#' @param intercept_penalty Numeric penalty weight for the intercept term. Default is 1.
#' @param w_dirichlet Numeric hyperparameter for the Dirichlet prior on mixture weights.
#'   Default is 1.
#' @param n_cores Integer specifying the number of CPU cores to use for parallel
#'   computation. Default is 1.
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#' @param init_control List of initialization control parameters including:
#'   \itemize{
#'     \item \code{lambda_init}: Initial lambda value (default: 1)
#'     \item \code{lambda_grid}: Grid of lambda values for calibration (default: seq(0.01, 5, length.out = 50))
#'     \item \code{n_init}: Number of initializations to try (default: 100)
#'     \item \code{init_iters}: Number of iterations for initialization (default: 20)
#'     \item \code{init_burn_in}: Burn-in for initialization (default: 10)
#'     \item \code{init_thin}: Thinning for initialization (default: 2)
#'   }
#'
#' @return A list containing:
#'   \item{model_data}{Processed model data structure}
#'   \item{args}{List of arguments passed to the function}
#'   \item{cluster_metrics}{Data frame of clustering metrics for each M value}
#'   \item{best_lambda}{List of optimal lambda values for each M}
#'   \item{fit}{List of fitted model objects for each M value}
#'   \item{run_time}{Total computation time}
#'
#' @export
fit_mixscat = function(M,
                       z,
                       id,
                       time,
                       n_basis,
                       init_list,
                       iters,
                       thin,
                       burn_in,
                       chains,
                       n_cores,
                       lambda,
                       pw,
                       w,
                       intercept_penalty,
                       dirichlet_param,
                       a_lambda,
                       b_lambda,
                       seed,
                       verbose) {

  init_time = Sys.time()

  if(!is.null(seed)) set.seed(seed)

  args = list(
    iters = iters,
    thin = thin,
    burn_in = burn_in,
    chains = chains,
    lambda = lambda,
    intercept_penalty = intercept_penalty,
    dirichlet_param = dirichlet_param,
    n_cores = n_cores
  )

  model_data = create_model_data(
    z = z,
    id = id,
    time = time,
    n_basis = n_basis,
    intercept_penalty = intercept_penalty,
    dirichlet_param = dirichlet_param
  )

  if(!is.null(init_list$w)) {

    stan_model = readRDS(
      system.file("stan", "stan_full_logpost.rds", package = "mixscat")
    )

    post_map = find_post_map(
      M = M,
      w = init_list$w,
      stan_model = stan_model,
      init_list = NULL,
      lambda = lambda,
      intercept_penalty = intercept_penalty,
      dirichlet_param = dirichlet_param,
      model_data = model_data

    )

    init_list$beta = cbind(post_map$par$beta, 0)

  }

  fit = pipeline(
    M = M,
    pw = pw,
    w = w,
    chains = chains,
    n_cores = n_cores,
    model_data = model_data,
    lambda = lambda,
    intercept_penalty = intercept_penalty,
    dirichlet_param = dirichlet_param,
    init_list = init_list,
    iters = iters,
    burn_in = burn_in,
    thin = thin,
    verbose = verbose,
    seed = NULL
  )

  run_time = Sys.time() - init_time

  out = list(
    args = args,
    model_data = model_data,
    fit = fit,
    run_time = run_time
  )

  return(out)


}
