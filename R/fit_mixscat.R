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
                       iters,
                       thin,
                       burn_in,
                       chains,
                       seed,
                       lambda = NULL,
                       intercept_penalty = 1,
                       w_dirichlet = 1,
                       n_cores = 1,
                       verbose = TRUE,
                       init_control = list(
                         lambda_init = 1,
                         lambda_grid = seq(from = 0.01, to = 5, length.out = 50),
                         n_init = 100,
                         init_iters = 20,
                         init_burn_in = 10,
                         init_thin = 2
                       )) {

  init_time = Sys.time()

  args = list(
    z = z,
    id = id,
    time = time,
    n_basis = n_basis,
    intercept_penalty = intercept_penalty,
    w_dirichlet = w_dirichlet,
    M = M,
    chains = chains,
    iters = iters,
    burn_in = burn_in,
    thin = thin,
    lambda_init = init_control$lambda_init,
    seed = seed
  )

  model_data = create_model_data(
    z = z,
    id = id,
    time = time,
    n_basis = n_basis,
    intercept_penalty = intercept_penalty,
    w_dirichlet = w_dirichlet
  )

  if(n_cores > 1) {
    future::plan(future::multisession, workers = n_cores)
  }else{
    future::plan(future::sequential, split = TRUE)
  }

  if(verbose) cat("Starting parallel runs with", n_cores, "cores.\n")

  runs = future.apply::future_lapply(
    M, function(m){

      if(verbose) cat("Running for M =", m, "\n")

      mixscat::pipeline(
        model_data = model_data,
        M = m,
        chains = chains,
        iters = iters,
        burn_in = burn_in,
        thin = thin,
        lambda_init = init_control$lambda_init,
        lambda_grid = init_control$lambda_grid,
        n_init = init_control$n_init,
        init_iters = init_control$init_iters,
        init_burn_in = init_control$init_burn_in,
        init_thin = init_control$init_thin,
        seed = seed,
        verbose = verbose
      )

    },
    future.seed = TRUE,
    future.packages = c("Rcpp", "dplyr", "tidyr", "purrr", "posterior"),
    future.stdout = TRUE
  )

  names(runs) = paste0("M=", M)
  cluster_metrics = lapply(runs, function(run) run$metrics) |> purrr::list_rbind()
  best_lambda = lapply(runs, function(run) run$lambda_opt$best_lambda)
  clusters = lapply(runs, function(run) run$w) |> (\(x) do.call(cbind, x))()
  colnames(clusters) = paste0("M=", M)

  endtime = Sys.time()
  run_time = endtime - init_time

  out = list(
    model_data = model_data,
    args = args,
    cluster_metrics = cluster_metrics,
    best_lambda = best_lambda,
    clusters = clusters,
    fit = runs,
    run_time = run_time
  )

  return(out)

}
