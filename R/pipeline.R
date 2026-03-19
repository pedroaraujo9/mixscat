#' Pipeline for fitting mixture of spline-based categorical distributions
#'
#' @description
#' Complete pipeline for fitting a Bayesian mixture model with spline-based
#' categorical distributions. This function orchestrates the entire fitting process,
#' including initialization, lambda calibration (if needed), running multiple MCMC chains,
#' convergence checking, and computing model evaluation metrics.
#'
#' @param model_data List containing the processed model data structure created by
#'   \code{\link{create_model_data}}. Should include response data, covariates,
#'   spline basis matrices, and prior specifications.
#' @param M Integer specifying the number of mixture components (groups/clusters).
#' @param chains Integer specifying the number of MCMC chains to run. Default is 2.
#'   Multiple chains are recommended for convergence assessment.
#' @param iters Integer specifying the total number of MCMC iterations per chain.
#'   Default is 1000.
#' @param burn_in Integer specifying the number of initial iterations to discard
#'   as burn-in. Default is 500.
#' @param thin Integer specifying the thinning interval. Only every \code{thin}-th
#'   iteration is retained. Default is 5.
#' @param lambda Numeric value for the ridge regularization parameter. If NULL
#'   (default), lambda will be automatically calibrated using cross-validation
#'   across the \code{lambda_grid} values in \code{init_control}.
#' @param init_control List of settings controlling initialization and lambda
#'   calibration. Elements:
#'   \describe{
#'     \item{\code{lambda_init}}{Numeric initial lambda value used during
#'       initialization when \code{lambda = NULL}. Default is 1.}
#'     \item{\code{n_init}}{Integer number of random initializations to try.
#'       The one with highest log-posterior is selected. Default is 5.}
#'     \item{\code{lambda_grid}}{Numeric vector of lambda values to evaluate
#'       during calibration. Default is \code{seq(0.01, 5, length.out = 30)}.}
#'     \item{\code{init_iters}}{Integer number of MCMC iterations per
#'       initialization attempt. Default is 10.}
#'     \item{\code{init_burn_in}}{Integer burn-in iterations for each
#'       initialization run. Default is 5.}
#'     \item{\code{init_thin}}{Integer thinning interval for initialization
#'       runs. Default is 2.}
#'     \item{\code{init_final_run}}{Integer number of iterations for the final
#'       initialization run used to obtain starting values. Default is 100.}
#'     \item{\code{verbose}}{Logical whether to print progress during
#'       initialization. Default is FALSE.}
#'   }
#' @param verbose Logical indicating whether to print progress messages during
#'   execution. Default is TRUE.
#' @param seed Integer seed for random number generation to ensure reproducibility.
#'
#' @details
#' The pipeline executes the following steps:
#' \enumerate{
#'   \item \strong{Initialization}: Runs \code{n_init} short MCMC chains with
#'         different random starting values and selects the one with highest
#'         log-posterior probability.
#'   \item \strong{Lambda Calibration}: If \code{lambda = NULL}, evaluates model
#'         fit across \code{lambda_grid} values and selects the optimal regularization
#'         parameter based on predictive performance.
#'   \item \strong{Main MCMC Sampling}: Runs \code{chains} independent MCMC chains
#'         starting from the best initialization, each with \code{iters} iterations.
#'   \item \strong{Convergence Assessment}: Checks convergence of log-posterior
#'         and regression coefficients (beta) across chains.
#'   \item \strong{Post-processing}: Combines chains, computes predicted probabilities
#'         across the spline basis, and evaluates model fit metrics.
#' }
#'
#' @return A list with the following components:
#' \describe{
#'   \item{\code{sample_list}}{List containing posterior samples combined across
#'         all chains after burn-in and thinning. Includes:
#'         \itemize{
#'           \item \code{beta}: Array of regression coefficient samples
#'           \item \code{w}: Matrix of cluster assignment samples
#'           \item \code{pw}: Matrix of mixing proportion samples
#'           \item \code{logpost}: Matrix of log-posterior values
#'           \item \code{spline_probs}: Array of predicted probabilities across spline basis
#'         }}
#'   \item{\code{lambda_opt}}{List with optimal lambda value. Contains:
#'         \itemize{
#'           \item \code{best_lambda}: The selected regularization parameter
#'         }}
#'   \item{\code{conv_check}}{List of convergence diagnostics:
#'         \itemize{
#'           \item \code{logpost}: Convergence statistics for log-posterior
#'           \item \code{beta}: Convergence statistics for regression coefficients
#'         }}
#'   \item{\code{metrics}}{Data frame of model evaluation metrics including
#'         information criteria, predictive performance, and clustering quality measures.}
#'   \item{\code{w}}{Vector of posterior mode cluster assignments for each observation.}
#' }
#'
#' @seealso
#' \code{\link{fit_mixscat}} for the high-level user interface,
#' \code{\link{create_model_data}} for data preparation,
#' \code{\link{find_init_w}} for initialization,
#' \code{\link{calibrate_lambda}} for lambda selection
#'
#' @examples
#' \dontrun{
#' # Prepare model data
#' model_data <- create_model_data(
#'   z = response_vector,
#'   id = subject_ids,
#'   time = time_points,
#'   n_basis = 10
#' )
#'
#' # Run pipeline for a 3-component mixture
#' result <- pipeline(
#'   model_data = model_data,
#'   M = 3,
#'   chains = 2,
#'   iters = 2000,
#'   burn_in = 1000,
#'   thin = 5,
#'   seed = 123
#' )
#'
#' # Check convergence
#' print(result$conv_check)
#'
#' # View model metrics
#' print(result$metrics)
#'
#' # Extract cluster assignments
#' clusters <- result$w
#' }
#'
#' @export
pipeline = function(model_data,
                    M,
                    chains,
                    init_list,
                    iters,
                    burn_in,
                    thin,
                    lambda,
                    n_cores,
                    seed,
                    verbose = TRUE) {

  if(!is.null(seed)) set.seed(seed)

  #### MCMC ####
  if(verbose) cat("Running chains \n")

  # set number of cores for parallel chains
  if(n_cores > 1) {
    future::plan(future::multisession,  workers = n_cores)
  }else{
    future::plan(future::sequential, split = TRUE)
  }

  if(verbose) cat("Number of workers:", future::nbrOfWorkers(), "\n")

  # run chains
  runs = future.apply::future_lapply(1:chains, function(chain){

        if(verbose) cat("Running for chain", chain, "\n")

        mixscat::single_run(
          M = M,
          w = NULL,
          model_data = model_data,
          lambda = lambda,
          init_list = init_list,
          iters = iters,
          burn_in = burn_in,
          thin = thin,
          verbose = verbose,
          seed = NULL
      )

      },

      future.seed = TRUE,
      future.packages = c("Rcpp", "dplyr", "tidyr", "purrr", "posterior"),
      future.stdout = TRUE

    )

  init_time = Sys.time()
  #### check label swtiching ####
  runs = relabel(runs)

  #### convergence check ####
  logpost_conv_check = check_conv_logpost(runs = runs)
  beta_conv_check = check_conv_param(runs = runs, param_name = "beta")

  conv_check = list(
    logpost = logpost_conv_check,
    beta = beta_conv_check
  )

  #### combine chains ####
  sample_list = combine_chains(runs = runs, model_data)

  for(i in 1:nrow(sample_list$logpost)) {
    sample_list$spline_probs[i,,] = predict_prob_cpp(
      M = M, w = 1:M, B = model_data$spline$B, beta = sample_list$beta[i,,]
    )
  }

  #### clustering quality metrics ####
  metrics = compute_metrics(
    sample_list = sample_list,
    lambda = lambda,
    model_data = model_data
  )

  end_time = Sys.time()
  runtime = end_time - init_time

  #### return ####
  out = list(
    sample_list = sample_list,
    conv_check = conv_check,
    metrics = metrics,
    w = sample_list$w |> comp_class(),
    runs = runs,
    runtime = runtime
  )

  return(out)

}
