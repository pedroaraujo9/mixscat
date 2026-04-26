#' Find the number of mixture components
#'
#' @description
#' Runs multiple short MCMC chains for an overfitted mixture model with
#' maximum number of components \code{M_max} and uses clustering stability
#' and entropy-based summaries to choose the most supported number of active
#' clusters.
#'
#' @param M_max Integer giving the maximum number of mixture components used
#'   in the overfitted model.
#' @param n_init Integer number of independent initial runs (short MCMC
#'   chains) to perform. Default is 50.
#' @param iters Integer number of MCMC iterations per initialization run.
#'   Default is 10.
#' @param lambda Numeric regularization parameter for the spline coefficients.
#' @param dirichlet_param Numeric hyperparameter for the Dirichlet prior on
#'   mixture weights.
#' @param a_lambda Numeric shape parameter of the Gamma prior for \code{lambda}.
#'   Default is 1.
#' @param b_lambda Numeric rate parameter of the Gamma prior for \code{lambda}.
#'   Default is 1.
#' @param intercept_penalty Numeric penalty weight for the intercept term in
#'   the spline basis. Default is 1.
#' @param model_data List containing the processed model data as returned by
#'   \code{create_model_data()}.
#' @param seed Optional integer seed for reproducibility of the
#'   initialization runs. If \code{NULL}, no seed is set.
#' @param verbose Logical indicating whether to print progress information for
#'   the initialization runs. Default is TRUE.
#'
#' @return A list with components
#'   \item{runs}{List of \code{n_init} fitted short-run objects returned by
#'     \code{single_run()}.}
#'   \item{clust_size}{Matrix with the number of active clusters at each
#'     iteration for each initialization run.}
#'   \item{clust_size_p}{Data frame with the empirical distribution of the
#'     number of active clusters at the last iteration across runs.}
#'   \item{best_size}{Integer giving the selected number of clusters.}
#'   \item{w_sample}{Matrix of cluster allocations at the last iteration of
#'     each run (rows correspond to runs, columns to individuals).}
#'   \item{w_sample_best}{Subset of \code{w_sample} restricted to runs with
#'     \code{best_size} active clusters, after relabeling.}
#'   \item{w_pear}{Vector of cluster labels given by the posterior expected
#'     partition (PEAR) clustering.}
#'   \item{w_post_prob}{Matrix of posterior cluster membership probabilities
#'     for each individual and cluster.}
#'   \item{entropy}{Data frame with entropy summaries across iterations and
#'     runs.}
#'   \item{prob_last}{Data frame with the probability that each individual
#'     keeps the same allocation as in the last iteration, across runs.}
#'
#' @export
find_number_clust = function(M_max,
                             z,
                             id,
                             time,
                             n_basis,
                             n_init = 50,
                             init_iters = 10,
                             lambda,
                             dirichlet_param = 0.01,
                             a_lambda = 1,
                             b_lambda = 1,
                             intercept_penalty = 1,
                             seed,
                             verbose = TRUE) {

  if(!is.null(seed)) set.seed(seed)

  model_data = create_model_data(
    z = z,
    id = id,
    time = time,
    n_basis = n_basis,
    intercept_penalty = intercept_penalty,
    dirichlet_param = dirichlet_param
  )

  n_id = model_data$dims$n_id
  id_unique = model_data$data$id_unique

  model_path = system.file(
    "stan", "stan_full_logpost.rds", package = "mixscat"
  )

  stan_model = readRDS(model_path)


  runs = lapply(1:n_init, function(i){

    if(verbose) cat("Init run:", i, "\r")

    single_run(
      M = M_max,
      w = NULL,
      model_data = model_data,
      intercept_penalty = intercept_penalty,
      dirichlet_param = dirichlet_param,
      lambda = lambda,
      a_lambda = a_lambda,
      b_lambda = b_lambda,
      init_list = NULL,
      iters = init_iters,
      burn_in = 0,
      thin = 1,
      verbose = FALSE,
      seed = NULL
    )
  })

  entropy = lapply(runs, function(run){

    lapply(1:init_iters, function(i){

      run$sample_list$w_post_prob[i, ,] %>% compute_entropy() %>% sum()

    }) %>% do.call(c, .)

  }) %>% do.call(rbind, .) %>%
    as.data.frame() %>%
    mutate(runs = 1:nrow(.)) %>%
    gather(iter, entropy, -runs) %>%
    mutate(iter = as.integer(gsub("V", "", iter)))

  clust_size = lapply(runs, function(run){
    lapply(1:init_iters, function(j){
      run$sample_list$w[j, ] %>% unique() %>% length()
    }) %>% do.call(c, .)
  }) %>%
    do.call(rbind, .)

  w_last = lapply(runs, function(run){
    run$sample_list$w[init_iters, ]
  }) %>%
    do.call(rbind, .)

  post_prob_list = lapply(runs, function(run){
    run$sample_list$w_post_prob[init_iters, ,]
  })

  prob_last = lapply(runs, function(run){
    eq = run$sample_list$w == matrix(run$sample_list$w[init_iters, ], nrow = init_iters, ncol = n_id, byrow = T)
    colMeans(eq)
  }) %>%
    do.call(rbind, .)

  prob_last = prob_last %>%
    as.data.frame() %>%
    mutate(runs = 1:nrow(.)) %>%
    gather(id, prob, -runs)

  post_prob = array(dim = c(n_init, nrow(post_prob_list[[1]]), ncol(post_prob_list[[1]])))

  for(i in 1:n_init) {

    post_prob[i,,] = post_prob_list[[i]]

  }


  #### number of active cluster distribution ####
  p = clust_size[, init_iters] %>% table() %>% prop.table()

  clust_size_p = data.frame(
    n_clust = as.numeric(names(p)),
    prop = as.numeric(p)
  )

  sizes = clust_size_p$n_clust

  w_size = apply(w_last, MARGIN = 1, FUN = function(x) {
    x %>% unique() %>% length()
  })

  post_modes = lapply(sizes, function(size){

    #### condionate on the chosen model #####
    w_m = rbind(w_last[w_size == size, ])
    n_best = nrow(w_m)

    for(i in 1:n_best){
      w_m[i, ]  = w_m[i, ] %>% factor(labels = 1:size) %>% as.integer()
    }

    psm = mcclust::comp.psm(cls = w_m)
    colnames(psm) = rownames(psm) = id_unique

    cl = mcclust::maxpear(psm, cls.draw = w_m, max.k = size, method = "draws")
    w_pear = cl$cl
    names(w_pear) = id_unique

    ls = label.switching::label.switching(
      method = "ECR",
      z = w_m,
      zpivot = w_pear,
      K = size
    )

    for(i in 1:n_best) {

      perm = ls$permutations$`ECR`[i, ]
      w_m[i, ] = order(perm)[w_m[i, ]]

    }

    w_post_prob = lapply(1:n_id, function(i){
      w_m[, i] %>%
        factor(levels = 1:size) %>%
        table() %>%
        prop.table() %>%
        as.numeric()
    }) %>%
      do.call(rbind, .)

    rownames(w_post_prob) = id_unique
    w = w_post_prob %>% apply(MARGIN = 1, FUN = which.max)

    post_map = find_post_map(
      M = M_max,
      w = w,
      stan_model = stan_model,
      init_list = NULL,
      lambda = lambda,
      intercept_penalty = intercept_penalty,
      dirichlet_param = dirichlet_param,
      model_data = model_data
    )

    beta_post_map = find_post_map(
      M = size,
      w = w,
      stan_model = stan_model,
      init_list = NULL,
      lambda = lambda,
      intercept_penalty = intercept_penalty,
      dirichlet_param = dirichlet_param,
      model_data = model_data
    )

    out = list(
      logpost_map = post_map$par$log_like,
      beta_map = cbind(beta_post_map$par$beta, 0),
      w = w,
      w_pear = w_pear,
      w_sample = w_m,
      w_post_prob = w_post_prob
    )

    return(out)

  })

  names(post_modes) = sizes

  out = list(
    model_data = model_data,
    post_modes = post_modes,
    clust_size = clust_size,
    clust_size_p = clust_size_p,
    entropy = entropy,
    prob_last = prob_last
  )

  return(out)

}
