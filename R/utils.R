gv = c(
  ".", "G", "M", "nvars", "w", "time", "log_penal",
  "li", "ui", "id", "logp_z", "lambda", "model", "mode_diff",
  "min_id", "max_id", "mid_point", "y", "variable"
)

utils::globalVariables(gv)

stick_breaking = function(x) {
  x * c(1, cumprod(1 - x[-length(x)]))
}

complete_dim = function(array, dimension) {
  if(length(dim(array)) == 0) {
    if(length(array) < dimension) {
      array = c(array, rep(NA, dimension - length(array)))
    }
  }else if(length(dim(array)) == 2) {
    if(dim(array)[2] < dimension) {
      array = cbind(array, matrix(NA, nrow(array), dimension - dim(array)[2]))
    }
  }

  return(array)
}


logit = function(x) {
  log(x/(1-x))
}

mse_score = function(y_true, y_pred) {
  mean((y_true - y_pred)^2)
}

rmse_score = function(y_true, y_pred) {
  sqrt(mse_score(y_true, y_pred))
}

corr_score = function(y_true, y_pred) {
  diag(cor(y_true, y_pred))
}

compute_theta = function(mu, Z, E, K) {
  n = nrow(Z)
  G = ncol(Z)
  K = ncol(cbind(E))

  if(K == 1 & G > 1) mu = cbind(mu)
  if(K > 1 & G == 1) mu = rbind(mu)

  mu_matrix = Z %*% mu
  theta = mu_matrix + E
  return(theta)
}

create_dummy = function(z, G) {

  if(G != 1) {

    Z = fast_dummy_dense(z, G)

  }else{
    Z = matrix(1, nrow = length(z), ncol = 1)
  }

  return(Z)
}

filter_chain = function(chain, thin, burn_in) {

  iters = dim(chain)[1]
  index = seq(burn_in+1, iters, thin)
  dims = length(dim(chain))

  if(dims == 3) {
    chain = chain[index,,]
  }else if(dims == 2){
    chain = chain[index,]
  }else if(dims == 1){
    chain = chain[index]
  }

  return(chain)

}


#' Generate Initial Sample Array
#'
#' Creates an empty array to store MCMC samples, initialized with optional initial values.
#'
#' @param iters Integer. The total number of iterations for the MCMC chain.
#' @param dimension Integer vector. The dimensions of the array for each iteration's sample
#'   (excluding the iteration dimension).
#' @param sampler Function. An optional function to generate initial values if `init` is `NULL`.
#' @param init Numeric array. Optional initial values for the first iteration.
#'
#' @return A numeric array of dimensions `c(iters, dimension)` with the first slice
#'   (corresponding to the first iteration) potentially filled with `init`.
#' @keywords internal
gen_sample_array = function(iters, dimension, sampler = NULL, init = NULL) {

  dim_len = length(dimension)
  sample_array = array(dim = c(iters, dimension))


  if(is.null(init)) {
    init = sampler(prod(dimension))
  }

  if(dim_len == 0) {
    sample_array[1] = init
  }else if(dim_len == 1) {
    sample_array[1,] = as.numeric(init)
  }else if(dim_len == 2) {
    sample_array[1,,] = init
  }

  return(sample_array)

}

compute_mode = function(x) {
  x |> table() %>% which.max() %>% names() %>% as.integer()
}

compute_post_stat = function(sample, stat_function = mean, ...) {

  dims = length(dim(sample))

  if(!is.null(sample)) {

    if(dims == 3) {
      est = lapply(1:dim(sample)[3], function(g){
        cbind(sample[, , g]) |> apply(MARGIN = 2, FUN = stat_function, ...)
      }) %>% do.call(cbind, .)

    }else if(dims == 2) {

      est = apply(sample, FUN = stat_function, MARGIN = 2,  ...)

    }else if(dims == 1) {

      est = stat_function(sample, ...)

    }

  }else{
    est = NULL
  }

  return(est)

}

compute_posterior_mean = function(sample_list) {
  purrr::map(sample_list, compute_post_stat)
}


get_rotation_list = function(post_sample, reference_matrix) {
  rotation_list = list()
  iters = dim(post_sample)[1]
  K = dim(post_sample)[3]

  if(K == 1) {
    for(i in 1:iters) {
      rotation_list[[i]] = vegan::procrustes(Y = cbind(post_sample[i,,]),
                                             X = cbind(reference_matrix),
                                             scale = F)
    }
  }else{
    for(i in 1:iters) {
      rotation_list[[i]] = vegan::procrustes(Y = post_sample[i,,],
                                             X = reference_matrix,
                                             scale = F)
    }
  }

  return(rotation_list)
}

apply_rotation = function(post_sample, rotation_list) {
  new_post_sample = post_sample
  iters = dim(post_sample)[1]
  K = dim(post_sample)[3]

  if(K == 1) {
    for(i in 1:iters) {
      r = rotation_list[[i]]$rotation
      new_post_sample[i,,] = cbind(post_sample[i,,]) %*% r
    }
  }else{
    for(i in 1:iters) {
      r = rotation_list[[i]]$rotation
      new_post_sample[i,,] = post_sample[i,,] %*% r
    }
  }

  return(new_post_sample)
}

get_adapted_rotation_list = function(post_sample, reference_matrix, H_vec) {

  rotation_list = list()
  iters = dim(post_sample)[1]
  K = dim(post_sample)[3]

  if(K == 1) {
    for(i in 1:iters) {
      rotation_list[[i]] = vegan::procrustes(
        Y = cbind(post_sample[i,, 1:H_vec[i]]),
        X = cbind(reference_matrix[, 1:(H_vec[i]-1)], 0),
        scale = F
      )
    }
  }else{
    for(i in 1:iters) {
      rotation_list[[i]] = vegan::procrustes(
        Y = cbind(post_sample[i,, 1:H_vec[i]]),
        X = cbind(reference_matrix[, 1:(H_vec[i]-1)], 0),
        scale = F
      )
    }
  }

  return(rotation_list)
}

apply_adapted_rotation = function(post_sample, rotation_list, H_vec) {
  new_post_sample = post_sample
  iters = dim(post_sample)[1]
  K = dim(post_sample)[3]

  if(K == 1) {
    for(i in 1:iters) {
      r = rotation_list[[i]]$rotation
      new_post_sample[i,, 1:H_vec[i]] = cbind(post_sample[i,, 1:H_vec[i]]) %*% r
    }
  }else{
    for(i in 1:iters) {
      r = rotation_list[[i]]$rotation
      new_post_sample[i,, 1:H_vec[i]] = post_sample[i,, 1:H_vec[i]] %*% r
    }
  }

  return(new_post_sample)
}


gv = c(
  ".", "G", "M", "nvars", "w", "time", "log_penal",
  "li", "ui", "id", "logp_z", "lambda", "model", "mode_diff",
  "min_id", "max_id", "mid_point"
)

utils::globalVariables(gv)

#' Indices of Non-Penalized Basis Functions
#'
#' Generates indices for non-penalized basis functions across groups or classes.
#'
#' @param n_basis Integer. Number of basis functions per group.
#' @param M Integer. Number of groups or classes.
#' @param order Integer. Order of difference for penalization. Default is 1.
#'
#' @return Integer vector of indices.
#' @keywords internal
#' @examples
#' \dontrun{
#' gen_notpen_index(5, 3, order = 2)
#' }
gen_notpen_index = function(n_basis, M, order = 1) {
  notpen_index = lapply(1:M, function(k){
    (k -1) * n_basis + 1:(order)
  }) %>% do.call(c, .)
  return(notpen_index)
}

#' Indices for Grouped Basis Functions
#'
#' Creates a list of index vectors for basis functions grouped by class.
#'
#' @param n_basis Integer. Number of basis functions per group.
#' @param M Integer. Number of groups.
#'
#' @return List of integer vectors, one per group.
#' @keywords internal
#' @examples
#' \dontrun{
#' gen_basis_index(5, 3)
#' }
gen_basis_index = function(n_basis, M) {
  lapply(1:M, function(k){
    (((k-1)*(n_basis) + 1):(k * (n_basis)))
  })
}

#' Generate Inverse Covariance Matrices for Penalized Coefficients
#'
#' Creates a list of diagonal inverse covariance matrices for each class (except the last).
#'
#' @param lambda Numeric matrix or vector of penalty parameters.
#' @param model_data List. Model structure (see create_model_data).
#' @param fixed_sd Numeric. Standard deviation for non-penalized coefficients. Default is 100.
#'
#' @return List of diagonal matrices.
#' @keywords internal
#' @examples
#' \dontrun{
#' gen_inv_cov(c(1, 2, 3), model_data = list(), fixed_sd = 100)
#' }
gen_inv_cov = function(lambda, model_data, fixed_sd = 100) {

  G = model_data$G
  M = model_data$M
  n_basis = model_data$n_basis
  notpen_index = model_data$notpen_index
  n_id = model_data$n_id
  nD = model_data$nD
  lambda = cbind(lambda)

  var_vec = rep((lambda), each = n_basis * M)
  D = rep(1/diag(nD), times = M)
  var_vec = var_vec * D
  var_vec[notpen_index] = (fixed_sd^2)
  prec_matrix = diag(1/var_vec)

  inv_cov_list = lapply(1:(G-1), function(g){prec_matrix})
  return(inv_cov_list)

}

#' Regularized Design Matrix for Cluster-Specific Effects
#'
#' Constructs a design matrix for cluster-specific basis coefficients.
#'
#' @param w_vec Integer vector. Cluster assignment for each row.
#' @param M Integer. Number of clusters.
#' @param B_expand Matrix. Basis expansion for all clusters.
#' @param intercept Logical. If TRUE, includes an intercept column. Default is FALSE.
#'
#' @return Numeric matrix, same shape as B_expand.
#' @keywords internal
#' @examples
#' \dontrun{
#' gen_reg_matrix(c(1, 2, 3), 3, matrix(1, 3, 3), intercept = TRUE)
#' }
gen_reg_matrix = function(w_vec, M, B_expand, intercept = FALSE) {

  n_basis = ncol(B_expand)/(M)

  W = gen_dummy(w_vec, n_cat = M, intercept = intercept)
  W_expand = W[, rep(1:M, each = n_basis)]

  X = B_expand * W_expand

  return(X)
}

#' Posterior Summary Statistic by Group
#'
#' Computes a summary statistic (e.g., mean) across posterior samples for each group.
#'
#' @param sample 3D numeric array: iterations x items x groups.
#' @param stat_function Function. Statistical function to apply (e.g., mean). Default is mean.
#' @param ... Additional arguments passed to the statistical function.
#'
#' @return Numeric matrix: items x groups.
#' @keywords internal
#' @examples
#' \dontrun{
#' comp_post_stat(array(1:27, dim = c(3, 3, 3)), stat_function = mean)
#' }
comp_post_stat = function(sample, stat_function = mean, ...) {
  if(!is.null(sample)) {
    est = lapply(1:dim(sample)[3], function(g){
      sample[, , g] |> apply(MARGIN = 2, FUN = stat_function, ...)
    }) %>%
      do.call(cbind, .)
  }else{
    est = NULL
  }

  return(est)
}

#' Most Frequent Class Assignment (Posterior Mode)
#'
#' Computes the most frequent class for each item across samples.
#'
#' @param sample Matrix or array: iterations x items.
#'
#' @return Integer vector of modal class assignments.
#' @keywords internal
#' @examples
#' \dontrun{
#' comp_class(matrix(sample(1:3, 30, replace = TRUE), nrow = 10))
#' }
comp_class = function(sample) {

  class = apply(
    sample,
    MARGIN = 2,
    FUN = function(x) x |> table() |> which.max() |> names() |> as.integer()
  )

  return(class)

}

#' Convert Posterior Summary Matrix to Long Format
#'
#' Internal: Convert summary matrix to tidy long-format tibble for plotting/analysis.
#'
#' @param summ_matrix Numeric matrix. Posterior summary (e.g., means).
#' @param col_name Character. Name for summary column.
#' @param time Numeric/integer vector. Time points.
#' @param w Integer/factor vector. Group labels.
#'
#' @return Tibble in long format: time, w, cat, summary.
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate
#' @importFrom tidyr gather
#' @importFrom stringr str_replace
#' @keywords internal
trans_summ_prob = function(summ_matrix, col_name, time, w) {
  summ_matrix |>
    as.data.frame() |>
    tibble::as_tibble() |>
    dplyr::mutate(time = time, w = w) |>
    tidyr::gather(cat, {{col_name}}, -time, -w) |>
    dplyr::mutate(cat = str_replace(cat, "V", ""))
}

#' Compute Model Selection Metrics from MCMC Runs
#'
#' Internal: Calculate AICM, DIC and BICM metrics based on penalized log-posterior samples.
#'
#' @param runs List of MCMC run objects. Each element must contain a `logpost` matrix
#'   with column `"penal_logpost"`.
#' @param n Integer. Number of observations (overwritten if `model_data_min$z` is present).
#' @param model_data_min List. Minimal model data; must include either `n` (if `z` is `NULL`)
#'   or `n_id` (if `z` is present).
#'
#' @return List with numeric vectors:
#' \describe{
#'   \item{l_mean}{Mean penalized log-posterior for each run.}
#'   \item{l_var}{Variance of penalized log-posterior for each run.}
#'   \item{AICM}{Akaike Information Criterion for MCMC.}
#'   \item{BICM}{Bayesian Information Criterion for MCMC.}
#' }
#' @importFrom purrr map_dbl
#' @importFrom stats var
#' @keywords internal
comp_metrics = function(runs, n, model_data_min) {

  if(is.null(model_data_min$z)) {
    n = model_data_min$n
  }else{
    n = model_data_min$n_id
  }

  l_mean = runs %>% purrr::map_dbl(~{mean(.x$logpost[, "penal_logpost"])})
  l_var = runs %>% purrr::map_dbl(~{stats::var(.x$logpost[, "penal_logpost"])})

  AICM = -2*(l_mean - l_var)
  BICM = -2*(l_mean - (log(n)-1)*l_var)

  l_post = runs %>% purrr::map_dbl(~{.x$logpost_est[1]})
  pp = l_post - l_mean
  DIC2 = -2*(l_mean - pp)
  DIC3 = -2*(l_mean - (log(n)-1)*pp)

  out = list(
    l_mean = l_mean,
    l_var = l_var,
    l_post = l_post,
    AICM = AICM,
    BICM = BICM,
    DIC2 = DIC2,
    DIC3 = DIC3,
    pp = pp
  )

  out

}

#' Format Model Selection Metrics
#'
#' Internal: Format the list of model selection metrics into a tibble.
#'
#' @param metrics_list List. Output from `comp_metrics`.
#'
#' @return A tibble with columns "G", "M", and the calculated metrics (l_mean, l_var, AICM, BICM).
#' @importFrom tibble as_tibble
#' @importFrom dplyr mutate
#' @importFrom tidyr separate
format_metrics = function(metrics_list) {

  metrics_df = metrics_list %>%
    as.data.frame() %>%
    mutate(model = rownames(.)) %>%
    tidyr::separate(model, into = c("G", "M"), sep = ",") |>
    mutate(G = stringr::str_extract(G, "\\d{1,10}"),
           M = stringr::str_extract(M, "\\d{1,10}"))
  rownames(metrics_df) = NULL

  return(metrics_df)
}

#' Filter MCMC Samples
#'
#' Internal: Filters an array of MCMC samples by applying a burn-in period and thinning rate.
#'
#' @param array A numeric array containing MCMC samples (first dimension is iterations).
#' @param burn_in Integer. The number of initial iterations to discard.
#' @param thin Integer. The thinning rate; keep every `thin` iteration after burn-in.
#'
#' @return A filtered numeric array, or `NULL` if the input array was `NULL`.
#'   The dimensions of the output array will depend on the original array's dimensions.
#'   For a 1D array, it returns a vector. For 2D and 3D arrays, it returns a matrix or array,
#'   respectively.
#' @keywords internal
filter_array = function(array, burn_in, thin) {

  if(!is.null(array)) {

    array_dim = dim(array)
    dim_len = length(array_dim)
    iters = array_dim[1]
    post_iters = seq(burn_in, iters, thin)

    if(dim_len == 1) {
      post_array = array[post_iters]
    }else if(dim_len == 2) {
      post_array = array[post_iters, ]
    }else if(dim_len == 3) {
      post_array = array[post_iters,,]
    }

    colnames(post_array)

  }else{
    post_array = NULL
  }

  return(post_array)
}

#' Generate Initial Sample Array
#'
#' Internal: Creates an empty array to store MCMC samples, initialized with optional initial values.
#'
#' @param iters Integer. The total number of iterations for the MCMC chain.
#' @param dimension Integer vector. The dimensions of the array for each iteration's sample
#'   (excluding the iteration dimension).
#' @param sampler Function. An optional function to generate initial values if `init` is `NULL`.
#' @param init Numeric array. Optional initial values for the first iteration.
#'
#' @return A numeric array of dimensions `c(iters, dimension)` with the first slice
#'   (corresponding to the first iteration) potentially filled with `init`.
#' @keywords internal
gen_sample_array = function(iters, dimension, sampler = NULL, init = NULL) {

  dim_len = length(dimension)
  sample_array = array(dim = c(iters, dimension))


  if(is.null(init)) {
    init = sampler(prod(dimension))
  }

  if(dim_len == 0) {
    sample_array[1] = init
  }else if(dim_len == 1) {
    sample_array[1,] = init
  }else if(dim_len == 2) {
    sample_array[1,,] = init
  }

  return(sample_array)

}

#' Create Minimal Model Data List
#'
#' Internal: Creates a reduced list containing essential elements from the full model data.
#' This is used to pass a smaller object between functions when the full `model_data` is not required.
#'
#' @param model_data List. The full model data object (see `create_model_data`).
#' @param M Integer. The number of clusters for the current model.
#' @param G Integer. The number of categories (outcome levels) for the current model.
#'
#' @return A list containing a subset of the elements from `model_data` (`M`, `G`, `w`, `z`, `x`, `B`, `P`, `n_vars`, `n_id`, `n_time`, `n`).
#' @keywords internal
create_model_data_min = function(model_data, M, G) {
  model_data_min = list(
    B = model_data$B_unique,
    P = model_data$nD,
    n_vars = model_data$n_vars,
    n_id = model_data$n_id,
    n_time = model_data$n_time,
    n = model_data$n
  )

  return(model_data_min)
}

#' Compute Hamming Distance
#'
#' Internal: Computes the Hamming distance between sequences represented in a vector `z`.
#'
#' @param z Integer vector. The sequence data, flattened (individuals x time points).
#' @param model_data List. The model data object, containing `n_time` and `n_id`.
#'
#' @return A `dist` object containing the pairwise Hamming distances between sequences.
#' @importFrom TraMineR seqdef seqdist
#' @keywords internal
compute_hamming = function(z, model_data) {
  n_time = model_data$n_time
  n_id = model_data$n_id

  z_seq_matrix = matrix(z, ncol = n_time, nrow = n_id, byrow = T)

  z_ham_dist = diss = suppressMessages(
    suppressWarnings(
      z_seq_matrix |>
        TraMineR::seqdef() |>
        TraMineR::seqdist(, method = "HAM") |>
        as.dist()
    )
  )

  return(z_ham_dist)
}

format_best_lambda = function(best_lambda_list) {
  best_lambda_list %>%
    map_dbl(~{.x}) %>%
    data.frame(model = names(.), lambda = .) %>%
    tidyr::separate(col = model, into = c("G", "M"), sep = "\\,") %>%
    mutate(G = stringr::str_remove(G, "G=") %>% as.integer(),
           M = stringr::str_remove(M, " M=") %>% as.integer()) %>%
    as_tibble()
}

comp_prior_mean = function(w, beta_theta, B_unique, M) {

  W = gen_dummy(w, M, intercept = FALSE)
  prior_mean = kronecker(W, B_unique) %*% beta_theta
  #prior_mean = scale(prior_mean, center = TRUE, scale = FALSE)
  return(prior_mean)
}

gen_normal_mat = function(n, k) {
  matrix(rnorm(n * k), nrow = n, ncol = k)
}

compute_probs = function(w, M, B, beta) {

  W = create_dummy(w, M)
  X = cbind(1, kronecker(W, B))
  prob = mclust::softmax(X %*% beta)

  return(prob)
}




get_w_ward = function(M, z_dist) {
  w_ward_init = hclust(as.dist(z_dist), method = "ward.D") |> cutree(k = M)
  w_ward_init
}

compute_model_precision = function(model_data,
                                   lambda,
                                   intercept_penalty,
                                   M) {
  n_basis = model_data$spline$n_basis
  S = model_data$spline$S
  S[1, 1] = 1

  if(length(lambda) == 1) {

    lambda_expand = rep(lambda, M * n_basis)

  }else{

    lambda_expand = rep(lambda, each = n_basis)

  }

  s = diag(S)
  diag_precision = rep(s, times = M) * lambda_expand

  intercept_idx = seq(1, M * n_basis, by = n_basis)
  diag_precision[intercept_idx] = intercept_penalty

  return(diag_precision)
}

create_beta_precision_matrix = function(model_data,
                                        lambda,
                                        intercept_penalty,
                                        M) {

  diag_precision = compute_model_precision(
    model_data = model_data,
    lambda = lambda,
    intercept_penalty = intercept_penalty,
    M = M
  )

  beta_precision_matrix = diag(diag_precision)

  return(beta_precision_matrix)

}

create_beta_sd_matrix = function(model_data,
                                 beta_precision_matrix = NULL,
                                 lambda,
                                 intercept_penalty,
                                 M) {

  G = model_data$dims$G
  n_basis = model_data$spline$n_basis

  if(is.null(beta_precision_matrix)) {

    diag_precision = compute_model_precision(
      model_data = model_data,
      lambda = lambda,
      intercept_penalty = intercept_penalty,
      M = M
    )

  }else{
    diag_precision = diag(beta_precision_matrix)
  }

  beta_sd = 1/sqrt(diag_precision)

  beta_sd = matrix(
    beta_sd,
    nrow = length(beta_sd),
    ncol = G - 1,
    byrow = F
  )

  return( beta_sd)

}

compute_entropy = function(pw, normalize = TRUE) {

  h = -rowSums(pw * log(pw), na.rm = T)
  M = ncol(pw)
  if(normalize == TRUE) h = h/log(M)
  return(h)

}



