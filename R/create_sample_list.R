create_sample_list = function(M,
                              iters = 1000,
                              burn_in = iters/2,
                              thin = 1,
                              model_data,
                              init_list = NULL,
                              seed = NULL) {

  n = model_data$dims$n
  G = model_data$dims$G
  n_id = model_data$dims$n_id
  n_basis = model_data$spline$n_basis
  B = model_data$spline$B
  n_time = model_data$dims$n_time

  id_unique = model_data$data$id_unique

  if(!is.null(seed)) set.seed(seed)

  iters_vec = seq(from = burn_in + 1, to = iters, by = thin)
  eff_iters = length(iters_vec)

  sample_list = list()

  #### global clustering parameters ####
  sample_list$beta = gen_sample_array(
    iters = eff_iters,
    dimension = c(M*n_basis, G),
    sampler = function(x) rnorm(x, sd = 0.01),
    init = init_list$beta
  )

  sample_list$beta[1,,G] = 0

  sample_list$w = gen_sample_array(
    iters = eff_iters,
    dimension = c(n_id),
    sampler = function(x) {
      w = sample(1:M, size = n_id, replace = TRUE)
      names(w) = id_unique
      return(w)
    },
    init = init_list$w
  )

  colnames(sample_list$w) = id_unique

  sample_list$w_post_prob = gen_sample_array(
    iters = eff_iters,
    dimension = c(n_id, M),
    sampler = function(x) 1/M,
    init = init_list$w_post_prob
  )

  sample_list$pw = gen_sample_array(
    iters = eff_iters,
    dimension = c(M),
    sampler = function(x) rep(1/M, M),
    init = init_list$pw
  )

  sample_list$spline_probs = gen_sample_array(
    iters = eff_iters,
    dimension = c(M * n_time, G),
    sampler = function(x) predict_prob_cpp(M = M, w = 1:M, B, sample_list$beta[1,,]),
    init = init_list$spline_probs
  )

  #### log posterior ####
  sample_list$logpost = gen_sample_array(
    iters = eff_iters,
    dimension = c(4),
    sampler = function(x) 0,
    init = NULL
  )

  colnames(sample_list$logpost) = c("logpost", "loglike", "logpenal", "beta_logprior")

  sample_list$iters_vec = iters_vec

  return(sample_list)

}

