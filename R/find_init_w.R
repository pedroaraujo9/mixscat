# find_init_w = function(M,
#                        model_data,
#                        lambda_init = 1,
#                        n_init = 30,
#                        init_iters = 30,
#                        init_burn_in = 15,
#                        init_thin = 2,
#                        verbose = TRUE,
#                        seed) {
#
#   set.seed(seed)
#
#   z_matrix = matrix(
#     model_data$data$z, nrow = model_data$dims$n_id, ncol = model_data$dims$n_time, byrow = T
#   )
#
#   d = z_matrix |> TraMineR::seqdef() |> TraMineR::seqdist(method = "HAM")
#   colnames(d) = rownames(d) = model_data$data$id_unique
#   w0 = hclust(as.dist(d)) |> cutree(k = M) |> sort()
#
#   opt = find_map(
#     M = M,
#     X = kronecker(create_dummy(w0, M), model_data$spline$B),
#     lambda = lambda_init,
#     stan_model = readRDS(system.file(
#       "stan", "stan_for_opt.rds", package = "bmixscat"
#     )),
#     model_data = model_data
#   )
#
#   beta0 = cbind(opt$opt_fit$par$beta, 0)
#
#
#   init_list = list(
#     w = w0,
#     pw = prop.table(table(w0)),
#     # beta = matrix(0, nrow = M * model_data$spline$n_basis, ncol = model_data$dims$G)
#     beta = beta0
#   )
#
#   init_runs = lapply(1:n_init, function(i){
#
#     if(verbose) cat(paste0("init iter: ", i, "\r"))
#
#     run = single_run(
#       M = M,
#       model_data = model_data,
#       lambda = lambda_init,
#       init_list = init_list,
#       iters = init_iters,
#       burn_in = init_burn_in,
#       thin = init_thin,
#       seed = NULL,
#       verbose = FALSE
#     )
#
#     logpost = mean(run$sample_list$logpost[, "logpost"])
#     run$logpost = logpost
#     return(run)
#
#   })
#
#   logpost = purrr::map_dbl(init_runs, ~{.x$logpost})
#   max_logpost = which.max(logpost)
#   w0 = init_runs[[max_logpost]]$sample_list$w |> comp_class()
#   beta0 = init_runs[[max_logpost]]$sample_list$beta |> compute_post_stat()
#   pw0 = init_runs[[max_logpost]]$sample_list$pw |> colMeans()
#
#   out = list(
#     w0 = w0,
#     beta0 = beta0,
#     pw0 = pw0
#   )
#
#   return(out)
#
# }
#
#
# find_init_w2 = function(M,
#                         model_data,
#                         lambda_init = 1,
#                         n_init = 30,
#                         init_iters = 30,
#                         init_burn_in = 15,
#                         init_thin = 2,
#                         verbose = TRUE,
#                         seed) {
#
#   set.seed(seed)
#
#   z_matrix = matrix(
#     model_data$data$z, nrow = model_data$dims$n_id, ncol = model_data$dims$n_time, byrow = T
#   )
#
#   d = z_matrix |> TraMineR::seqdef() |> TraMineR::seqdist(method = "HAM")
#   colnames(d) = rownames(d) = model_data$data$id_unique
#   w0 = hclust(as.dist(d)) |> cutree(k = M) |> sort()
#
#   opt = find_map(
#     M = M,
#     X = kronecker(create_dummy(w0, M), model_data$spline$B),
#     lambda = lambda_init,
#     stan_model = readRDS(system.file(
#       "stan", "stan_for_opt.rds", package = "bmixscat"
#     )),
#     model_data = model_data
#   )
#
#   beta0 = cbind(opt$opt_fit$par$beta, 0)
#
#   out = list(
#     w0 = w0,
#     pw0 = prop.table(table(w0)),
#     beta0 = beta0
#   )
#
#   return(out)
#
# }
#
# find_init_w3 = function(M,
#                         model_data,
#                         lambda_init = 1,
#                         n_init = 30,
#                         init_iters = 30,
#                         init_burn_in = 15,
#                         init_thin = 2,
#                         verbose = TRUE,
#                         seed) {
#
#   set.seed(seed)
#
#   z_matrix = matrix(
#     model_data$data$z, nrow = model_data$dims$n_id, ncol = model_data$dims$n_time, byrow = T
#   )
#
#   d = z_matrix |> TraMineR::seqdef() |> TraMineR::seqdist(method = "HAM")
#   colnames(d) = rownames(d) = model_data$data$id_unique
#   w0 = hclust(as.dist(d)) |> cutree(k = M) |> sort()
#
#   init_runs = lapply(1:n_init, function(i){
#
#     if(verbose) cat(paste0("init iter: ", i, "\r"))
#
#     w = sample(1:M, size = model_data$dims$n_id, replace = TRUE)
#
#     run = find_map(
#       M = M,
#       X = kronecker(create_dummy(w, M), model_data$spline$B),
#       lambda = lambda_init,
#       stan_model = readRDS(system.file(
#         "stan", "stan_for_opt.rds", package = "bmixscat"
#       )),
#       model_data = model_data
#     )
#
#     run$logpost = run$opt_fit$par$log_posterior
#     run$w = w
#     run$beta = run$opt_fit$par$beta
#     return(run)
#
#   })
#
#   logpost = purrr::map_dbl(init_runs, ~{.x$logpost})
#   max_logpost = which.max(logpost)
#   w0 = init_runs[[max_logpost]]$w
#   names(w0) = model_data$data$id_unique
#   beta0 = cbind(init_runs[[max_logpost]]$beta, 0)
#   pw0 = init_runs[[max_logpost]]$w |> table() |> prop.table() |> as.numeric()
#
#   out = list(
#     w0 = w0,
#     pw0 = pw0,
#     beta0 = beta0
#   )
#
#   return(out)
#
# }


find_init_w = function(M,
                       model_data,
                       lambda_init = 1,
                       n_init = 30,
                       init_iters = 30,
                       init_burn_in = 15,
                       init_thin = 2,
                       verbose = TRUE,
                       seed) {

  set.seed(seed)

  init_runs = lapply(1:n_init, function(i){

    if(verbose) cat(paste0("init iter: ", i, "\r"))

    run = single_run(
      M = M,
      model_data = model_data,
      lambda = lambda_init,
      init_list = NULL,
      iters = init_iters,
      burn_in = init_burn_in,
      thin = init_thin,
      seed = NULL,
      verbose = FALSE
    )

    w = run$sample_list$w |> comp_class()

    opt_w = find_map(
      M = M,
      X = kronecker(create_dummy(w, M), model_data$spline$B),
      lambda = lambda_init,
      stan_model = readRDS(system.file(
        "stan", "stan_for_opt.rds", package = "mixscat"
      )),
      model_data = model_data
    )

    run$logpost = opt_w$opt_fit$par$log_posterior
    run$beta = cbind(opt_w$opt_fit$par$beta, 0)
    run$w = w

    return(run)

  })

  logpost = purrr::map_dbl(init_runs, ~{.x$logpost})
  max_logpost = which.max(logpost)
  w0 = init_runs[[max_logpost]]$w
  beta0 = init_runs[[max_logpost]]$beta
  pw0 = w0 |> factor(levels = 1:M) |> table() |> prop.table()
  pw0[pw0 == 0] = 1e-300

  out = list(
    w0 = w0,
    beta0 = beta0,
    pw0 = pw0
  )

  return(out)

}











