check_conv_logpost = function(runs) {

  logpost_chains = lapply(runs, function(run){
    run$sample_list$logpost[, "logpost"]
  })

  logpost_chains = do.call(cbind, logpost_chains)
  rhat = posterior::rhat(logpost_chains)
  ess = posterior::ess_basic(logpost_chains)

  out = data.frame(rhat = rhat, ess = ess)
  return(out)

}

check_conv_param = function(runs, param_name) {

  dimension = dim(runs[[1]]$sample_list[[param_name]])[-1]

  conv_metrics = lapply(1:dimension[1], function(i){

    lapply(1:dimension[2], function(j){

      param_chains = lapply(runs, function(run){
        run$sample_list[[param_name]][, i, j]
      })

      param_chains = do.call(cbind, param_chains)
      rhat = posterior::rhat(param_chains)
      ess = posterior::ess_basic(param_chains)

      out = data.frame(rhat = rhat, ess = ess, i = i, j = j)
      return(out)

    }) %>% purrr::list_rbind()
  }) %>% purrr::list_rbind()

  rhat = conv_metrics |>
    dplyr::select(-ess) |>
    tidyr::spread(j, rhat) |>
    dplyr::select(-i)

  ess = conv_metrics |>
    dplyr::select(-rhat) |>
    tidyr::spread(j, ess) |>
    dplyr::select(-i)

  out = list(
    rhat = rhat,
    ess = ess
  )

  return(out)

}

