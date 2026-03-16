apply_relabel = function(sample_list, pivot) {

  M = length(sample_list$pw[1, ])

  if(M > 1) {

    ls = label.switching::label.switching(
      method = "ECR",
      z = sample_list$w,
      K = M,
      p = sample_list$w_post_prob, zpivot = pivot
    )

    perm = ls$permutations$ECR
    n_basis = model_data$spline$n_basis

    for(i in 1:nrow(perm)) {

      if(any(perm[i, ] != (1:M))) {

        sample_list$pw[i, ] = sample_list$pw[i, perm[i, ]]
        sample_list$w[i, ] = perm[i, sample_list$w[i, ]]
        sample_list$beta[i,,] = build_block_perm_mat(perm[i, ], n_basis) %*% sample_list$beta[i,,]
        sample_list$w_post_prob[i,,] = sample_list$w_post_prob[i,,perm[i, ]]
      }

    }
  }

  return(sample_list)

}


relabel = function(runs) {

  logpost = purrr::map_dbl(runs, ~{
    mean(.x$sample_list$logpost[, "logpost"])
  })

  pivot = runs[[which.max(logpost)]]$sample_list$w |> comp_class()

  for(chain in 1:length(runs)) {
    runs[[chain]]$sample_list = apply_relabel(
      runs[[chain]]$sample_list, pivot = pivot
    )
  }

  return(runs)

}
