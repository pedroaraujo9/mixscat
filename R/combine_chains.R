build_block_perm_mat = function(perm, n_basis) {

  M = length(perm)
  P_sigma = diag(M)[perm, ]
  T_mat = kronecker(P_sigma, diag(n_basis))

  return(T_mat)
}

combine_chains = function(runs, model_data) {

  sample_list = runs[[1]]$sample_list
  M = ncol(sample_list$pw)

  if(length(runs) > 1) {

    for(i in 2:length(runs)) {

      for(param in c("beta", "pw", "w", "w_post_prob", "logpost", "spline_probs")) {

        sample_list[[param]] = abind::abind(
          sample_list[[param]], runs[[i]]$sample_list[[param]], along = 1
        )

      }

    }

  }

  if(M > 1) {

    ls = label.switching::label.switching(
      method = "STEPHENS",
      z = sample_list$w,
      K = M,
      p = sample_list$w_post_prob,
    )

    perm = ls$permutations$STEPHENS
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

