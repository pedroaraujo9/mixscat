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

  return(sample_list)

}

