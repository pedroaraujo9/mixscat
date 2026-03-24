build_block_perm_mat = function(perm, n_basis) {

  M = length(perm)
  P_sigma = diag(M)[perm, ]
  T_mat = kronecker(P_sigma, diag(n_basis))

  return(T_mat)
}

apply_relabel = function(sample_list, pivot) {

  iters = nrow(sample_list$w_post_prob)
  M = dim(sample_list$w_post_prob)[3]

  ls = label.switching::label.switching(
    method = "ECR",
    z = sample_list$w,
    zpivot = pivot,
    K = M
  )


  for(i in 1:iters) {
    perm = labels = ls$permutations$`ECR`[i, ]
    sample_list$w[i, ] = order(perm)[sample_list$w[i, ]]
    if(!is.null(sample_list$pw)) sample_list$pw[i, ] = sample_list$pw[i, perm]
    sample_list$beta[i,,] = build_block_perm_mat(perm, n_basis) %*% sample_list$beta[i,,]
    sample_list$w_post_prob[i,,] = sample_list$w_post_prob[i,, perm]
  }

  sample_list = sample_list

  return(sample_list)

}

relabel = function(runs, pivot = NULL) {

  if(is.null(pivot)) {

    pivot = runs[[1]]$sample_list$w[nrow(runs[[1]]$sample_list$w), ]

  }

  for(chain in 1:length(runs)) {
    runs[[chain]]$sample_list = apply_relabel(
      runs[[chain]]$sample_list, pivot = pivot
    )
  }

  return(runs)

}
