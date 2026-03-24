update_lambda = function(beta, a_lambda = 1, b_lambda = 1, model_data) {

  M = model_data$dims$M
  G = model_data$dims$G
  S = model_data$spline$S
  n_basis = model_data$spline$n_basis
  s_matrix = model_data$temp$s_matrix[, -G]
  idx_list = model_data$temp$basis_idx_list
  lambda = numeric(M)

  lambda = update_lambda_loop_cpp(
    beta = beta[, -G],
    idx_list = idx_list,
    s_matrix = s_matrix[-1, ],
    a_lambda = a_lambda,
    b_lambda = b_lambda,
    n_basis = n_basis,
    G = G
  ) |> as.numeric()

  return(lambda)

}
