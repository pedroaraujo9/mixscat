augment_beta = function(g,
                        Z,
                        X,
                        linear_pred,
                        precision_matrix) {

  max_other = apply(linear_pred[, -g, drop = FALSE], 1, max)
  C = max_other + log(rowSums(exp(linear_pred[, -g, drop = FALSE] - max_other)))
  n = nrow(Z)

  omega = BayesLogit::rpg(n, h = 1, z = linear_pred[, g] - C)

  beta = sample_beta(
    X = X,
    omega = omega,
    inv_cov = precision_matrix,
    z = Z[, g],
    C = C
  )



  return(as.numeric(beta))
}

update_beta = function(beta,
                       w,
                       lambda,
                       fixed_lambda = FALSE,
                       beta_precision_matrix = NULL,
                       model_data) {

  G = model_data$dims$G
  M = model_data$dims$M
  B = model_data$spline$B
  intercept_penalty = model_data$spline$intercept_penalty

  beta = rbind(beta)

  W = create_dummy(w, M)
  X = kronecker(W, B)
  X[, 1] = 1
  Z = model_data$data$Z

  if(fixed_lambda == FALSE) {

    beta_precision_matrix = create_beta_precision_matrix(
      model_data = model_data,
      lambda = lambda,
      intercept_penalty = intercept_penalty,
      M = M
    )

  }

  for(g in 1:(G-1)) {

    linear_pred = X %*% beta

    beta[, g] = augment_beta(
      g = g,
      Z = Z,
      X = X,
      linear_pred = linear_pred,
      precision_matrix = beta_precision_matrix
    )
  }

  beta[, G] = 0

  return(beta)

}






