augment_beta = function(g,
                        Z,
                        X,
                        linear_pred,
                        precision_matrix) {

  max_other = apply(linear_pred[, -g, drop = FALSE], 1, max)
  C = max_other + log(rowSums(exp(linear_pred[, -g, drop = FALSE] - max_other)))
  n = nrow(Z)

  omega = BayesLogit::rpg(n, h = 1, z = linear_pred[, g, drop = FALSE] - C)

  beta = sample_beta(
    X = X,
    omega = omega,
    inv_cov = precision_matrix,
    z = Z[, g, drop = FALSE],
    C = C
  )

  return(as.numeric(beta))
}

update_beta = function(beta,
                       w,
                       model_data) {

  G = model_data$dims$G
  M = model_data$dims$M

  B = model_data$spline$B
  S = model_data$spline$S_expand

  beta = rbind(beta)

  W = create_dummy(w, M)
  X = kronecker(W, B)
  X = cbind(1, X)
  Z = model_data$data$Z

  for(g in 1:(G-1)) {

    linear_pred = X %*% beta

    beta[, g] = augment_beta(
      g = g,
      Z = Z,
      X = X,
      linear_pred = linear_pred,
      precision_matrix = S
    )
  }

  beta[, G] = 0

  return(beta)

}






