# sample_beta_r = function(X,
#                          omega,
#                          inv_cov,
#                          z,
#                          C) {
#
#   X_omega = X * matrix(omega, nrow = nrow(X), ncol = ncol(X), byrow = F)
#   Xt_omega_X = crossprod(X_omega, X)
#   V = solve(Xt_omega_X + inv_cov)
#
#   m = V %*% (crossprod(X, cbind(z - 0.5 + omega * C)))
#   beta = m + t(chol(V)) %*% rnorm(ncol(X))
#   return(beta)
#
# }

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
                       beta_precision_matrix,
                       model_data) {

  G = model_data$dims$G
  M = model_data$dims$M
  B = model_data$spline$B

  beta = rbind(beta)

  W = create_dummy(w, M)
  X = kronecker(W, B)
  Z = model_data$data$Z

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






