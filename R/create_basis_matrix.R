create_basis_matrix = function(n_basis, time) {

  B = splines::bs(
    x = time, df = max(c(4, n_basis)), degree = 3, intercept = TRUE
  )

  # QR decomposition of B
  qrB = qr(B)
  R = qr.R(qrB)
  Rinv = solve(R)
  Q = qr.Q(qrB)

  D = diff(diag(max(c(n_basis, 4))), differences = 1)
  S = crossprod(D)

  P = t(Rinv) %*% S %*% Rinv

  # Eigen decomposition
  eP = eigen(P)
  evals = eP$values
  evecs = eP$vectors

  # Threshold for numerical zero
  tol = 1e-10
  idx_null = which(evals < tol)
  idx_pos  = which(evals >= tol)

  # Fixed and random components
  U0 = evecs[, idx_null, drop = FALSE]
  Up = evecs[, idx_pos, drop = FALSE]
  lam_pos = evals[idx_pos]

  X = Q %*% U0
  Z = Q %*% Up

  var_z = apply(Z, MARGIN = 2, FUN = var)

  model_matrix = cbind(1, scale(Z[, rev(1:ncol(Z))], center = F, scale = T))
  model_matrix = model_matrix[, 1:n_basis]

  rownames(model_matrix) = as.character(time)

  out = list(model_matrix = model_matrix, nD = diag(c(0, rev(evals[idx_pos])/var_z)))

  return(out)
}
