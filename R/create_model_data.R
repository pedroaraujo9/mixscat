create_model_data = function(z,
                             id,
                             time,
                             n_basis = 10,
                             intercept_penalty = 1,
                             dirichlet_param = 1) {

  z_levels = sort(unique(z))
  G = length(z_levels)

  data = data.frame(
    z = z,
    id = id,
    time = time
  ) |>
    mutate(id = factor(id, levels = unique(id))) |>
    dplyr::arrange(id, time)

  data$z_label = factor(z, levels = z_levels)
  data$z = as.numeric(data$z_label)

  z = data$z
  Z = create_dummy(z, G)

  id = data$id
  time = data$time
  n = length(z)

  id_unique = unique(id)
  time_unique = unique(time)

  n_id = length(id_unique)
  n_time = length(time_unique)

  Z_matrix = matrix(z, nrow = n_id, ncol = n_time, byrow = TRUE)
  z_dist = TraMineR::seqdef(Z_matrix) |> TraMineR::seqdist(method = "DHD")
  colnames(z_dist) = rownames(z_dist) =  id_unique

  time_seq = time - min(time) + 1
  time_seq_unique = unique(time_seq)

  id_time_df = tibble::tibble(
    id = id,
    time_seq = time_seq,
    time = time
  ) |>
    dplyr::mutate(id = factor(id, levels = id_unique))

  basis_funcions = create_basis_matrix(
    n_basis = n_basis,
    time = time_seq_unique
  )

  B = basis_funcions$model_matrix
  S = basis_funcions$nD
  S[1, 1] = intercept_penalty

  out = list()

  out$data = list(
    z = z,
    z_label = data$z_label,
    Z = Z,
    z_levels = z_levels,
    z_dist = z_dist,
    id = id,
    time = time,
    time_seq = time_seq,
    time_unique = time_unique,
    id_unique = id_unique
  )

  out$dims = list(
    n = n,
    n_id = n_id,
    n_time = n_time,
    G = G
  )

  out$spline = list(
    dirichlet_param = dirichlet_param,
    n_basis = n_basis,
    intercept_penalty = intercept_penalty,
    B = B,
    S = S
  )

  return(out)

}
