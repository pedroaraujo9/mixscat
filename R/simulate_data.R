simulate_data = function(seed) {

  set.seed(seed)

  #### define dimensions and hyperparameters ####
  M = 3
  G = 3
  n_id = 30
  n_time = 60
  n = n_id * n_time
  id_unique = paste0("id", 1:n_id)
  id = rep(id_unique, each = n_time)
  time = rep(1:n_time, times = n_id)

  #### sample global group ####
  w = rep(1:M, each = n_id / M)
  names(w) = id_unique

  #### sample local group ####
  ZM = matrix(nrow = n_id, ncol = n_time)

  for(i in 1:n_id) {

    if(w[i] == 1) {
      # 1 -- > 2
      thr = sample(1:10, size = 2)
      ZM[i,] = c(
        rep(1, 30 - thr[1]),
        rep(2, 60 - 30 + thr[1])
      )

    }else if(w[i] == 2) {
      # 1 -- > 3
      thr = sample(1:10, size = 2)
      ZM[i,] = c(
        rep(1, 30 - thr[1]),
        rep(3, 60 - 30 + thr[1])
      )

    }else if(w[i] == 3) {
      # 1 --> 2 --> 2
      thr = sample(1:5, size = 2)
      ZM[i,] = c(
        rep(1, 20 - thr[1]),
        rep(2, 20 - thr[2]),
        rep(3, 60 - 40 + thr[1] + thr[2])
      )
    }
  }

  z = as.vector(t(ZM))

  #### return ####
  data_sim = list(
    z = z,
    id = id,
    time = time,
    w = w
  )

  return(data_sim)

}



