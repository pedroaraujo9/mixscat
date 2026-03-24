update_pw = function(w, dirichlet_param = 1, model_data) {

  M = model_data$dims$M

  if(M > 1) {

    w = factor(w, levels = 1:M, ordered = T)
    nw = w |> table()
    pw = extraDistr::rdirichlet(1, alpha = as.numeric(nw + dirichlet_param))
    pw[1, pw[1, ] < 1e-300] = 1e-300

  }else{

    pw = matrix(1)

  }

  return(pw)
}
