update_pw = function(w, dirichlet_param = 1, model_data) {

  M = model_data$dims$M

  if(M > 1) {

    w = factor(w, levels = 1:M, ordered = T)
    nw = w |> table()
    pw = extraDistr::rdirichlet(1, alpha = as.numeric(nw + dirichlet_param)) |> as.numeric()
    pw[pw < 1e-200] = 1e-200

  }else{

    pw = 1

  }

  return(pw)
}
