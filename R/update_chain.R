update_chain = function(beta,
                        w,
                        pw,
                        model_data,
                        update_w_iter = TRUE,
                        S,
                        temperature = 1,
                        w_temp = 1) {

  #model_data$spline$B = model_data$spline$B * sqrt(1/temperature)

  beta = update_beta(
    beta = beta,
    w = w,
    model_data = model_data,
    S = S
  )

  pw = update_pw(
    w = w,
    dirichlet_param = model_data$spline$dirichlet_param,
    model_data = model_data
  )

  if(update_w_iter == TRUE) {

    w_out = update_w(
      beta = beta,
      pw = pw,
      model_data = model_data,
      w_temp = w_temp
    )

    w = w_out$w
    w_post_prob = w_out$w_post_prob

  } else {

    w_post_prob = create_dummy(w, G = length(pw))

  }

  #model_data$spline$B = model_data$spline$B * sqrt(temperature)

  out = list(
    beta = beta,
    w = w,
    w_post_prob = w_post_prob,
    pw = pw
  )

  return(out)

}


