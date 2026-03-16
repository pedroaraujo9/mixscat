update_chain = function(beta,
                        w,
                        pw,
                        model_data,
                        update_w_iter = TRUE,
                        temperature = 1) {

  beta = update_beta(
    beta = beta,
    w = w,
    model_data = model_data
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
      temperature = temperature
    )

    w = w_out$w
    w_post_prob = w_out$w_post_prob

  } else {

    w_post_prob = create_dummy(w, G = length(pw))

  }


  out = list(
    beta = beta,
    w = w,
    w_post_prob = w_post_prob,
    pw = pw
  )

  return(out)

}


