update_chain = function(beta,
                        w,
                        pw,
                        model_data,
                        dirichlet_param,
                        beta_precision_matrix,
                        update_w_iter = TRUE) {

  beta = update_beta(
    beta = beta,
    w = w,
    beta_precision_matrix = beta_precision_matrix,
    model_data = model_data
  )

  pw = update_pw(
    w = w,
    dirichlet_param = dirichlet_param,
    model_data = model_data
  )

  if(update_w_iter == TRUE) {

    w_out = update_w(
      beta = beta,
      pw = pw,
      model_data = model_data
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


