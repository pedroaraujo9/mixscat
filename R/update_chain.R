update_chain = function(beta,
                        w,
                        pw,
                        lambda,
                        model_data,
                        dirichlet_param,
                        a_lambda = 1,
                        b_lambda = 1,
                        beta_precision_matrix = NULL,
                        fixed_lambda = FALSE,
                        fixed_w = FALSE,
                        fixed_pw = FALSE) {

  beta = update_beta(
    beta = beta,
    w = w,
    lambda = lambda,
    beta_precision_matrix = beta_precision_matrix,
    model_data = model_data,
    fixed_lambda = fixed_lambda
  )

  if(fixed_lambda == FALSE) {

    lambda = update_lambda(
      beta = beta,
      a_lambda = a_lambda,
      b_lambda = b_lambda,
      model_data = model_data
    )

  }

  if(fixed_pw == FALSE) {

    pw = update_pw(
      w = w,
      dirichlet_param = dirichlet_param,
      model_data = model_data
    )

  }


  if(fixed_w == FALSE) {

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
    pw = pw,
    lambda = lambda
  )

  return(out)

}


