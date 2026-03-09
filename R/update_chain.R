update_chain = function(beta,
                        w,
                        pw,
                        model_data) {

  beta = update_beta(
    beta = beta,
    w = w,
    model_data = model_data
  )

  pw = update_pw(
    w = w,
    model_data = model_data
  )

  w_out = update_w(
    beta = beta,
    pw = pw,
    model_data = model_data
  )

  w = w_out$w
  w_post_prob = w_out$w_post_prob

  out = list(
    beta = beta,
    w = w,
    w_post_prob = w_post_prob,
    pw = pw
  )

  return(out)

}


