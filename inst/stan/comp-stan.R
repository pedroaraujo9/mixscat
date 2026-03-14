library(rstan)

model = rstan::stan_model(
  file = "inst/stan/stan_for_opt.stan",
  model_name = "mixscat_stan",
  save_dso = T,
  auto_write = T
)

model2 = rstan::stan_model(
  file = "inst/stan/stan_full_logpost.stan",
  model_name = "mixscat_full_post",
  save_dso = T,
  auto_write = T
)

