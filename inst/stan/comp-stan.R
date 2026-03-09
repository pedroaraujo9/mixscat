library(rstan)

model = rstan::stan_model(
  file = "inst/stan/stan_for_opt.stan",
  model_name = "mixscat_stan",
  save_dso = T,
  auto_write = T
)

