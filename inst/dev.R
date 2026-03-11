library(tidyverse)
devtools::load_all()
#### data ####
data = readRDS("~/Documents/GitHub/clustering-mortality-two-step/data/input-data-cluster.rds")

x_male = data %>%
  filter(sex == "Male") %>%
  select(-country, -year, -sex, -(Z2:Z6)) %>%
  as.matrix()

id = data %>% select(country, year) %>% distinct() %>% .$country
time = data %>% select(country, year) %>% distinct() %>% .$year

z_male = data %>%
  filter(sex == "Male") %>%
  as.data.frame() %>%
  .[, paste0("Z", 5)]

male_levels = z_male %>%
  unique() %>%
  sort() %>%
  as.character()

z = factor(z_male, levels = male_levels)

M = 1:20
n_basis = 10
iters = 100
burn_in = 50
thin = 2
chains = 2
seed = 1

lambda = NULL
intercept_penalty = 1
w_dirichlet = 1

init_control = list(
  lambda_init = 1,
  n_init = 5,
  lambda_grid = seq(from = 0.01, to = 5, length.out = 30),
  init_iters = 10,
  init_burn_in = 5,
  init_thin = 2,
  init_final_run = 100,
  verbose = FALSE
)


verbose = TRUE
n_cores = 1

devtools::load_all()
fit = fit_mixscat(
  M = M,
  z = z,
  id = id,
  time = time,
  n_basis = n_basis,
  iters = iters,
  burn_in = burn_in,
  thin = thin,
  chains = chains,
  seed = seed,
  lambda = NULL,
  intercept_penalty = intercept_penalty,
  w_dirichlet = w_dirichlet,
  verbose = verbose,
  init_control = init_control,
  n_cores = n_cores
)

fit$cluster_metrics
fit$clusters[, 1] %>% unique()
fit$best_lambda %>% do.call(c, .)
fit$cluster_metrics$ASW %>% plot()
devtools::load_all()
model_data = create_model_data(
  z = z, id = id, time = time
)

chains = 2
iters = 500
burn_in = 200
thin = 5
lambda_init = 1
lambda_limits = c(0.01, 5)
lambda_grid_n_points = 30
n_init = 50
init_iters = 10
init_burn_in = 5
init_thin = 1

devtools::load_all()
fit = lapply(1:2, function(seed) {
  pipeline(
    model_data = model_data,
    M = M,
    chains = chains,
    iters = iters,
    thin = thin,
    burn_in = burn_in,
    lambda_init = lambda_init,
    lambda_limits = lambda_limits,
    lambda_grid_n_points = lambda_grid_n_points,
    n_init = n_init,
    init_iters = init_iters,
    init_burn_in = init_burn_in,
    init_thin = init_thin,
    seed = seed
  )
})


fit[[1]]$conv_check
fit[[2]]$conv_check
fit[[1]]$metrics$BICM
fit[[2]]$metrics$BICM

mclust::adjustedRandIndex(fit[[1]]$w, fit[[2]]$w)

table(fit[[1]]$w, fit[[2]]$w)
fit[[2]]$w %>% sort()
