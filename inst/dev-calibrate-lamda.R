library(tidyverse)
devtools::load_all()

ari = mclust::adjustedRandIndex
compute_asw = function(w, z_dist) {
  asw = cluster::silhouette(w, dist = z_dist)
  mean(asw[, 3])
}

#### data ####
data = readRDS("~/Documents/GitHub/clustering-mortality-two-step/data/input-data-cluster.rds")
sex_str = "Female"
G = 3

x = data %>%
  filter(sex == sex_str) %>%
  select(-country, -year, -sex, -(Z2:Z6)) %>%
  as.matrix()

id = data %>% select(country, year) %>% distinct() %>% .$country
time = data %>% select(country, year) %>% distinct() %>% .$year

z = data %>%
  filter(sex == sex_str) %>%
  as.data.frame() %>%
  .[, paste0("Z", G)]

levels = z %>%
  unique() %>%
  sort() %>%
  as.character()

z = factor(z, levels = levels)
n_id = length(unique(id))
n_time = length(unique(time))
id_unique = id |> unique() |> sort()

#### data model ####
M = 2
n_basis = 10
iters = 100
burn_in = 50
thin = 2
chains = 2
seed = 1


n_cores = 1
verbose = TRUE

lambda = 1
intercept_penalty = 1
dirichlet_param = 0.001
submodels = 1:M

model_data = create_model_data(
  z = z,
  id = id,
  time = time,
  n_basis = n_basis,
  intercept_penalty = intercept_penalty,
  dirichlet_param = dirichlet_param
)


opt_lambda = calibrate_lambda(
  w = get_w_ward(M = 1, model_data$data$z_dist),
  M = 10,
  intercept_penalty = 1,
  model_data = model_data,
  lambda_grid = seq(1, 120, 5)
)

opt_lambda$lambda_plot
opt_lambda$best_lambda

opt_lambda$log_posterior$log_like %>% plot()
opt_lambda$log_posterior$log_prior %>% plot()

