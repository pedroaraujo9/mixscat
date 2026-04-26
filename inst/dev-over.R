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

data.frame(
  id = id,
  time = time,
  z = z
) %>%
  ggplot(aes(x=time, y=id, fill=z)) +
  geom_tile(color="black") +
  viridis::scale_fill_viridis(discrete = T)

M_max = 15
iters = 4000
burn_in = 2000
thin = 10
chains = 3
n_cores = 1

n_init = 200
init_iters = 20

n_basis = 10
intercept_penalty = 1
dirichlet_param = 0.01
a_lambda = 1
b_lambda = 1
seed = 1

w = NULL
lambda = 1

devtools::load_all()
clust_number = find_number_clust(
  M_max = M_max,
  z = z,
  id = id,
  time = time,
  n_basis = n_basis,
  n_init = n_init,
  init_iters = init_iters,
  lambda = lambda,
  dirichlet_param = dirichlet_param,
  a_lambda = a_lambda,
  b_lambda = b_lambda,
  intercept_penalty = intercept_penalty,
  seed = seed,
  verbose = TRUE
)

clust_number$clust_size_p

M_best = 10
pw = clust_number$post_modes[[as.character(M_best)]]$w_post_prob
init_list = list(
  w = clust_number$post_modes[[as.character(M_best)]]$w,
  beta = clust_number$post_modes[[as.character(M_best)]]$beta_map
)

devtools::load_all()
fit = fit_mixscat(
  M = M_best,
  z = z,
  id = id,
  time = time,
  pw = pw,
  n_basis = n_basis,
  init_list = init_list,
  iters = iters,
  thin = thin,
  burn_in = burn_in,
  chains = chains,
  n_cores = n_cores,
  lambda = lambda,
  w = NULL,
  intercept_penalty = intercept_penalty,
  dirichlet_param = dirichlet_param,
  a_lambda = a_lambda,
  b_lambda = b_lambda,
  seed = seed,
  verbose = TRUE
)

w1 = fit$fit$runs[[1]]$sample_list$w %>% comp_class()
w2 = fit$fit$runs[[2]]$sample_list$w %>% comp_class()
w3 = fit$fit$runs[[3]]$sample_list$w %>% comp_class()

w1 %>% unique() %>% sort()
w2 %>% unique() %>% sort()
w3 %>% unique() %>% sort()

ari(w1, w2)
ari(w1, w3)
ari(w2, w3)

which(w1 != w2)
fit$fit$runs[[1]]$sample_list$w[, "Brazil"]
fit$fit$runs[[2]]$sample_list$w[, "Brazil"]
fit$fit$runs[[3]]$sample_list$w[, "Brazil"]


