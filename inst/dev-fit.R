library(tidyverse)
devtools::load_all()

ari = mclust::adjustedRandIndex
compute_asw = function(w, z_dist) {
  asw = cluster::silhouette(w, dist = z_dist)
  mean(asw[, 3])
}

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
n_id = length(unique(id))
n_time = length(unique(time))
id_unique = id |> unique() |> sort()

z_seq = z |>
  matrix(nrow = n_id, ncol = n_time, byrow = T) |>
  TraMineR::seqdef()

rownames(z_seq) = id_unique
z_dist = TraMineR::seqdist(z_seq, method = "HAM")


data.frame(
  id = id,
  time = time,
  z = z
) %>%
  ggplot(aes(x=time, y=id, fill=z)) +
  geom_tile(color="black") +
  viridis::scale_fill_viridis(discrete = T)

#### fit ####
chains = 2
iters = 400
burn_in = 200
thin = 5
lambda = 1
dirichlet_param = 1
intercept_penalty = 1
seed = 1
n_init = 5
init_mcmc_iters = 100
verbose = TRUE
n_basis = 10
n_cores = 1

M = c(19)

devtools::load_all()
fit = fit_mixscat(
  M = M,
  z = z,
  id = id,
  time = time,
  n_basis = n_basis,
  iters = iters,
  thin = thin,
  burn_in = burn_in,
  chains = chains,
  lambda = lambda,
  intercept_penalty = intercept_penalty,
  dirichlet_param = dirichlet_param,
  n_init = n_init,
  init_mcmc_iters = init_mcmc_iters,
  n_cores = n_cores,
  seed = seed,
  verbose = TRUE
)

2 * log(71) * 82.34755
fit$cluster_metrics
fit$cluster_metrics$penal_BICM %>% plot(x = M, y = .)
fit$cluster_metrics$loglike_complex
fit$cluster_metrics$ASW %>% plot()
fit$cluster_metrics$ASW %>% which.max()
fit$cluster_metrics$BICM %>% plot()

fit$fit$`M=18`$conv_check$logpost
fit$fit$`M=19`$sample_list$logpost[, "loglike"] %>% var()
fit$cluster_metrics$logpenal_complex %>% plot()
fit$cluster_metrics$avg_loglike %>%  plot()


