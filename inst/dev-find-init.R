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

data %>%
  filter(sex == "Male", country %in% c("Argentina","Bulgaria","Malaysia","Serbia","Uruguay")) %>%
  ggplot(aes(x=year, y=country, fill=factor(Z5))) +
  geom_tile(color="black")

M = 2
n_basis = 10
iters = 100
burn_in = 50
thin = 2
chains = 2
seed = 2

lambda = NULL
intercept_penalty = 1
dirichlet_param = 1 #1/M

init_control = list(
  lambda_init = 1,
  n_init = 10,
  lambda_grid = seq(from = 0.01, to = 5, length.out = 30),
  init_iters = 100,
  init_burn_in = 0,
  init_thin = 1,
  init_final_run = 50,
  verbose = TRUE
)

verbose = TRUE
n_cores = 1

model_data = create_model_data(
  z = z,
  id = id,
  time = time,
  n_basis = n_basis,
  intercept_penalty = intercept_penalty,
  dirichlet_param = dirichlet_param
)

init1 = find_init_w(
  M = M,
  model_data = model_data,
  seed = 1,
  init_control = init_control
)

init2 = find_init_w(
  M = M,
  model_data = model_data,
  seed = 2,
  init_control = init_control
)

init1$logpost
init2$logpost

w1 = init1$w
w2 = init2$w

ari(w1, w2)
w1 %>% sort()
w2 %>% sort()

table(w1, w2)

devtools::load_all()
init_runs = lapply(1, function(seed) {
  find_init_w(
    M = M,
    model_data = model_data,
    seed = seed,
    init_control = init_control
  )
})


init_runs[[1]]$logpost
init_runs[[1]]$w %>% unique() %>% length()
init_runs[[1]]$w %>% sort()

init_runs[[1]]$sample_list$w %>% comp_class() %>% sort()


init_runs[[1]]$logpost
init_runs[[2]]$logpost
init_runs[[3]]$logpost

init_runs[[1]]$w %>% unique() %>% length()
init_runs[[2]]$w %>% unique() %>% length()
init_runs[[3]]$w %>% unique() %>% length()

mclust::adjustedRandIndex(init_runs[[1]]$w, init_runs[[2]]$w)
mclust::adjustedRandIndex(init_runs[[1]]$w, init_runs[[3]]$w)
mclust::adjustedRandIndex(init_runs[[2]]$w, init_runs[[3]]$w)

init_runs[[1]]$w %>% sort()


run1 = single_run(
  M = M,
  w = NULL,
  model_data = model_data,
  lambda = init_control$lambda_init,
  init_list = NULL,
  iters = init_control$init_iters,
  burn_in = init_control$init_burn_in,
  thin = init_control$init_thin,
  verbose = TRUE,
  seed = 1
)

run2 = single_run(
  M = M,
  w = NULL,
  model_data = model_data,
  lambda = init_control$lambda_init,
  init_list = NULL,
  iters = init_control$init_iters,
  burn_in = init_control$init_burn_in,
  thin = init_control$init_thin,
  verbose = TRUE,
  seed = 4
)

w1 = run1$sample_list$w[nrow(run1$sample_list$w), ]
w2 = run2$sample_list$w[nrow(run2$sample_list$w), ]
mclust::adjustedRandIndex(w1, w2)










