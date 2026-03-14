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
  n_init = 100,
  lambda_grid = seq(from = 0.01, to = 5, length.out = 30),
  init_iters = 20,
  init_burn_in = 1,
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


ws = lapply(1:200, function(i){

  cat(i, "\r")

  run1 = single_run(
    M = 10,
    w = NULL,
    model_data = model_data,
    lambda = init_control$lambda_init,
    init_list = NULL,
    iters = 10, #init_control$init_iters,
    burn_in = 0, #init_control$init_burn_in,
    thin = 1, #init_control$init_thin,
    verbose = FALSE,
    seed = NULL
  )

  w = run1$sample_list$w[nrow(run1$sample_list$w), ]
  w

}) %>% do.call(rbind, .)

devtools::load_all()



M = 20
n_basis = 10
iters = 100
burn_in = 50
thin = 2
chains = 2
seed = 2

lambda = NULL
intercept_penalty = 1
dirichlet_param = 0.01

init_control = list(
  lambda_init = 1,
  n_init = 100,
  lambda_grid = seq(from = 0.01, to = 5, length.out = 30),
  init_iters = 20,
  init_burn_in = 1,
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

iters = 500
cut_point = 200
temp_seq = seq(from = 1000, to = 1, length.out = cut_point)
temp_seq = c(temp_seq, rep(1, times = iters - cut_point))
plot(temp_seq)
temp_seq = rep(1000, iters)

devtools::load_all()
run1 = single_run(
  M = M,
  w = NULL,
  model_data = model_data,
  lambda = 1,
  init_list = NULL,
  iters = iters, #init_control$init_iters,
  burn_in = 0, #init_control$init_burn_in,
  thin = 1, #init_control$init_thin,
  verbose = TRUE,
  seed = NULL,
  temperature = temp_seq
)


Md = lapply(1:nrow(run1$sample_list$w), function(i){
  run1$sample_list$w[i, ] %>% unique() %>% length()
}) %>% do.call(c, .)

Md %>% tail()

plot(Md, type = "l")

run1$sample_list$w[nrow(run1$sample_list$w), ] %>% sort()
unique(run1$sample_list$w[nrow(run1$sample_list$w), ])
run1$sample_list$w[, "Argentina"]

w2 = run1$sample_list$w[nrow(run1$sample_list$w), ]
ari(w1, w2)


w1 = run1$sample_list$w %>% comp_class()
w2 = run1$sample_list$w %>% comp_class()
ari(w1, w2)

unique(w1)
unique(w2)

run1$sample_list$w[, 1]

run1$S_expand %>% diag() %>% round(5)
run1$sample_list$beta[10,,] %>% round(5)
run1$sample_list$w[, 10]
run1$sample_list$w %>% comp_class() %>% sort() %>% unique()

ws

ari = mclust::adjustedRandIndex
ari_matrix = lapply(1:nrow(ws), function(j){
  lapply(1:nrow(ws), function(i){
    ari(ws[i, ], ws[j, ])
  }) %>% do.call(c, .)
}) %>% do.call(rbind, .)

ari_matrix[lower.tri(ari_matrix)] %>% summary()
sum(ari_matrix[lower.tri(ari_matrix)] == 1)

i = 20
im = ari_matrix[i, ][-i] %>% which.max()
ari_matrix[i, im]

w1 = ws[i, ]
w2 = ws[im, ]

table(w1, w2)








